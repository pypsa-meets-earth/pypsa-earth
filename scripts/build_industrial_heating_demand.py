# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import os
import random
from copy import deepcopy

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _helpers import mock_snakemake
from _industry_heat_helpers import coords_to_relative_utm
from build_egs_potentials import tif_to_gdf
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
from scipy.spatial import distance, distance_matrix
from shapely.geometry import Point
from tqdm import tqdm


def prepare_demand_data(fn):

    df = gpd.read_file(fn)

    id_cols = [
        # 'names_primary', 'geometry',
        # 'industrial_type', 'industrial_sub_type'
        "FACILITY_ID",
        "FACILITY_NAME",
        "geometry",
    ]

    heating_cols = [
        col
        for col in df.columns
        if "Heating (BBtu)" in col and "total" not in col and "<250" not in col
    ]

    # Melt the dataframe to unpivot heating columns
    df_melted = df.melt(
        id_vars=id_cols,
        value_vars=heating_cols,
        var_name="Temperature Range",
        value_name="Heating Demand (BBtu)",
    )

    df_melted = df_melted[df_melted["Heating Demand (BBtu)"] != 0]

    df_melted["Temperature Range"] = df_melted["Temperature Range"].str.replace(
        "Heating (BBtu), ", "", regex=False
    )

    # df_melted['Heating Demand (MWh)'] = df_melted['Heating Demand (BBtu)'] * 293.07107
    df_melted["Heating Demand (BBtu)"] = (
        df_melted["Heating Demand (BBtu)"] * 1e9
    )  # returns new data to old scaling
    df_melted["Heating Demand (MWh)"] = df_melted["Heating Demand (BBtu)"] * 2.9307e-7

    df_melted = df_melted.drop(columns=["Heating Demand (BBtu)"])

    df_melted = df_melted[id_cols + ["Temperature Range", "Heating Demand (MWh)"]]

    cutoff = 5000  # MWh; cutoff of 5GWh removes around 1% of total demand, and 50% of datapoints

    df_melted = df_melted.loc[df_melted["Heating Demand (MWh)"] > cutoff]

    df_melted = df_melted.reset_index(drop=True)

    temp_mapper = {
        "0 to 49°C": 50,  # this dict implies choices about the temperature bands
        "50 to 79°C": 79,  # only demands up to 249C are relevant
        "80 to 149°C": 149,
        "150 to 199°C": 199,
        "200 to 249°C": 249,
        "250 to 299°C": 299,
        "300 to 349°C": 350,
        "350 to 399°C": 400,
        "400 to 449°C": 450,
        ">450°C": 500,
    }

    for old, new in temp_mapper.items():
        df_melted.replace(old, new, inplace=True)

    df_melted.rename(
        columns={
            "Temperature Range": "temperature",
            "Heating Demand (MWh)": "total_demand",
        },
        inplace=True,
    )
    df_melted["avg_demand"] = df_melted["total_demand"] / 8760

    df_melted = df_melted.loc[df_melted["temperature"] <= 250]

    df_melted["x"] = df_melted["geometry"].x
    df_melted["y"] = df_melted["geometry"].y

    return df_melted


def compute_mst_cost_and_diameter(cluster_points):
    """
    Given an array of points in 2D, compute:
      1) MST total edge length
      2) MST diameter (the longest path between any two points along MST edges)

    We'll use scipy.sparse.csgraph.minimum_spanning_tree for MST edges,
    then compute the "shortest path" distances between all pairs of points *within the MST*
    to find the MST diameter.

    Returns
    -------
    mst_length : float
        sum of MST edges
    mst_diameter : float
        the longest distance between any two points along the MST paths
    """
    if len(cluster_points) == 1:
        return 0.0, 0.0  # Single point => no edges, diameter=0

    # Pairwise distance matrix
    dist_matrix = distance_matrix(cluster_points, cluster_points)
    # Compute MST (returns a sparse matrix of MST edges)
    mst_csr = minimum_spanning_tree(dist_matrix)
    mst_length = mst_csr.sum()

    # Convert MST edges to an adjacency matrix so we can compute the MST diameter
    mst_adj = mst_csr.toarray()
    mst_adj = mst_adj + mst_adj.T  # symmetrical for undirected MST

    # shortest path distances within MST (Floyd-Warshall or BFS if MST is a tree)
    # For small clusters, Floyd-Warshall is fine, or use BFS from each node in MST.
    from scipy.sparse.csgraph import shortest_path

    dist_all_pairs = shortest_path(csgraph=mst_adj, directed=False)

    # MST diameter = max distance between any two vertices in MST
    mst_diameter = np.nanmax(dist_all_pairs)

    return mst_length, mst_diameter


def greedy_clustering(points, values, T, min_sum_value=10.0, max_sum_value=np.inf):
    """
    Given a diameter threshold T, build clusters greedily:
      1. Sort points by descending value
      2. For each unassigned point:
         - start a cluster and try to add other unassigned points if it doesn't break MST diameter <= T
         - stop once the sum of values >= min_sum_value or can't add more points

    Returns
    -------
    clusters : list of lists
        each cluster is a list of point indices
    cluster_mst_info : list of tuples
        parallel structure to clusters, each element is (mst_length, mst_diameter, sum_values)
    """

    n_points = len(points)
    assigned = [False] * n_points
    idx_sorted = np.argsort(values)[::-1]  # descending by value

    clusters = []
    cluster_mst_info = []

    for idx in idx_sorted:
        if assigned[idx]:
            continue

        # Start a cluster with just this point
        candidate_indices = [idx]
        current_sum = values[idx]
        assigned[idx] = True

        # Try adding more points to this cluster if it doesn't break the MST diameter T
        # We'll do a naive incremental approach:
        #   - Look for the nearest unassigned point that might still keep diameter <= T
        #   - If the MST diameter goes beyond T, revert adding that point.
        #   - Stop once sum >= min_sum_value or no points can be added without violation.
        improved = True
        # while (current_sum < min_sum_value) and improved:
        while (current_sum < max_sum_value) and improved:

            improved = False

            best_point = None

            best_length = T

            for j in range(n_points):
                if assigned[j]:
                    continue

                test_indices = candidate_indices + [j]
                cluster_pts = points[test_indices]

                mst_length, mst_diam = compute_mst_cost_and_diameter(cluster_pts)

                if mst_diam <= T:

                    if mst_length < best_length:

                        best_point = j
                        best_length = mst_length

            if best_point is not None:
                candidate_indices.append(best_point)
                assigned[best_point] = True
                current_sum += values[best_point]
                improved = True

        if current_sum < min_sum_value:
            # Not enough points in the cluster
            for pt in candidate_indices:
                assigned[pt] = False
            continue

        final_pts = points[candidate_indices]
        mst_length, mst_diam = compute_mst_cost_and_diameter(final_pts)
        clusters.append(candidate_indices)
        cluster_mst_info.append((mst_length, mst_diam, current_sum))

    return clusters, cluster_mst_info


def evaluate_solution(cluster_mst_info, all_values):
    """
    Given the MST info for each cluster: (mst_length, mst_diameter, sum_values)
    and all values of points that could have been in a cluster
    compute:
      - fraction_of_values_in_clusters
      - cluster_total_value_weighted_mst_cost

    Returns
    -------
    frac_in_clusters : float
        fraction of total point values that ended up in clusters
    weighted_mst_cost : float
        (sum of (mst_length_i * cluster_value_i)) / sum_of_cluster_values
    """

    df = pd.DataFrame(cluster_mst_info).rename(
        columns={0: "cluster_length", 1: "cluster_diameter", 2: "cluster_value"}
    )

    sum_cluster_values = sum(info[2] for info in cluster_mst_info)
    if sum_cluster_values == 0.0:
        # No valid clusters
        return 0.0, 0.0

    frac_in_clusters = df["cluster_value"].sum() / sum(all_values)

    weighted_sum = 0.0

    for mst_length, _, cluster_sum_value in cluster_mst_info:
        weighted_sum += mst_length * cluster_sum_value

    weighted_mst_cost = weighted_sum / sum_cluster_values
    return frac_in_clusters, weighted_mst_cost


def compute_centroid(points):
    """
    Given an (m x 2) array of points, return the centroid (mean x, mean y).
    """
    if len(points) == 0:
        return np.array([np.nan, np.nan])
    return np.mean(points, axis=0)


def build_adjacency_matrix(points):
    """
    Build an adjacency matrix for 'points' in 2D using pairwise distances.
    This adjacency matrix can be used in scipy.sparse.csgraph.shortest_path.

    points: (m x 2) array.
    """
    # Calculate pairwise Euclidean distances
    dist_mat = distance.cdist(points, points, metric="euclidean")
    return dist_mat  # shape (m, m)


def cluster_shortest_path_sum(points):
    """
    Compute the sum of distances along the *shortest path connecting the cluster*.
    We can interpret 'shortest path connecting the cluster' in different ways:
      - The sum of pairwise shortest paths,
      - The MST length,
      - Or other metrics.

    Here, we'll sum the shortest distances from each point to every other point
    and then possibly take a representative measure. Alternatively, you might want
    just the MST. For demonstration, we'll sum the all-pairs shortest paths.
    """
    if len(points) <= 1:
        return 0.0

    dist_mat = build_adjacency_matrix(points)
    # Compute the all-pairs shortest paths
    sp = shortest_path(dist_mat, directed=False)  # shape (m, m)
    # Sum of all shortest path distances
    # Potentially you might want a different approach (like MST or diameter).
    sum_of_shortest_paths = np.sum(sp)
    return sum_of_shortest_paths


def cluster_ratio(points, values):
    """
    Compute the ratio = (sum of values of cluster points) / (shortest_path connecting cluster).
    Returns ratio and sum_of_values.
    """
    sum_values = np.sum(values)
    sp_sum = cluster_shortest_path_sum(points)
    if sp_sum == 0:
        # If only 1 point in cluster or sp_sum is zero, define ratio carefully
        # to avoid division by zero. Let's say ratio = sum_values for a single point.
        return sum_values, sum_values, 0.0
    ratio_val = sum_values / sp_sum
    return ratio_val, sum_values, sp_sum


def reassign_points(points, values, assignments, max_iterations=10, size_threshold=10):
    """
    Repeatedly try to reassign points to clusters if beneficial.

    Parameters:
    -----------
    points : (N x 2) array of point coordinates (all points).
    values : length N array of numeric values for each point.
    assignments : length N array (cluster ID for each point).
    max_iterations : maximum times to iterate the improvement step.

    Returns:
    --------
    assignments : updated cluster assignments
    """
    # Unique cluster labels
    cluster_ids = np.unique(assignments)

    # Track changes to know if we need to keep iterating
    for iteration in range(max_iterations):
        changes_made = False

        # 1. Compute centroids for all clusters
        centroids = []
        for cid in cluster_ids:
            cluster_mask = assignments == cid
            cluster_points = points[cluster_mask]
            ctr = compute_centroid(cluster_points)
            centroids.append(ctr)
        centroids = np.array(centroids)  # shape (num_clusters, 2)

        # 2. For each cluster, find the 8 nearest other clusters
        dist_centroids = distance.cdist(
            centroids, centroids, metric="euclidean"
        )  # shape (num_clusters, num_clusters)
        # For each cluster, get the sorted neighbor cluster indices (excluding itself)
        neighbor_indices = []
        for i in range(len(cluster_ids)):
            # Sort cluster IDs by centroid distance, ignoring i itself
            sorted_neigh = np.argsort(dist_centroids[i])
            # Exclude itself (the first if it is zero distance)
            sorted_neigh = sorted_neigh[sorted_neigh != i]
            # Take the 8 nearest
            closest_8 = sorted_neigh[:8]
            neighbor_indices.append(closest_8)

        # 3. Iterate over clusters, for each point in cluster,
        #    check beneficial reassignments with each of the 8 neighbors
        for i, cid in enumerate(cluster_ids):
            cluster_mask = assignments == cid
            cluster_points = points[cluster_mask]
            cluster_values = values[cluster_mask]

            if len(cluster_points) == 0:
                continue

            # Precompute cluster ratio & sum_of_values
            orig_ratio, orig_sum_values, orig_sum_sp = cluster_ratio(
                cluster_points, cluster_values
            )

            # For each neighbor cluster
            for neighbor_cid_idx in neighbor_indices[i]:
                neighbor_cid = cluster_ids[neighbor_cid_idx]
                neighbor_mask = assignments == neighbor_cid
                neighbor_points = points[neighbor_mask]
                neighbor_values = values[neighbor_mask]

                # Precompute neighbor ratio & sum_of_values
                neighbor_orig_ratio, neighbor_orig_sum_values, neighbor_orig_sum_sp = (
                    cluster_ratio(neighbor_points, neighbor_values)
                )

                # We only need to consider reassigning points in EITHER cluster
                # to the other. Let’s focus on reassigning points from
                # the neighbor to the current cluster as well as from
                # the current cluster to the neighbor.

                # (A) Check points in neighbor cluster for reassignment to current cluster
                for idx_neighbor, pt_val_neighbor in zip(
                    np.where(neighbor_mask)[0], neighbor_values
                ):
                    point_coords = points[idx_neighbor]
                    val = values[idx_neighbor]

                    # If we remove this point from neighbor cluster
                    new_neighbor_points = np.delete(
                        neighbor_points,
                        np.where(neighbor_points == point_coords)[0],
                        axis=0,
                    )
                    new_neighbor_values = np.delete(
                        neighbor_values,
                        np.where(neighbor_points == point_coords)[0],
                        axis=0,
                    )

                    # And add it to the current cluster
                    new_cluster_points = np.vstack((cluster_points, point_coords))
                    new_cluster_values = np.concatenate((cluster_values, [val]))

                    # Recompute ratio & sums
                    new_cluster_ratio, new_cluster_sum, new_cluster_sum_sp = (
                        cluster_ratio(new_cluster_points, new_cluster_values)
                    )
                    new_neighbor_ratio, new_neighbor_sum, new_neighbor_sum_sp = (
                        cluster_ratio(new_neighbor_points, new_neighbor_values)
                    )

                    # Condition 1: ratio must INCREASE for both
                    # condition1 = (new_cluster_ratio > orig_ratio) and (new_neighbor_ratio > neighbor_orig_ratio)

                    # New Condition 1: Total sum of shortest paths must decrease
                    condition1 = (
                        new_cluster_sum_sp + new_neighbor_sum_sp
                        < orig_sum_sp + neighbor_orig_sum_sp
                    )

                    # Condition 2: both sums of values must stay > size_threshold
                    condition2 = (new_cluster_sum > size_threshold) and (
                        new_neighbor_sum > size_threshold
                    )

                    if condition1 and condition2:
                        # Reassign point
                        assignments[idx_neighbor] = cid
                        changes_made = True
                        # Update cluster_points / neighbor_points in memory
                        cluster_points = new_cluster_points
                        cluster_values = new_cluster_values
                        neighbor_points = new_neighbor_points
                        neighbor_values = new_neighbor_values
                        # Update orig ratio/sum for subsequent checks
                        orig_ratio, orig_sum_values = new_cluster_ratio, new_cluster_sum
                        neighbor_orig_ratio, neighbor_orig_sum_values = (
                            new_neighbor_ratio,
                            new_neighbor_sum,
                        )

        if not changes_made:
            # No improvements found, break early
            break

    return assignments


def find_heat_exchanger_capacity(data, distance_threshold=1.0):
    """
    Check pairwise for connections between entries where the centroid of a lower temperature band
    is within a specified distance of a higher temperature band.

    The mapping is defined as:
      - For '150-250C', the lower band is '80-150C' (for the high-temp heat exchanger)
      - For '80-150C', the lower band is '50-80C'   (for the low-temp heat exchanger)

    When the Euclidean distance between centroids (in km) is <= distance_threshold,
    the candidate lower entry's capacity (size_mw) is added to the corresponding capacity.

    Parameters:
      data (list of dict): List of entries where each dict must incl)ude:
                           'band', 'centroid' (numpy array), 'region', 'size_mw'
      distance_threshold (float): Maximum allowable distance (in km) for the connection.

    Returns:
      dict: A dictionary with the following keys:
            - 'high_temp_heat_exchanger_max_capacity': Sum of size_mw from lower entries in the
              '80-150C' band that are connected to a '150-250C' entry.
            - 'low_temp_heat_exchanger_max_capacity': Sum of size_mw from lower entries in the
              '50-80C' band that are connected to an '80-150C' entry.
    """
    # Define mapping from a given band to the lower temperature band.
    mapping = {"150-250C": "80-150C", "80-150C": "50-80C"}

    # Initialize accumulators for each heat exchanger category.
    high_temp_capacity = 0.0  # from connections: 150-250C -> 80-150C
    low_temp_capacity = 0.0  # from connections: 80-150C -> 50-80C

    for i, upper in enumerate(data):
        for j, lower in enumerate(
            data[i + 1 :], start=i + 1
        ):  # Start from i+1 to avoid duplicates
            higher_band = upper.get("band")

            if higher_band not in mapping:
                continue

            higher_centroid = upper.get("centroid")
            lower_centroid = lower.get("centroid")

            if mapping[higher_band] == lower.get("band"):

                distance = np.linalg.norm(higher_centroid - lower_centroid)
                if distance <= distance_threshold:
                    if higher_band == "150-250C":
                        high_temp_capacity += lower.get("size_mw", 0)
                    elif higher_band == "80-150C":
                        low_temp_capacity += lower.get("size_mw", 0)

            # Also check the reverse direction
            higher_band = lower.get("band")
            if higher_band in mapping and mapping[higher_band] == upper.get("band"):
                distance = np.linalg.norm(higher_centroid - lower_centroid)
                if distance <= distance_threshold:
                    if higher_band == "150-250C":
                        high_temp_capacity += upper.get("size_mw", 0)
                    elif higher_band == "80-150C":
                        low_temp_capacity += upper.get("size_mw", 0)

    return {
        "high_temp_heat_exchanger_max_capacity": high_temp_capacity,
        "low_temp_heat_exchanger_max_capacity": low_temp_capacity,
    }


def process_regional_supply_curves(
        data,
        supply_curve_step_number=3,
        demand_column="heat_demand[MW]"
        ):
    """
    Process regional geothermal data to create discretized supply curves.

    For each region and technology combination:
    - If fewer than supply_curve_step_number entries exist, just sort by capex
    - If more entries exist, discretize the supply curve by splitting into equal-sized bins
      based on cumulative demand

    Parameters:
    -----------
    data : pandas.DataFrame
        DataFrame containing regional geothermal data
    supply_curve_step_number : int, default=3
        Number of steps to discretize the supply curve into

    Returns:
    --------
    pandas.DataFrame
        Processed data with discretized supply curves
    """
    import numpy as np

    data = data.loc[data['lcoe[USD/MWh]'] < 1000.] # can occur at the edge of suitable regions and can distort the supply curve

    # Identify share columns
    share_cols = [col for col in data.columns if "share" in col]

    # Create empty list to store processed data
    processed_data = []

    # Iterate over unique regions
    for region in data.region.unique():
        region_data = data[data.region == region]

        # Iterate over unique technologies in this region
        for tech in region_data.tech.unique():
            subset = region_data[region_data.tech == tech].copy()

            # Sort by capex
            subset = subset.sort_values("capex[USD/MW]")

            # If fewer entries than supply_curve_step_number, just add numbering
            if len(subset) <= supply_curve_step_number:
                # Add numbering to tech - starting at 1 for each tech
                for i, idx in enumerate(subset.index):
                    subset.loc[idx, "tech"] = f"{tech} {i+1}"

                processed_data.append(subset)

            # If more entries, discretize the supply curve
            else:
                # Create supply curve from heat_demand and capex
                subset["cumulative_demand"] = subset[demand_column].cumsum()

                # Calculate demand-weighted average capex for original data
                original_weighted_capex = (
                    subset["capex[USD/MW]"] * subset[demand_column]
                ).sum() / subset[demand_column].sum()

                # Calculate total demand
                total_demand = subset[demand_column].sum()

                # Create equal-sized bins based on cumulative demand
                bin_size = total_demand / supply_curve_step_number
                bin_boundaries = [
                    i * bin_size for i in range(supply_curve_step_number + 1)
                ]

                # Initialize discretized data
                discretized_data = []

                # Create discretized steps based on bin boundaries
                for step in range(supply_curve_step_number):
                    step_start = bin_boundaries[step]
                    step_end = bin_boundaries[step + 1]

                    # Find rows that fall within this step
                    if step == 0:
                        # For the first step, ensure we include the first point
                        mask = subset["cumulative_demand"] <= step_end
                    else:
                        mask = (subset["cumulative_demand"] > step_start) & (
                            subset["cumulative_demand"] <= step_end
                        )

                    step_subset = subset[mask]

                    # Ensure the first bin always has data
                    if step == 0 and step_subset.empty:
                        # If first bin is empty, include at least the first point
                        step_subset = subset.iloc[[0]]

                    if not step_subset.empty:
                        # Create a new row for this step
                        new_row = {}

                        # Calculate demand-weighted average capex for this step
                        new_row["capex[USD/MW]"] = (
                            step_subset["capex[USD/MW]"]
                            * step_subset[demand_column]
                        ).sum() / step_subset[demand_column].sum()

                        # Sum heat demand for this step
                        new_row[demand_column] = step_subset[demand_column].sum()

                        # Sum total output
                        new_row["total_output[MWh]"] = step_subset[
                            "total_output[MWh]"
                        ].sum()

                        # Average for opex and share columns
                        new_row["opex[USD/MWh]"] = step_subset["opex[USD/MWh]"].mean()
                        new_row["lcoe[USD/MWh]"] = step_subset["lcoe[USD/MWh]"].mean()

                        for col in share_cols:
                            if col in step_subset.columns:
                                new_row[col] = step_subset[col].mean()

                        # Add region and tech with numbering - starting at 1 for each tech
                        new_row["region"] = region
                        new_row["tech"] = f"{tech} {step+1}"

                        discretized_data.append(new_row)

                # Convert discretized data to DataFrame
                if discretized_data:
                    discretized_df = pd.DataFrame(discretized_data)

                    # Calculate demand-weighted average capex for discretized data
                    discretized_weighted_capex = (
                        discretized_df["capex[USD/MW]"]
                        * discretized_df[demand_column]
                    ).sum() / discretized_df[demand_column].sum()

                    # Adjust capex values to ensure demand-weighted average matches original
                    if discretized_weighted_capex != 0:  # Avoid division by zero
                        adjustment_factor = (
                            original_weighted_capex / discretized_weighted_capex
                        )
                        discretized_df["capex[USD/MW]"] = (
                            discretized_df["capex[USD/MW]"] * adjustment_factor
                        )

                    processed_data.append(discretized_df)
        # break

    # Combine all processed data
    if processed_data:
        result = pd.concat(processed_data, ignore_index=True)

        # Set multi-index for region and tech
        result = result.set_index(["region", "tech"])

        return result
    else:
        return pd.DataFrame()


def process_techno_economic_data(df):
    """
    Process the techno-economic data for different technologies.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing technology parameters

    Returns:
    --------
    pandas.DataFrame
        DataFrame with processed techno-economic parameters
    """
    # Define the standardized technology names mapping
    tech_mapping = {
        "directheat_100": "directheat100degC",
        "directheat_200": "directheat200degC",
        "power_residheat_egs": "pwr_residheat80degC_egs",
        "power_residheat_hs": "pwr_residheat80degC_hs",
        "steam150_egs": "steam150degC_egs",
        "steam150_hs": "steam150degC_hs",
        "steam175_power_residheat80_egs": "steam175degC_power_residheat80degC_egs",
        "steam175_power_residheat80_hs": "steam175degC_power_residheat80degC_hs",
        "steam200_power_residheat80_egs": "steam200degC_power_residheat80degC_egs",
        "steam200_power_residheat80_hs": "steam200degC_power_residheat80degC_hs",
        "steam225_power_residheat80_egs": "steam225degC_power_residheat80degC_egs",
        "steam225_power_residheat80_hs": "steam225degC_power_residheat80degC_hs",
    }

    # Mapper for output types and their temperature bands
    tech_output_mapping = {
        "directheat100degC": {"heat": "80-150C"},
        "directheat200degC": {"heat": "150-250C"},
        "pwr_residheat80degC_egs": {"heat": "50-80C", "power": "AC"},
        "pwr_residheat80degC_hs": {"heat": "50-80C", "power": "AC"},
        "steam150degC_egs": {"steam": "80-150C", "heat": "50-80C"},
        "steam150degC_hs": {"steam": "80-150C", "heat": "50-80C"},
        "steam175degC_power_residheat80degC_egs": {
            "power": "AC",
            "heat": "50-80C",
            "steam": "80-150C",
        },
        "steam175degC_power_residheat80degC_hs": {
            "power": "AC",
            "heat": "50-80C",
            "steam": "80-150C",
        },
        "steam200degC_power_residheat80degC_egs": {
            "steam": "150-250C",
            "power": "AC",
            "heat": "50-80C",
        },  # , 'steam': '80-150C'},
        "steam200degC_power_residheat80degC_hs": {
            "steam": "150-250C",
            "power": "AC",
            "heat": "50-80C",
        },  # , 'steam': '80-150C'},
        "steam225degC_power_residheat80degC_egs": {
            "steam": "150-250C",
            "power": "AC",
            "heat": "50-80C",
        },  # , 'steam': '80-150C'},
        "steam225degC_power_residheat80degC_hs": {
            "steam": "150-250C",
            "power": "AC",
            "heat": "50-80C",
        },  # , 'steam': '80-150C'}
    }

    # Create a new DataFrame to store processed data
    result_data = {}
    result_data["geometry"] = df["geometry"]

    lifetimes = {"steam": 20, "power": 25, "directheat": 30}

    # Process each technology type
    for tech_key, std_tech_name in tech_mapping.items():
        # Extract technology type and mode
        if "_egs" in tech_key:
            mode = "egs"
        elif "_hs" in tech_key:
            mode = "hs"
        else:
            mode = None

        # Calculate original LCOE from raw data
        if tech_key.startswith("directheat"):
            temp = tech_key.split("_")[1]
            prefix = f"directheat_{temp}"
            lifetime = lifetimes["directheat"]
        elif tech_key.startswith("power_residheat"):
            prefix = f"power_residheat_{mode}"
            lifetime = lifetimes["power"]
        elif "steam" in tech_key and "power_residheat" not in tech_key:
            temp = tech_key.split("_")[0].replace("steam", "")
            prefix = f"steam{temp}_{mode}"
            lifetime = lifetimes["steam"]
        elif "steam" in tech_key and "power_residheat" in tech_key:
            parts = tech_key.split("_")
            temp = parts[0].replace("steam", "")
            prefix = f"steam{temp}_{mode}"
            lifetime = lifetimes["power"]

        # Get all capex, opex, and sales columns for this technology
        capex_cols = [
            col for col in df.columns if col.startswith(prefix) and "capex" in col
        ]
        opex_cols = [
            col for col in df.columns if col.startswith(prefix) and "opex" in col
        ]
        sales_cols = [
            col for col in df.columns if col.startswith(prefix) and "sales" in col
        ]

        if capex_cols and opex_cols and sales_cols:
            # Sum all capex, opex, and sales
            total_capex = sum(df[col] for col in capex_cols)
            total_opex = sum(df[col] for col in opex_cols)
            total_sales = sum(df[col] for col in sales_cols)

            # Calculate CRF (Capital Recovery Factor)
            discount_rate = 0.07
            crf = (
                discount_rate
                * (1 + discount_rate) ** lifetime
                / ((1 + discount_rate) ** lifetime - 1)
            )

            # Calculate original LCOE
            original_lcoe = (total_capex * crf + total_opex) / (total_sales * 8760)

            # Store the original LCOE
            # result_data[(std_tech_name, 'original_lcoe[USD/MWh]')] = original_lcoe * 1e6

        # Process based on technology type
        if tech_key.startswith("directheat"):
            # Direct heat technologies (only heat output)
            temp = tech_key.split("_")[1]

            # Column prefixes
            prefix = f"directheat_{temp}"

            # For directheat, sales, capex and opex are at the end of column names
            sales_col = f"{prefix}_sales"

            if sales_col not in df.columns:
                error_msg = (
                    f"Sales column {sales_col} not found for technology {tech_key}"
                )
                print(error_msg)
                raise ValueError(error_msg)

            # Get the appropriate lifetime for this technology
            lifetime = lifetimes["directheat"]

            # Calculate total lifetime output
            lifetime_output = df[sales_col] * 8760 * lifetime
            result_data[(std_tech_name, "total_output[MWh]")] = lifetime_output

            # Get CAPEX
            capex_col = f"{prefix}_capex"
            if capex_col not in df.columns:
                error_msg = (
                    f"CAPEX column {capex_col} not found for technology {tech_key}"
                )
                print(error_msg)
                raise ValueError(error_msg)
            result_data[(std_tech_name, "capex[USD/MW]")] = (
                df[capex_col].div(df[sales_col]).mul(1e6)
            )

            # Get OPEX
            opex_cols = [
                col for col in df.columns if col.startswith(prefix) and "opex" in col
            ]
            if not opex_cols:
                error_msg = f"No OPEX columns found for technology {tech_key} with prefix {prefix}"
                print(error_msg)
                raise ValueError(error_msg)
            # Sum all OPEX columns instead of just using the first one
            total_opex = sum(df[col] for col in opex_cols)
            result_data[(std_tech_name, "opex[USD/MWh]")] = total_opex.div(
                lifetime_output
            ).mul(1e6)

            # Only set the heat share where opex is not NaN
            heat_share_key = (
                std_tech_name,
                f"{tech_output_mapping[std_tech_name]['heat']}_share",
            )
            mask = ~total_opex.isna()
            result_data[heat_share_key] = pd.Series(np.nan, index=total_opex.index)
            result_data[heat_share_key].loc[mask] = 1.0

        elif tech_key.startswith("power_residheat"):
            # Power with residual heat (power and heat outputs)
            prefix = f"power_residheat_{mode}"

            sales_cols = [
                col
                for col in df.columns
                if col.startswith(prefix) and col.split("_")[-2] == "sales"
            ]

            output_types = [col.split("_")[-1] for col in sales_cols]

            outputs = {}
            for col, output_type in zip(sales_cols, output_types):
                if col in df.columns:
                    outputs[output_type] = df[col]

            if not outputs:
                error_msg = f"No output types found for technology {tech_key} with prefix {prefix}"
                # Show all columns that start with the tech prefix
                matching_cols = [col for col in df.columns if col.startswith(prefix)]
                print(f"Columns starting with '{prefix}':", matching_cols)
                print(error_msg)
                raise ValueError(error_msg)

            # Calculate total output across all types
            total_output = sum(outputs.values())

            # Calculate output shares relative to the total output (sum to 1)
            for output_type, value in outputs.items():
                result_data[
                    (
                        std_tech_name,
                        f"{tech_output_mapping[std_tech_name][output_type]}_share",
                    )
                ] = value.div(total_output)

            # Get the appropriate lifetime for this technology
            lifetime = lifetimes["power"]

            # Calculate total lifetime output based on the total output
            lifetime_output = total_output * 8760 * lifetime
            result_data[(std_tech_name, f"total_output[MWh]")] = lifetime_output

            # Sum CAPEX components
            capex_cols = [
                col for col in df.columns if col.startswith(prefix) and "capex" in col
            ]
            total_capex = sum(df[col] for col in capex_cols)
            result_data[(std_tech_name, f"capex[USD/MW]")] = total_capex.div(
                total_output
            ).mul(1e6)

            # Sum OPEX components
            opex_cols = [
                col for col in df.columns if col.startswith(prefix) and "opex" in col
            ]
            total_opex = sum(df[col] for col in opex_cols)
            result_data[(std_tech_name, f"opex[USD/MWh]")] = total_opex.div(
                lifetime_output
            ).mul(1e6)

        elif "steam" in tech_key and "power_residheat" not in tech_key:
            # Steam only (heat output only)
            temp = tech_key.split("_")[0].replace("steam", "")
            prefix = f"steam{temp}_{mode}"

            # Find all sales columns for this technology
            sales_cols = [
                col
                for col in df.columns
                if col.startswith(prefix) and col.split("_")[-2] == "sales"
            ]

            # Extract output types from sales column names
            output_types = [col.split("_")[-1] for col in sales_cols]

            # Create a dictionary of output types to their sales values
            outputs = {}
            for col, output_type in zip(sales_cols, output_types):
                if col in df.columns:
                    outputs[output_type] = df[col]

            if not outputs:
                error_msg = f"No output types found for technology {tech_key} with prefix {prefix}"
                # Show all columns that start with the tech prefix
                matching_cols = [col for col in df.columns if col.startswith(prefix)]
                print(f"Columns starting with '{prefix}':", matching_cols)
                print(error_msg)
                raise ValueError(error_msg)

            # Calculate total output across all types
            total_output = sum(outputs.values())

            # Calculate output shares relative to the total output (sum to 1)
            for output_type, value in outputs.items():
                # result_data[(std_tech_name, f'{output_type}_share')] = value.div(total_output)
                result_data[
                    (
                        std_tech_name,
                        f"{tech_output_mapping[std_tech_name][output_type]}_share",
                    )
                ] = value.div(total_output)

            # Get the appropriate lifetime for this technology
            lifetime = lifetimes["steam"]

            # Calculate total lifetime output based on the total output
            lifetime_output = total_output * 8760 * lifetime
            result_data[(std_tech_name, f"total_output[MWh]")] = lifetime_output

            # Sum CAPEX components
            capex_cols = [
                col for col in df.columns if col.startswith(prefix) and "capex" in col
            ]
            total_capex = sum(df[col] for col in capex_cols)
            result_data[(std_tech_name, f"capex[USD/MW]")] = total_capex.div(
                total_output
            ).mul(1e6)

            # Sum OPEX components
            opex_cols = [
                col for col in df.columns if col.startswith(prefix) and "opex" in col
            ]
            total_opex = sum(df[col] for col in opex_cols)
            result_data[(std_tech_name, f"opex[USD/MWh]")] = total_opex.div(
                lifetime_output
            ).mul(1e6)

        elif "steam" in tech_key and "power_residheat" in tech_key:
            # Combined steam, power, and residual heat
            parts = tech_key.split("_")
            temp = parts[0].replace("steam", "")
            prefix = f"steam{temp}_{mode}"

            # Find all sales columns for this technology
            sales_cols = [
                col
                for col in df.columns
                if col.startswith(prefix) and col.split("_")[-2] == "sales"
            ]

            # Extract output types from sales column names
            output_types = [col.split("_")[-1] for col in sales_cols]

            # Create a dictionary of output types to their sales values
            outputs = {}
            for col, output_type in zip(sales_cols, output_types):
                if col in df.columns:
                    outputs[output_type] = df[col]

            if not outputs:
                error_msg = f"No output types found for technology {tech_key} with prefix {prefix}"
                # Show all columns that start with the tech prefix
                matching_cols = [col for col in df.columns if col.startswith(prefix)]
                print(f"Columns starting with '{prefix}':", matching_cols)
                print(error_msg)
                raise ValueError(error_msg)

            # Calculate total output across all types
            total_output = sum(outputs.values())

            # Calculate output shares relative to the total output (sum to 1)
            for output_type, value in outputs.items():
                result_data[
                    (
                        std_tech_name,
                        f"{tech_output_mapping[std_tech_name][output_type]}_share",
                    )
                ] = value.div(total_output)

            # Get the appropriate lifetime for this technology (using power as primary output)
            lifetime = lifetimes["power"]

            # Calculate total lifetime output based on the total output
            lifetime_output = total_output * 8760 * lifetime
            result_data[(std_tech_name, f"total_output[MWh]")] = lifetime_output

            # Sum CAPEX components
            capex_cols = [
                col for col in df.columns if col.startswith(prefix) and "capex" in col
            ]
            total_capex = sum(df[col] for col in capex_cols)
            result_data[(std_tech_name, f"capex[USD/MW]")] = total_capex.div(
                total_output
            ).mul(1e6)

            # Sum OPEX components
            opex_cols = [
                col for col in df.columns if col.startswith(prefix) and "opex" in col
            ]
            total_opex = sum(df[col] for col in opex_cols)
            result_data[(std_tech_name, f"opex[USD/MWh]")] = total_opex.div(
                lifetime_output
            ).mul(1e6)

    # Create DataFrame with MultiIndex columns
    result_df = pd.DataFrame(result_data)

    return result_df


if __name__ == "__main__":

    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_industrial_heating_demands",
            simpl="",
            clusters=10,
        )

    #############     OVERLAYS DIFFERENT GEOTHERMAL TECHS TO HAVE FOR EACH REGION THE
    #############     COST-OPTIMAL TECH FOR EACH TEMPERATURE BAND OF DEMAND

    # these are the geothermal input files, as passed in the data transfer pypsa_inputs_draft_20250403
    tif_files = {
        name: fn for name, fn in snakemake.input.items() if fn.endswith(".tif")
    }
    file_name_transformer = lambda x: "-".join(str(x).split("/")[-2:]).replace(
        ".tif", ""
    )

    gdf = tif_to_gdf(tif_files.values(), name_transformer=file_name_transformer)
    gdf = gdf.rename(
        columns={file_name_transformer(item): key for key, item in tif_files.items()}
    )

    gdf = process_techno_economic_data(gdf)

    gdf.rename(columns={"geometry": ("geometry", "")}, inplace=True)
    gdf.columns = pd.MultiIndex.from_tuples(gdf.columns)

    tech_levels = gdf.columns.get_level_values(0).unique()

    new_parts = [gdf[["geometry"]]]

    idx = pd.IndexSlice

    for tech in tech_levels:
        if tech == "geometry":
            continue

        new_parts.append(gdf.loc[:, idx[tech, :]].dropna())

    gdf = pd.concat(new_parts, axis=1)

    for tech in gdf.columns.get_level_values(0).unique()[::-1]:
        # Skip geometry column
        if tech == "geometry":
            continue

        # Get all columns for this technology
        idx = pd.IndexSlice
        tech_data = gdf.loc[:, idx[tech, :, :]].dropna()

        if tech_data.empty:
            continue

        # Calculate statistics, ignoring NaN values
        tech_averages = tech_data.mean(skipna=True)
        tech_medians = tech_data.median(skipna=True)
        tech_mins = tech_data.min(skipna=True)
        tech_maxs = tech_data.max(skipna=True)

        # Define discount rate and lifetimes for LCOE calculation
        discount_rate = 0.07

        # Calculate LCOE for each technology type
        lcoe_values = {}

        # Determine which type this technology is
        tech_type = None
        if "steam" in tech.lower():
            tech_type = "steam"
        elif "power" in tech.lower():
            tech_type = "power"
        elif "heat" in tech.lower():
            tech_type = "directheat"

        lifetimes = {"steam": 20, "power": 25, "directheat": 30}

        # Find capex and opex columns for this technology
        capex_cols = [
            col
            for col in tech_data.columns.get_level_values(1)
            if "capex[USD/MW" in col
        ]
        opex_cols = [
            col
            for col in tech_data.columns.get_level_values(1)
            if "opex[USD/MWh" in col
        ]
        shares_cols = [
            col for col in tech_data.columns.get_level_values(1) if "share" in col
        ]

        # Get the first capex and opex columns (they should be for the same output type)
        assert len(capex_cols) == 1
        assert len(opex_cols) == 1

        capex_col = capex_cols[0]
        opex_col = opex_cols[0]

        # Extract the output type from the column name
        if not shares_cols:
            output_type = [tech]
        else:
            output_type = [col.split("_")[0] for col in shares_cols]

        lifetime = lifetimes.get(tech_type, 25)  # Default to 25 if type not found

        # Calculate CRF (Capital Recovery Factor)
        crf = (
            discount_rate
            * (1 + discount_rate) ** lifetime
            / ((1 + discount_rate) ** lifetime - 1)
        )

        # Calculate LCOE for each row: LCOE = (CAPEX * CRF + OPEX)
        if shares_cols:
            shares = tech_data.loc[:, idx[:, shares_cols]].sum(axis=1)
            assert np.allclose(shares, 1)
        else:
            shares = 1
        lcoe = (
            tech_data.loc[:, idx[:, capex_col]].values.flatten() * crf / 8760
            + tech_data.loc[:, idx[:, opex_col]].values.flatten()
        )

        lcoe = pd.Series(lcoe, index=tech_data.index)

        # Add LCOE values to the gdf dataframe
        gdf.loc[lcoe.index, idx[tech, f"lcoe[USD/MWh]"]] = lcoe

        # Calculate statistics for LCOE
        lcoe_values = {
            "mean": lcoe.mean(skipna=True),
            "median": lcoe.median(skipna=True),
            "min": lcoe.min(skipna=True),
            "max": lcoe.max(skipna=True),
        }

        print(f"\n{tech}:")
        for metric in tech_averages.index.get_level_values(1):

            # Format the output based on the metric type
            if "capex" in metric:
                print(
                    "  {}: Mean=${:.2f}, Median=${:.2f}, Min=${:.2f}, Max=${:.2f}".format(
                        metric,
                        float(tech_averages.loc[idx[:, metric]].iloc[0]),
                        float(tech_medians.loc[idx[:, metric]].iloc[0]),
                        float(tech_mins.loc[idx[:, metric]].iloc[0]),
                        float(tech_maxs.loc[idx[:, metric]].iloc[0]),
                    )
                )
            elif "opex" in metric:
                print(
                    "  {}: Mean=${:.4f}, Median=${:.4f}, Min=${:.4f}, Max=${:.4f}".format(
                        metric,
                        float(tech_averages.loc[idx[:, metric]].iloc[0]),
                        float(tech_medians.loc[idx[:, metric]].iloc[0]),
                        float(tech_mins.loc[idx[:, metric]].iloc[0]),
                        float(tech_maxs.loc[idx[:, metric]].iloc[0]),
                    )
                )
            elif "share" in metric:
                print(
                    "  {}: Mean={:.4f}, Median={:.4f}, Min={:.4f}, Max={:.4f}".format(
                        metric,
                        float(tech_averages.loc[idx[:, metric]].iloc[0]),
                        float(tech_medians.loc[idx[:, metric]].iloc[0]),
                        float(tech_mins.loc[idx[:, metric]].iloc[0]),
                        float(tech_maxs.loc[idx[:, metric]].iloc[0]),
                    )
                )
            elif "output" in metric or "max" in metric:
                continue
            else:
                continue

        # Print LCOE if calculated
        if lcoe_values:
            print(
                "  LCOE[USD/MWh] ({}, {} years): Mean=${:.2f}, Median=${:.2f}, Min=${:.2f}, Max=${:.2f}".format(
                    tech_type,
                    lifetimes[tech_type],
                    float(lcoe_values["mean"]),
                    float(lcoe_values["median"]),
                    float(lcoe_values["min"]),
                    float(lcoe_values["max"]),
                )
            )

            # Print original LCOE if available
            if "original_lcoe[USD/MWh]" in tech_data.columns.get_level_values(1):
                original_lcoe = tech_data.loc[:, idx[:, "original_lcoe[USD/MWh]"]]
                original_lcoe_values = {
                    "mean": original_lcoe.mean(skipna=True).iloc[0],
                    "median": original_lcoe.median(skipna=True).iloc[0],
                    "min": original_lcoe.min(skipna=True).iloc[0],
                    "max": original_lcoe.max(skipna=True).iloc[0],
                }
                print(
                    "  Original LCOE[USD/MWh]: Mean=${:.2f}, Median=${:.2f}, Min=${:.2f}, Max=${:.2f}".format(
                        float(original_lcoe_values["mean"]),
                        float(original_lcoe_values["median"]),
                        float(original_lcoe_values["min"]),
                        float(original_lcoe_values["max"]),
                    )
                )

    print("-" * 80)

    print(gdf.columns.get_level_values(0).unique())

    demand_gdf = prepare_demand_data(snakemake.input["demand_data"])
    logger.warning(
        "Inconsistency between temperature ranges in the data and this scripts."
        "This script has 50-80, 80-150, 150-250."
        "The data has 0-49, 50-99, 100-149, 150-199, 200-249, 250-299, 300-349, 350-399, 400-449, >450."
    )

    regions = gpd.read_file(snakemake.input["regions"]).set_index("name")

    egs_params = snakemake.params["enhanced_geothermal"]
    max_network_diameter = egs_params["max_network_diameter"]  # km, default 20.
    min_network_average_capacity = egs_params[
        "min_network_average_capacity"
    ]  # MWh, default 10.
    max_network_average_capacity = egs_params[
        "max_network_average_capacity"
    ]  # MWh, default 30.

    n_cost_steps = egs_params["industrial_heating_n_cost_steps"]  # 2

    supply_curve_columns = [
        # "capex[$/MW]",
        # "avail_capacity[MW]",
        # "opex[$/MWh]",
        # "temperature",
    ]

    final_demands = pd.DataFrame(
        index=regions.index,
        columns=["demand(50-80C)[MW]", "demand(80-150C)[MW]", "demand(150-250C)[MW]"],
    )

    # Create a mapping from technology to applicable demand types
    # Create a mapping from temperature bands to applicable technologies
    demand_to_tech_mapping = {
        "50-80C": [
            "directheat100degC",
            "pwr_residheat80degC_egs",
            "pwr_residheat80degC_hs",
        ],
        "80-150C": [
            "directheat100degC",
            "steam150degC_egs",
            "steam150degC_hs",
        ],
        "150-250C": [
            "directheat200degC",
            "steam175degC_power_residheat80degC_egs",
            "steam175degC_power_residheat80degC_hs",
            "steam200degC_power_residheat80degC_egs",
            "steam200degC_power_residheat80degC_hs",
            "steam225degC_power_residheat80degC_egs",
            "steam225degC_power_residheat80degC_hs",
        ],
    }

    logger.info("preprocessed EGS data!")
    #############     CLUSTERING OF DEMANDS

    regional_supplies = list()
    heat_exchanger_max_capacities = list()

    # Convert string geometries to shapely objects
    if isinstance(gdf["geometry"].iloc[0], str):
        from shapely import wkt

        gdf["geometry"] = gdf["geometry"].apply(lambda x: wkt.loads(x))
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")

    regional_supply_shapes = pd.Series(index=regions.index)

    total_results = []
    total_clusters = []

    for region, geometry in tqdm(regions["geometry"].items()):

        # regional_supply = pd.DataFrame(columns=supply_curve_columns)
        regional_supply = []

        ss = demand_gdf.loc[demand_gdf["geometry"].within(geometry)]

        temp_bands = [
            ("50-80C", lambda x: (x >= 50) & (x <= 80)),
            ("80-150C", lambda x: (x > 80) & (x <= 150)),
            ("150-250C", lambda x: (x > 150) & (x <= 250)),
        ]

        band_centroids = []  # track centroids for heat exchanger capacity

        for band, filter in temp_bands:

            band_data = ss.loc[filter(ss["temperature"])]
            final_demands.loc[region, f"demand({band})[MW]"] = band_data[
                "avg_demand"
            ].sum()

            if band_data.empty:
                logger.warning(f"No data for {region}")
                continue

            lonlat_coords = band_data[["y", "x"]].values
            utm_coords = coords_to_relative_utm(lonlat_coords)
            avg_demand = band_data["avg_demand"].values

            clusters, cluster_info = greedy_clustering(
                utm_coords,
                avg_demand,
                max_network_diameter,
                min_sum_value=min_network_average_capacity,
                max_sum_value=max_network_average_capacity,
            )

            frac_in_clusters, w_mst_cost = evaluate_solution(cluster_info, avg_demand)

            logger.info(
                f"Region: {region}, band: {band}, Fraction of values in clusters: {frac_in_clusters:.2f}, Weighted MST cost: {w_mst_cost:.2f}"
            )
            logger.info(
                f"Total EGS applicable demand: {sum(band_data['total_demand'] * frac_in_clusters) / 1e3:.2f} GWh"
            )

            assignments = np.ones(len(utm_coords), dtype=int) * -1

            for i, cluster in enumerate(clusters):
                for idx in cluster:
                    assignments[idx] = i

            if len(utm_coords[assignments != -1]):

                new_clustering = reassign_points(
                    utm_coords[assignments != -1],
                    utm_coords[assignments != -1],
                    assignments[assignments != -1],
                    max_iterations=10,
                    size_threshold=egs_params[
                        "min_network_average_capacity"
                    ],  # MWh, default 10.
                )

                new_clusters = []

                for label in sorted(pd.Series(new_clustering).value_counts().index):
                    new_clusters.append(np.where(new_clustering == label)[0].tolist())

                # map indices back to original
                old_indexes = pd.Series(np.where(assignments != -1)[0])
                fixed_clustering = []

                for c in new_clusters:

                    fixed_c = []
                    for entry in c:
                        fixed_c.append(old_indexes.loc[entry])

                    fixed_clustering.append(fixed_c)

                new_clusters = fixed_clustering

            else:
                new_clusters = clusters

            piping_cost = egs_params["piping_cost"]  # $/km, default 1_300_000$/km

            total_cluster_size = 0
            for cluster in new_clusters:

                # Create MultiPoint from cluster coordinates
                from shapely.geometry import (
                    LineString,
                    MultiLineString,
                    MultiPoint,
                    Point,
                )

                cluster_points = MultiPoint(
                    [Point(coord) for coord in lonlat_coords[cluster]]
                )

                # Create minimum spanning tree as MultiLineString

                # Calculate pairwise distances between all points in the cluster
                dist_matrix = distance_matrix(utm_coords[cluster], utm_coords[cluster])

                # Get MST as a CSR matrix
                mst_csr = minimum_spanning_tree(dist_matrix)

                # Convert to array and find edges
                mst_array = mst_csr.toarray()
                mst_edges = []
                rows, cols = np.where(mst_array > 0)
                for i, j in zip(rows, cols):
                    if i < j:  # Only add each edge once (avoid duplicates)
                        mst_edges.append((i, j))

                # Create lines for the MST
                mst_lines = [
                    LineString([lonlat_coords[cluster][i], lonlat_coords[cluster][j]])
                    for i, j in mst_edges
                ]

                # If no edges (single point), create empty MultiLineString
                if len(mst_lines) == 0:
                    cluster_mst = MultiLineString([])
                else:
                    cluster_mst = MultiLineString(mst_lines)

                # Append both to total_clusters with labels
                total_clusters.append(
                    {
                        "points": cluster_points,
                        "mst": cluster_mst,
                        "region": region,
                        "band": band,
                    }
                )

                cluster_supply = pd.Series()
                cluster_sp, _ = compute_mst_cost_and_diameter(utm_coords[cluster])
                cluster_size = sum(avg_demand[cluster])
                total_cluster_size += cluster_size

                query_point = Point(lonlat_coords[cluster].mean(axis=0)[::-1])

                # Create a buffer of 0.5 degrees around the query point
                buffer_distance = 0.5  # in degrees
                buffered_point = query_point.buffer(buffer_distance)

                # Select entries where the geometry intersects with the buffered point
                geothermal_subset = gdf.loc[gdf.geometry.intersects(buffered_point)]

                ### This is where the most cost-optimal geothermal supply is selected
                if len(geothermal_subset) == 0:
                    logger.warning(f"No geothermal data for {region}, band: {band}")
                    continue
                else:
                    geothermal_subset = geothermal_subset.iloc[0]

                idx = pd.IndexSlice

                cluster_supply = (
                    geothermal_subset.loc[idx[demand_to_tech_mapping[band], :]]
                    .dropna()
                    .unstack()
                    .replace(np.nan, 0)
                )

                if cluster_supply.empty:
                    logger.info(f'Temperature band {band} cannot be met in {region}')
                    continue
            
                # The current method only clusters heat demands of the same temperature band,
                # and therefore it can not be assumed that heat outputs on different bands can be used at the same site.
                # The following code therefore either assumes that heat output of higher temperature could also be used for a lower temperature, or removes that share of output and adjusts the LCOE accordingly.
                for other in ['80-150C', '150-250C']:
                    if band == '50-80C' and f'{other}_share' in cluster_supply.columns:
                        for idx in cluster_supply.index:
                            if '50-80C_share' in cluster_supply.columns:
                                cluster_supply.at[idx, '50-80C_share'] += cluster_supply.at[idx, f'{other}_share']
                                cluster_supply.at[idx, f'{other}_share'] = 0.
                            else:
                                cluster_supply.at[idx, '50-80C_share'] = cluster_supply.at[idx, f'{other}_share']
                                cluster_supply.at[idx, f'{other}_share'] = 0.

                if (
                    cluster_supply.empty
                    or f"{band}_share" not in cluster_supply.columns
                    ):
                    logger.info(f'Temperature band {band} cannot be met in {region}')
                    continue

                heat_bands = ['50-80C', '80-150C', '150-250C']
                others = [b for b in heat_bands if b != band]

                for rowname in cluster_supply.index:
                    removed_share = 0
                    for o in others:
                        if f'{o}_share' in cluster_supply.columns:
                            removed_share += cluster_supply.loc[rowname, f'{o}_share']
                            cluster_supply.loc[rowname, f'{o}_share'] = 0.
                    
                    if removed_share > 0:
                        cluster_supply.loc[rowname, 'lcoe[USD/MWh]'] *= 1 / (1 - removed_share)
                
                logger.info(cluster_supply)

                cluster_supply = cluster_supply.sort_values(
                    by=["lcoe[USD/MWh]"]
                ).iloc[[0]]

                cluster_supply.loc[cluster_supply.index[0], ["heat_demand[MW]"]] = (
                    cluster_size
                )
                cluster_supply.loc[cluster_supply.index[0], "capex[USD/MW]"] = (
                    cluster_supply.loc[cluster_supply.index[0], "capex[USD/MW]"]
                    * cluster_size
                    + piping_cost * cluster_sp
                ) / cluster_size

                regional_supply.append(cluster_supply)

                cluster_points = utm_coords[cluster]
                centroid = compute_centroid(cluster_points)

                appendage = {
                    "region": region,
                    "centroid": centroid,
                    "size_mw": cluster_size,
                    "band": band,
                }
                band_centroids.append(appendage)

        if not regional_supply:
            continue

        regional_supply = pd.concat(regional_supply).replace(np.nan, 0)

        regional_supply_shapes.loc[region] = len(regional_supply)

        os.makedirs("hold", exist_ok=True)

        regional_supply.to_csv(f"hold/regional_supply_{region}.csv")

        regional_supply.loc[:, "region"] = region
        total_results.append(regional_supply)

        # heat_exchanger_capacity = pd.Series(
        #     find_heat_exchanger_capacity(band_centroids, distance_threshold=2),
        #     name=region,
        # )

    total_results = (
        pd.concat(total_results).reset_index().rename(columns={"index": "tech"})
    )

    ### discretizes supply curves
    total_results = process_regional_supply_curves(total_results).replace(np.nan, 0)

    # df = pd.DataFrame((total_clusters))
    # df.to_csv("clusters.csv")

    total_results.to_csv(snakemake.output["industrial_heating_egs_supply_curves"])
    final_demands.to_csv(snakemake.output["industrial_heating_demands"])

    # pd.concat(heat_exchanger_max_capacities, axis=1).to_csv(
    #     snakemake.output["heat_exchanger_capacity"]
    # )
