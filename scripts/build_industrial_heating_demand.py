
import logging

logger = logging.getLogger(__name__)

import random
import numpy as np
import pandas as pd
from tqdm import tqdm
import geopandas as gpd
import matplotlib.pyplot as plt

from copy import deepcopy
from scipy.spatial import distance
from scipy.spatial import distance_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.csgraph import minimum_spanning_tree

from _industry_heat_helpers import coords_to_relative_utm


def prepare_demand_data(fn):

    df = gpd.read_file(fn)

    id_cols = [
        # 'names_primary', 'geometry',
        # 'industrial_type', 'industrial_sub_type'
        'FACILITY_ID', 'FACILITY_NAME', 'geometry',
    ]

    heating_cols = [
        col for col in df.columns
        if 'Heating (BBtu)' in col and 'total' not in col and '<250' not in col
    ]

    # Melt the dataframe to unpivot heating columns
    df_melted = df.melt(
        id_vars=id_cols,
        value_vars=heating_cols,
        var_name='Temperature Range',
        value_name='Heating Demand (BBtu)'
    )

    df_melted = df_melted[df_melted['Heating Demand (BBtu)'] != 0]

    df_melted['Temperature Range'] = df_melted['Temperature Range'].str.replace(
        'Heating (BBtu), ', '', regex=False
    )

    # df_melted['Heating Demand (MWh)'] = df_melted['Heating Demand (BBtu)'] * 293.07107
    df_melted['Heating Demand (MWh)'] = df_melted['Heating Demand (BBtu)'] * 2.9307e-7

    df_melted = df_melted.drop(columns=['Heating Demand (BBtu)'])

    df_melted = df_melted[
        id_cols + ['Temperature Range', 'Heating Demand (MWh)']
    ]

    cutoff = 5000 # MWh; cutoff of 5GWh removes around 1% of total demand, and 50% of datapoints

    df_melted = df_melted.loc[df_melted['Heating Demand (MWh)'] > cutoff]

    df_melted = df_melted.reset_index(drop=True)

    temp_mapper = {
        '0 to 49°C': 50,
        '50 to 99°C': 100,
        '100 to 149°C': 150,
        '150 to 199°C': 200,
        '200 to 249°C': 250,
        '250 to 299°C': 300,
        '300 to 349°C': 350,
        '350 to 399°C': 400,
        '400 to 449°C': 450,
        '>450°C': 500,
    }

    for old, new in temp_mapper.items():
        df_melted.replace(old, new, inplace=True)

    df_melted.rename(columns={
        'Temperature Range': 'temperature',
        'Heating Demand (MWh)': 'total_demand',
        }, inplace=True)
    df_melted['avg_demand'] = df_melted['total_demand'] / 8760

    df_melted = df_melted.loc[df_melted['temperature'] <= 250]

    df_melted['x'] = df_melted['geometry'].x
    df_melted['y'] = df_melted['geometry'].y

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


def greedy_clustering(
        points,
        values,
        T,
        min_sum_value=10.0,
        max_sum_value=np.inf
        ):
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
    assigned = [False]*n_points
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

    df = (
        pd.DataFrame(cluster_mst_info)
        .rename(columns={0: 'cluster_length', 1: 'cluster_diameter', 2: 'cluster_value'})
    )

    sum_cluster_values = sum(info[2] for info in cluster_mst_info)
    if sum_cluster_values == 0.0:
        # No valid clusters
        return 0.0, 0.0

    frac_in_clusters = df['cluster_value'].sum() / sum(all_values)

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
    dist_mat = distance.cdist(points, points, metric='euclidean')
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
        return sum_values, sum_values, 0.
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
            cluster_mask = (assignments == cid)
            cluster_points = points[cluster_mask]
            ctr = compute_centroid(cluster_points)
            centroids.append(ctr)
        centroids = np.array(centroids)  # shape (num_clusters, 2)
        
        # 2. For each cluster, find the 8 nearest other clusters
        dist_centroids = distance.cdist(centroids, centroids, metric='euclidean')  # shape (num_clusters, num_clusters)
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
            cluster_mask = (assignments == cid)
            cluster_points = points[cluster_mask]
            cluster_values = values[cluster_mask]
            
            if len(cluster_points) == 0:
                continue

            # Precompute cluster ratio & sum_of_values
            orig_ratio, orig_sum_values, orig_sum_sp = cluster_ratio(cluster_points, cluster_values)
            
            # For each neighbor cluster
            for neighbor_cid_idx in neighbor_indices[i]:
                neighbor_cid = cluster_ids[neighbor_cid_idx]
                neighbor_mask = (assignments == neighbor_cid)
                neighbor_points = points[neighbor_mask]
                neighbor_values = values[neighbor_mask]
                
                # Precompute neighbor ratio & sum_of_values
                neighbor_orig_ratio, neighbor_orig_sum_values, neighbor_orig_sum_sp = cluster_ratio(neighbor_points, neighbor_values)
                
                # We only need to consider reassigning points in EITHER cluster
                # to the other. Let’s focus on reassigning points from 
                # the neighbor to the current cluster as well as from 
                # the current cluster to the neighbor.

                # (A) Check points in neighbor cluster for reassignment to current cluster
                for idx_neighbor, pt_val_neighbor in zip(np.where(neighbor_mask)[0], neighbor_values):
                    point_coords = points[idx_neighbor]
                    val = values[idx_neighbor]
                    
                    # If we remove this point from neighbor cluster
                    new_neighbor_points = np.delete(neighbor_points, 
                                                    np.where(neighbor_points == point_coords)[0], axis=0)
                    new_neighbor_values = np.delete(neighbor_values, 
                                                    np.where(neighbor_points == point_coords)[0], axis=0)
                    
                    # And add it to the current cluster
                    new_cluster_points = np.vstack((cluster_points, point_coords))
                    new_cluster_values = np.concatenate((cluster_values, [val]))
                    
                    # Recompute ratio & sums
                    new_cluster_ratio, new_cluster_sum, new_cluster_sum_sp = cluster_ratio(new_cluster_points, new_cluster_values)
                    new_neighbor_ratio, new_neighbor_sum, new_neighbor_sum_sp = cluster_ratio(new_neighbor_points, new_neighbor_values)
                    
                    # Condition 1: ratio must INCREASE for both
                    # condition1 = (new_cluster_ratio > orig_ratio) and (new_neighbor_ratio > neighbor_orig_ratio)

                    # New Condition 1: Total sum of shortest paths must decrease
                    condition1 = new_cluster_sum_sp + new_neighbor_sum_sp < orig_sum_sp + neighbor_orig_sum_sp

                    # Condition 2: both sums of values must stay > size_threshold
                    condition2 = (new_cluster_sum > size_threshold) and (new_neighbor_sum > size_threshold)
                    
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
                        neighbor_orig_ratio, neighbor_orig_sum_values = new_neighbor_ratio, new_neighbor_sum


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
      data (list of dict): List of entries where each dict must include:
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
    mapping = {
        '150-250C': '80-150C',
        '80-150C': '50-80C'
    }
    
    # Initialize accumulators for each heat exchanger category.
    high_temp_capacity = 0.0  # from connections: 150-250C -> 80-150C
    low_temp_capacity = 0.0   # from connections: 80-150C -> 50-80C
    
    for i, upper in enumerate(data):
        for j, lower in enumerate(data[i+1:], start=i+1):  # Start from i+1 to avoid duplicates
            higher_band = upper.get('band')
            
            if higher_band not in mapping:
                continue

            higher_centroid = upper.get('centroid')
            lower_centroid = lower.get('centroid')

            if mapping[higher_band] == lower.get('band'):

                distance = np.linalg.norm(higher_centroid - lower_centroid)
                if distance <= distance_threshold:
                    if higher_band == '150-250C':
                        high_temp_capacity += lower.get('size_mw', 0)
                    elif higher_band == '80-150C':
                        low_temp_capacity += lower.get('size_mw', 0)

            # Also check the reverse direction
            higher_band = lower.get('band')
            if higher_band in mapping and mapping[higher_band] == upper.get('band'):
                distance = np.linalg.norm(higher_centroid - lower_centroid)
                if distance <= distance_threshold:
                    if higher_band == '150-250C':
                        high_temp_capacity += upper.get('size_mw', 0)
                    elif higher_band == '80-150C':
                        low_temp_capacity += upper.get('size_mw', 0)

    return {
        'high_temp_heat_exchanger_max_capacity': high_temp_capacity,
        'low_temp_heat_exchanger_max_capacity': low_temp_capacity
    }



if __name__ == '__main__':

    gdf = prepare_demand_data(snakemake.input['demand_data'])
    logger.warning(
        'Inconcistency between temperature ranges in the data and this scripts.'
        'This script has 50-80, 80-150, 150-250.'
        'The data has 0-49, 50-99, 100-149, 150-199, 200-249, 250-299, 300-349, 350-399, 400-449, >450.'
    )

    regions = gpd.read_file(snakemake.input['regions']).set_index('name')

    egs_params = snakemake.params['enhanced_geothermal']
    max_network_diameter = egs_params['max_network_diameter']
    min_network_average_capacity = egs_params['min_network_average_capacity']
    max_network_average_capacity = egs_params['max_network_average_capacity']

    n_cost_steps = egs_params['industrial_heating_n_cost_steps']

    supply_curve_columns = ['capex[$/MW]', 'avail_capacity[MW]', 'opex[$/MWh]', 'temperature']
    
    final_demands = pd.DataFrame(
        index=regions.index,
        columns=['demand(50-80C)[MW]', 'demand(80-150C)[MW]', 'demand(150-250C)[MW]']
        )

    regional_supplies = list()
    heat_exchanger_max_capacities = list()

    for region, geometry in tqdm(regions['geometry'].items()):

        regional_supply = pd.DataFrame(columns=supply_curve_columns)
        ss = gdf.loc[gdf['geometry'].within(geometry)]

        temp_bands = [
            ('50-80C', lambda x: (x >= 50) & (x <= 80)),
            ('80-150C', lambda x: (x > 80) & (x <= 150)),
            ('150-250C', lambda x: (x > 150) & (x <= 250))
        ]

        band_centroids = [] # track centroids for heat exchanger capacity

        for band, filter in temp_bands:

            band_data = ss.loc[filter(ss['temperature'])]
            final_demands.loc[region, f'demand({band})[MW]'] = band_data['avg_demand'].sum()

            if band_data.empty:
                logger.warning(f"No data for {region}")
                continue

            utm_coords = coords_to_relative_utm(band_data[['y', 'x']].values)
            avg_demand = band_data['avg_demand'].values

            clusters, cluster_info = greedy_clustering(
                utm_coords,
                avg_demand,
                max_network_diameter,
                min_sum_value=min_network_average_capacity,
                max_sum_value=max_network_average_capacity
                )

            frac_in_clusters, w_mst_cost = evaluate_solution(cluster_info, avg_demand)

            logger.info(f"Region: {region}, band: {band}, Fraction of values in clusters: {frac_in_clusters:.2f}, Weighted MST cost: {w_mst_cost:.2f}")
            logger.info(f"Total EGS applicable demand: {sum(band_data['total_demand'] * frac_in_clusters) / 1e3:.2f} GWh")

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
                    size_threshold=egs_params['min_network_average_capacity']
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

            piping_cost = egs_params['piping_cost']

            def dummy_egs_query(cluster):
                return 0 # $/MWth

            total_cluster_size = 0
            for cluster in new_clusters:

                cluster_supply = pd.Series()
                cluster_sp, _ = compute_mst_cost_and_diameter(utm_coords[cluster])
                cluster_size = sum(avg_demand[cluster]) 
                total_cluster_size += cluster_size

                cluster_supply.loc['capex[$/MW]'] = (
                    dummy_egs_query(cluster) * cluster_size +
                    piping_cost * cluster_sp) / cluster_size

                cluster_supply.loc['avail_capacity[MW]'] = cluster_size            
                cluster_supply.loc['opex[$/MWh]'] = 0.0
                cluster_supply.loc['temperature'] = band

                if cluster_supply.empty:
                    continue

                cluster_points = utm_coords[cluster]
                centroid = compute_centroid(cluster_points)
                
                appendage = {
                    'region': region,
                    'centroid': centroid,
                    'size_mw': cluster_size,
                    'band': band
                    }
                band_centroids.append(appendage)
            
                regional_supply = pd.concat(
                    (regional_supply, cluster_supply.to_frame().T),
                    ignore_index=True
                    )

        heat_exchanger_capacity = pd.Series(
            find_heat_exchanger_capacity(band_centroids, distance_threshold=2), name=region
        )

        regional_supply.loc[:,'capex[$/MW]'] = regional_supply.loc[:,'capex[$/MW]'].round(3)

        # Process each temperature band separately
        temperature_bands = regional_supply['temperature'].unique()
        processed_supplies = []

        for temp_band in temperature_bands:

            temp_supply = regional_supply[regional_supply['temperature'] == temp_band].copy()
            
            temp_supply = temp_supply.dropna()

            temp_supply = temp_supply.groupby('capex[$/MW]').agg(
                {'avail_capacity[MW]': 'sum', 'opex[$/MWh]': 'mean', 'temperature': 'first'}
                ).reset_index()

            if temp_supply.empty:
                continue

            if len(temp_supply) <= n_cost_steps: # nothing needs to be done
                temp_supply['level'] = temp_supply['capex[$/MW]'].values
                pass

            else:
                logger.warning('Check this method once real data is used!!')
                bins = pd.Series(
                        np.linspace(
                            temp_supply['capex[$/MW]'].min(),
                            temp_supply['capex[$/MW]'].max(),
                            n_cost_steps+1
                            )
                        )

                labels = (
                    bins
                    .rolling(2)
                    .mean()
                    .dropna()
                    .tolist()
                )

                temp_supply['level'] = pd.cut(
                    temp_supply['capex[$/MW]'],
                    bins=bins,
                    labels=labels,
                    duplicates='drop'
                    )

            temp_supply = temp_supply.dropna()

            temp_supply = (
                temp_supply
                .groupby('level', observed=False)
                [['avail_capacity[MW]', 'opex[$/MWh]', 'temperature']]
                .agg({
                    'avail_capacity[MW]': 'sum',
                    'opex[$/MWh]': 'mean',
                    'temperature': 'first'
                    })
                .reset_index()
                .rename(columns={'level': 'capex[$/MW]'})
            )

            temp_supply = temp_supply.dropna()
            temp_supply.index = pd.MultiIndex.from_arrays(
                [
                    [region] * len(temp_supply),
                    temp_supply['temperature'].values,
                    range(len(temp_supply))
                ],
                names=['region', 'temperature', 'cost_step']
            )

            temp_supply.drop(columns=['temperature'], inplace=True)
            processed_supplies.append(temp_supply)
        
        if processed_supplies:
            regional_supply = pd.concat(processed_supplies)
            regional_supplies.append(regional_supply)

        heat_exchanger_max_capacities.append(
            pd.Series(  
                heat_exchanger_capacity, name=region
                )
        )

    logger.warning('Assumes that EGS drilling cost are added somewhere else to the model.')
    (
        pd.concat(regional_supplies)
        .dropna()
        .to_csv(snakemake.output['industrial_heating_egs_supply_curves'])
    )

    final_demands.to_csv(snakemake.output['industrial_heating_demands'])
    pd.concat(heat_exchanger_max_capacities, axis=1).to_csv(snakemake.output['heat_exchanger_capacity'])
