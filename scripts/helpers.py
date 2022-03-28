from pathlib import Path
import os
import numpy as np
import pandas as pd
from vresutils.costdata import annuity
from pypsa.descriptors import Dict
from pypsa.components import components, component_attrs


def sets_path_to_root(root_directory_name):  # Imported from pypsa-africa
    """
    Search and sets path to the given root directory (root/path/file).

    Parameters
    ----------
    root_directory_name : str
        Name of the root directory.
    n : int
        Number of folders the function will check upwards/root directed.

    """
    import os

    repo_name = root_directory_name
    n = 8  # check max 8 levels above. Random default.
    n0 = n

    while n >= 0:
        n -= 1
        # if repo_name is current folder name, stop and set path
        if repo_name == os.path.basename(os.path.abspath(".")):
            repo_path = os.getcwd()  # os.getcwd() = current_path
            os.chdir(repo_path)  # change dir_path to repo_path
            print("This is the repository path: ", repo_path)
            print("Had to go %d folder(s) up." % (n0 - 1 - n))
            break
        # if repo_name NOT current folder name for 5 levels then stop
        if n == 0:
            print("Cant find the repo path.")
        # if repo_name NOT current folder name, go one dir higher
        else:
            upper_path = os.path.dirname(
                os.path.abspath("."))  # name of upper folder
            os.chdir(upper_path)


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    workflow = sm.Workflow(snakefile, overwrite_configfiles=[])
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake


def prepare_costs(cost_file, USD_to_EUR, discount_rate, Nyears, lifetime):

    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= USD_to_EUR

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = (costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(
        min_count=1))
    costs = costs.fillna({
        "CO2 intensity": 0,
        "FOM": 0,
        "VOM": 0,
        "discount rate": discount_rate,
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": lifetime,
    })

    def annuity_factor(v):
        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    costs["fixed"] = [
        annuity_factor(v) * v["investment"] * Nyears
        for i, v in costs.iterrows()
    ]

    return costs


def create_network_topology(n,
                            prefix,
                            like="ac",
                            connector=" <-> ",
                            bidirectional=True):
    """
    Create a network topology like the power transmission network.

    Parameters
    ----------
    n : pypsa.Network
    prefix : str
    connector : str
    bidirectional : bool, default True
        True: one link for each connection
        False: one link for each connection and direction (back and forth)

    Returns
    -------
    pd.DataFrame with columns bus0, bus1 and length
    """

    ln_attrs = ["bus0", "bus1", "length"]
    lk_attrs = ["bus0", "bus1", "length", "underwater_fraction"]

    candidates = pd.concat(
        [n.lines[ln_attrs], n.links.loc[n.links.carrier == "DC",
                                        lk_attrs]]).fillna(0)

    positive_order = candidates.bus0 < candidates.bus1
    candidates_p = candidates[positive_order]
    swap_buses = {"bus0": "bus1", "bus1": "bus0"}
    candidates_n = candidates[~positive_order].rename(columns=swap_buses)
    candidates = pd.concat([candidates_p, candidates_n])

    def make_index(c):
        return prefix + c.bus0 + connector + c.bus1

    topo = candidates.groupby(["bus0", "bus1"], as_index=False).mean()
    topo.index = topo.apply(make_index, axis=1)

    if not bidirectional:
        topo_reverse = topo.copy()
        topo_reverse.rename(columns=swap_buses, inplace=True)
        topo_reverse.index = topo_reverse.apply(make_index, axis=1)
        topo = topo.append(topo_reverse)

    return topo


def create_dummy_data(n, sector, carriers):
    ind = n.buses_t.p.index
    ind = n.buses.index[n.buses.carrier == "AC"]

    if sector == "industry":
        col = [
            "electricity",
            "coal",
            "coke",
            "solid biomass",
            "methane",
            "hydrogen",
            "low-temperature heat",
            "naphtha",
            "process emission",
            "process emission from feedstock",
            "current electricity",
        ]
    else:
        raise Exception("sector not found")
    data = np.random.randint(10, 500, size=(len(ind), len(col)))

    return pd.DataFrame(data, index=ind, columns=col)


def create_transport_data_dummy(pop_layout,
                                transport_data,
                                cars=4000000,
                                average_fuel_efficiency=0.7):

    for country in pop_layout.ct.unique():

        country_data = pd.DataFrame(
            data=[[cars, average_fuel_efficiency]],
            columns=transport_data.columns,
            index=[country],
        )
        transport_data = pd.concat([transport_data, country_data], axis=0)

    transport_data_dummy = transport_data

    return transport_data_dummy


def create_temperature_dummy(pop_layout, temperature):

    temperature_dummy = pd.DataFrame(index=temperature.index)

    for index in pop_layout.index:
        temperature_dummy[index] = temperature["ES0 0"]

    return temperature_dummy


def create_energy_totals_dummy(pop_layout, energy_totals):
    """
    Function to add additional countries specified in pop_layout.index to energy_totals, these countries take the same values as Spain
    """
    # All countries in pop_layout get the same values as Spain
    for country in pop_layout.ct.unique():
        energy_totals.loc[country] = energy_totals.loc["ES"]

    return energy_totals


def cycling_shift(df, steps=1):
    """Cyclic shift on index of pd.Series|pd.DataFrame by number of steps"""
    df = df.copy()
    new_index = np.roll(df.index, steps)
    df.values[:] = df.reindex(index=new_index).values
    return df


def override_component_attrs(directory):
    """Tell PyPSA that links can have multiple outputs by
    overriding the component_attrs. This can be done for
    as many buses as you need with format busi for i = 2,3,4,5,....
    See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs

    Parameters
    ----------
    directory : string
        Folder where component attributes to override are stored 
        analogous to ``pypsa/component_attrs``, e.g. `links.csv`.

    Returns
    -------
    Dictionary of overriden component attributes.
    """

    attrs = Dict({k: v.copy() for k, v in component_attrs.items()})

    for component, list_name in components.list_name.items():
        fn = f"{directory}/{list_name}.csv"
        if os.path.isfile(fn):
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs
