# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import logging
import os
import pathlib
import subprocess
import sys

import country_converter as coco
import geopandas as gpd
import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# list of recognised nan values (NA and na excluded as may be confused with Namibia 2-letter country code)
NA_VALUES = ["NULL", "", "N/A", "NAN", "NaN", "nan", "Nan", "n/a", "null"]

REGION_COLS = ["geometry", "name", "x", "y", "country"]

# filename of the regions definition config file
REGIONS_CONFIG = "regions_definition_config.yaml"


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Customise errors traceback.
    """
    tb = exc_traceback
    while tb.tb_next:
        tb = tb.tb_next
    flname = tb.tb_frame.f_globals.get("__file__")
    funcname = tb.tb_frame.f_code.co_name

    if issubclass(exc_type, KeyboardInterrupt):
        logger.error(
            "Manual interruption %r, function %r: %s",
            flname,
            funcname,
            exc_value,
        )
    else:
        logger.error(
            "An error happened in module %r, function %r: %s",
            flname,
            funcname,
            exc_value,
            exc_info=(exc_type, exc_value, exc_traceback),
        )


def create_logger(logger_name, level=logging.INFO):
    """
    Create a logger for a module and adds a handler needed to capture in logs
    traceback from exceptions emerging during the workflow.
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    handler = logging.StreamHandler(stream=sys.stdout)
    logger.addHandler(handler)
    sys.excepthook = handle_exception
    return logger


def read_osm_config(*args):
    """
    Read values from the regions config file based on provided key arguments.

    Parameters
    ----------
    *args : str
        One or more key arguments corresponding to the values to retrieve
        from the config file. Typical arguments include "world_iso",
        "continent_regions", "iso_to_geofk_dict", and "osm_clean_columns".

    Returns
    -------
    tuple or str or dict
        If a single key is provided, returns the corresponding value from the
        regions config file. If multiple keys are provided, returns a tuple
        containing values corresponding to the provided keys.

    Examples
    --------
    >>> values = read_osm_config("key1", "key2")
    >>> print(values)
    ('value1', 'value2')

    >>> world_iso = read_osm_config("world_iso")
    >>> print(world_iso)
    {"Africa": {"DZ": "algeria", ...}, ...}
    """
    if "__file__" in globals():
        base_folder = get_dirname_path(__file__)
        if not pathlib.Path(base_folder, "configs").exists():
            base_folder = get_dirname_path(base_folder)
    else:
        base_folder = get_current_directory_path()
    osm_config_path = get_path(base_folder, "configs", REGIONS_CONFIG)
    with open(osm_config_path, "r") as f:
        osm_config = yaml.safe_load(f)
    if len(args) == 0:
        return osm_config
    elif len(args) == 1:
        return osm_config[args[0]]
    else:
        return tuple([osm_config[a] for a in args])


def sets_path_to_root(root_directory_name):
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
        if repo_name == get_basename_abs_path("."):
            repo_path = get_current_directory_path()  # current_path
            os.chdir(repo_path)  # change dir_path to repo_path
            print("This is the repository path: ", repo_path)
            print("Had to go %d folder(s) up." % (n0 - 1 - n))
            break
        # if repo_name NOT current folder name for 5 levels then stop
        if n == 0:
            print("Can't find the repo path.")
        # if repo_name NOT current folder name, go one directory higher
        else:
            change_to_script_dir(".")  # change to the upper folder


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """
    import logging

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = get_path(
            get_dirname_path(__file__), "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the "python" log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs, force=True)


def load_network(import_name=None, custom_components=None):
    """
    Helper for importing a pypsa.Network with additional custom components.

    Parameters
    ----------
    import_name : str
        As in pypsa.Network(import_name)
    custom_components : dict
        Dictionary listing custom components.
        For using ``snakemake.params.override_components"]``
        in ``config.yaml`` define:

        .. code:: yaml

            override_components:
                ShadowPrice:
                    component: ["shadow_prices","Shadow price for a global constraint.",np.nan]
                    attributes:
                    name: ["string","n/a","n/a","Unique name","Input (required)"]
                    value: ["float","n/a",0.,"shadow value","Output"]

    Returns
    -------
    pypsa.Network
    """
    import pypsa
    from pypsa.descriptors import Dict

    override_components = None
    override_component_attrs = None

    if custom_components is not None:
        override_components = pypsa.components.components.copy()
        override_component_attrs = Dict(
            {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
        )
        for k, v in custom_components.items():
            override_components.loc[k] = v["component"]
            override_component_attrs[k] = pd.DataFrame(
                columns=["type", "unit", "default", "description", "status"]
            )
            for attr, val in v["attributes"].items():
                override_component_attrs[k].loc[attr] = val

    return pypsa.Network(
        import_name=import_name,
        override_components=override_components,
        override_component_attrs=override_component_attrs,
    )


def pdbcast(v, h):
    return pd.DataFrame(
        v.values.reshape((-1, 1)) * h.values, index=v.index, columns=h.index
    )


def load_network_for_plots(
    fn, tech_costs, cost_config, elec_config, combine_hydro_ps=True
):
    import pypsa
    from add_electricity import load_costs, update_transmission_costs

    n = pypsa.Network(fn)

    n.loads["carrier"] = n.loads.bus.map(n.buses.carrier) + " load"
    n.stores["carrier"] = n.stores.bus.map(n.buses.carrier)

    n.links["carrier"] = (
        n.links.bus0.map(n.buses.carrier) + "-" + n.links.bus1.map(n.buses.carrier)
    )
    n.lines["carrier"] = "AC line"
    n.transformers["carrier"] = "AC transformer"

    n.lines["s_nom"] = n.lines["s_nom_min"]
    n.links["p_nom"] = n.links["p_nom_min"]

    if combine_hydro_ps:
        n.storage_units.loc[
            n.storage_units.carrier.isin({"PHS", "hydro"}), "carrier"
        ] = "hydro+PHS"

    # if the carrier was not set on the heat storage units
    # bus_carrier = n.storage_units.bus.map(n.buses.carrier)
    # n.storage_units.loc[bus_carrier == "heat","carrier"] = "water tanks"

    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = load_costs(tech_costs, cost_config, elec_config, Nyears)
    update_transmission_costs(n, costs)

    return n


def update_p_nom_max(n):
    # if extendable carriers (solar/onwind/...) have capacity >= 0,
    # e.g. existing assets from the OPSD project are included to the network,
    # the installed capacity might exceed the expansion limit.
    # Hence, we update the assumptions.

    n.generators.p_nom_max = n.generators[["p_nom_min", "p_nom_max"]].max(1)


def aggregate_p_nom(n):
    return pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum(),
            n.storage_units.groupby("carrier").p_nom_opt.sum(),
            n.links.groupby("carrier").p_nom_opt.sum(),
            n.loads_t.p.groupby(n.loads.carrier, axis=1).sum().mean(),
        ]
    )


def aggregate_p(n):
    return pd.concat(
        [
            n.generators_t.p.sum().groupby(n.generators.carrier).sum(),
            n.storage_units_t.p.sum().groupby(n.storage_units.carrier).sum(),
            n.stores_t.p.sum().groupby(n.stores.carrier).sum(),
            -n.loads_t.p.sum().groupby(n.loads.carrier).sum(),
        ]
    )


def aggregate_e_nom(n):
    return pd.concat(
        [
            (n.storage_units["p_nom_opt"] * n.storage_units["max_hours"])
            .groupby(n.storage_units["carrier"])
            .sum(),
            n.stores["e_nom_opt"].groupby(n.stores.carrier).sum(),
        ]
    )


def aggregate_p_curtailed(n):
    return pd.concat(
        [
            (
                (
                    n.generators_t.p_max_pu.sum().multiply(n.generators.p_nom_opt)
                    - n.generators_t.p.sum()
                )
                .groupby(n.generators.carrier)
                .sum()
            ),
            (
                (n.storage_units_t.inflow.sum() - n.storage_units_t.p.sum())
                .groupby(n.storage_units.carrier)
                .sum()
            ),
        ]
    )


def aggregate_costs(n, flatten=False, opts=None, existing_only=False):
    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
        Transformer=("s_nom", None),
    )

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components.keys(), skip_empty=False), components.values()
    ):
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs)

    if flatten:
        assert opts is not None
        conv_techs = opts["conv_techs"]

        costs = costs.reset_index(level=0, drop=True)
        costs = costs["capital"].add(
            costs["marginal"].rename({t: t + " marginal" for t in conv_techs}),
            fill_value=0.0,
        )

    return costs


def progress_retrieve(
    url, file, data=None, headers=None, disable_progress=False, roundto=1.0
):
    """
    Function to download data from a url with a progress bar progress in
    retrieving data.

    Parameters
    ----------
    url : str
        Url to download data from
    file : str
        File where to save the output
    data : dict
        Data for the request (default None), when not none Post method is used
    disable_progress : bool
        When true, no progress bar is shown
    roundto : float
        (default 0) Precision used to report the progress
        e.g. 0.1 stands for 88.1, 10 stands for 90, 80
    """
    import urllib

    from tqdm import tqdm

    pbar = tqdm(total=100, disable=disable_progress)

    def dlProgress(count, blockSize, totalSize, roundto=roundto):
        pbar.n = round(count * blockSize * 100 / totalSize / roundto) * roundto
        pbar.refresh()

    if data is not None:
        data = urllib.parse.urlencode(data).encode()

    if headers:
        opener = urllib.request.build_opener()
        opener.addheaders = headers
        urllib.request.install_opener(opener)

    urllib.request.urlretrieve(url, file, reporthook=dlProgress, data=data)


def get_aggregation_strategies(aggregation_strategies):
    """
    Default aggregation strategies that cannot be defined in .yaml format must
    be specified within the function, otherwise (when defaults are passed in
    the function's definition) they get lost when custom values are specified
    in the config.
    """
    import numpy as np

    # to handle the new version of PyPSA.
    try:
        from pypsa.clustering.spatial import _make_consense
    except Exception:
        # TODO: remove after new release and update minimum pypsa version
        from pypsa.clustering.spatial import _make_consense

    bus_strategies = dict(country=_make_consense("Bus", "country"))
    bus_strategies.update(aggregation_strategies.get("buses", {}))

    generator_strategies = {"build_year": lambda x: 0, "lifetime": lambda x: np.inf}
    generator_strategies.update(aggregation_strategies.get("generators", {}))

    return bus_strategies, generator_strategies


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the "scripts"-directory of "
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards**.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """

    import snakemake as sm
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = pathlib.Path(__file__).parent.resolve()
    assert (
        pathlib.Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if pathlib.Path(p).exists():
            snakefile = p
            break
    workflow = sm.Workflow(
        snakefile, overwrite_configfiles=[], rerun_triggers=[]
    )  # overwrite_config=config
    workflow.include(snakefile)
    workflow.global_resources = {}
    try:
        rule = workflow.get_rule(rulename)
    except Exception as exception:
        print(
            exception,
            f"The {rulename} might be a conditional rule in the Snakefile.\n"
            f"Did you enable {rulename} in the config?",
        )
        raise
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = get_abs_path(io[i])

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
    snakemake.benchmark = job.benchmark

    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        build_directory(path)

    os.chdir(script_dir)
    return snakemake


def two_2_three_digits_country(two_code_country):
    """
    Convert 2-digit to 3-digit country code:

    Parameters
    ----------
    two_code_country: str
        2-digit country name

    Returns
    ----------
    three_code_country: str
        3-digit country name
    """
    if two_code_country == "SN-GM":
        return f"{two_2_three_digits_country('SN')}-{two_2_three_digits_country('GM')}"

    three_code_country = coco.convert(two_code_country, to="ISO3")
    return three_code_country


def three_2_two_digits_country(three_code_country):
    """
    Convert 3-digit to 2-digit country code:

    Parameters
    ----------
    three_code_country: str
        3-digit country name

    Returns
    ----------
    two_code_country: str
        2-digit country name
    """
    if three_code_country == "SEN-GMB":
        return f"{three_2_two_digits_country('SN')}-{three_2_two_digits_country('GM')}"

    two_code_country = coco.convert(three_code_country, to="ISO2")
    return two_code_country


def two_digits_2_name_country(two_code_country, nocomma=False, remove_start_words=[]):
    """
    Convert 2-digit country code to full name country:

    Parameters
    ----------
    two_code_country: str
        2-digit country name
    nocomma: bool (optional, default False)
        When true, country names with comma are extended to remove the comma.
        Example CD -> Congo, The Democratic Republic of -> The Democratic Republic of Congo
    remove_start_words: list (optional, default empty)
        When a sentence starts with any of the provided words, the beginning is removed.
        e.g. The Democratic Republic of Congo -> Democratic Republic of Congo (remove_start_words=["The"])

    Returns
    ----------
    full_name: str
        full country name
    """
    if two_code_country == "SN-GM":
        return f"{two_digits_2_name_country('SN')}-{two_digits_2_name_country('GM')}"

    full_name = coco.convert(two_code_country, to="name_short")

    if nocomma:
        # separate list by delim
        splits = full_name.split(", ")

        # reverse the order
        splits.reverse()

        # return the merged string
        full_name = " ".join(splits)

    # when list is non empty
    if remove_start_words:
        # loop over every provided word
        for word in remove_start_words:
            # when the full_name starts with the desired word, then remove it
            if full_name.startswith(word):
                full_name = full_name.replace(word, "", 1)

    return full_name


def country_name_2_two_digits(country_name):
    """
    Convert full country name to 2-digit country code.

    Parameters
    ----------
    country_name: str
        country name

    Returns
    ----------
    two_code_country: str
        2-digit country name
    """
    if (
        country_name
        == f"{two_digits_2_name_country('SN')}-{two_digits_2_name_country('GM')}"
    ):
        return "SN-GM"

    full_name = coco.convert(country_name, to="ISO2")
    return full_name


def read_csv_nafix(file, **kwargs):
    "Function to open a csv as pandas file and standardize the na value"
    if "keep_default_na" not in kwargs:
        kwargs["keep_default_na"] = False
    if "na_values" not in kwargs:
        kwargs["na_values"] = NA_VALUES

    if get_path_size(file) > 0:
        return pd.read_csv(file, **kwargs)
    else:
        return pd.DataFrame()


def to_csv_nafix(df, path, **kwargs):
    if "na_rep" in kwargs:
        del kwargs["na_rep"]
    # if len(df) > 0:
    if not df.empty or not df.columns.empty:
        return df.to_csv(path, **kwargs, na_rep=NA_VALUES[0])
    else:
        with open(path, "w") as fp:
            pass


def save_to_geojson(df, fn):
    pathlib.Path(fn).unlink(missing_ok=True)  # remove file if it exists

    # save file if the (Geo)DataFrame is non-empty
    if df.empty:
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass
    else:
        # save file
        df.to_file(fn, driver="GeoJSON")


def read_geojson(fn, cols=[], dtype=None, crs="EPSG:4326"):
    """
    Function to read a geojson file fn. When the file is empty, then an empty
    GeoDataFrame is returned having columns cols, the specified crs and the
    columns specified by the dtype dictionary it not none.

    Parameters:
    ------------
    fn : str
        Path to the file to read
    cols : list
        List of columns of the GeoDataFrame
    dtype : dict
        Dictionary of the type of the object by column
    crs : str
        CRS of the GeoDataFrame
    """
    # if the file is non-zero, read the geodataframe and return it
    if get_path_size(fn) > 0:
        return gpd.read_file(fn)
    else:
        # else return an empty GeoDataFrame
        df = gpd.GeoDataFrame(columns=cols, geometry=[], crs=crs)
        if isinstance(dtype, dict):
            for k, v in dtype.items():
                df[k] = df[k].astype(v)
        return df


def create_country_list(input, iso_coding=True):
    """
    Create a country list for defined regions..

    Parameters
    ----------
    input : str
        Any two-letter country name, regional name, or continent given in the regions config file.
        Country name duplications won't distort the result.
        Examples are:
        ["NG","ZA"], downloading osm data for Nigeria and South Africa
        ["africa"], downloading data for Africa
        ["NAR"], downloading data for the North African Power Pool
        ["TEST"], downloading data for a customized test set.
        ["NG","ZA","NG"], won't distort result.

    Returns
    -------
    full_codes_list : list
        Example ["NG","ZA"]
    """
    import logging

    _logger = logging.getLogger(__name__)
    _logger.setLevel(logging.INFO)

    def filter_codes(c_list, iso_coding=True):
        """
        Filter list according to the specified coding.

        When iso code are implemented (iso_coding=True), then remove the
        geofabrik-specific ones. When geofabrik codes are
        selected(iso_coding=False), ignore iso-specific names.
        """
        if (
            iso_coding
        ):  # if country lists are in iso coding, then check if they are 2-string
            # 2-code countries
            ret_list = [c for c in c_list if len(c) == 2]

            # check if elements have been removed and return a working if so
            if len(ret_list) < len(c_list):
                _logger.warning(
                    "Specified country list contains the following non-iso codes: "
                    + ", ".join(list(set(c_list) - set(ret_list)))
                )

            return ret_list
        else:
            return c_list  # [c for c in c_list if c not in iso_to_geofk_dict]

    full_codes_list = []

    world_iso, continent_regions = read_osm_config("world_iso", "continent_regions")

    for value1 in input:
        codes_list = []
        # extract countries in world
        if value1 == "Earth":
            for continent in world_iso.keys():
                codes_list.extend(list(world_iso[continent]))

        # extract countries in continent
        elif value1 in world_iso.keys():
            codes_list = list(world_iso[value1])

        # extract countries in regions
        elif value1 in continent_regions.keys():
            codes_list = continent_regions[value1]

        # extract countries
        else:
            codes_list.extend([value1])

        # create a list with all countries
        full_codes_list.extend(codes_list)

    # Removing duplicates and filter outputs by coding
    full_codes_list = filter_codes(list(set(full_codes_list)), iso_coding=iso_coding)

    return full_codes_list


def get_last_commit_message(path):
    """
    Function to get the last PyPSA-Earth Git commit message.

    Returns
    -------
    result : string
    """
    _logger = logging.getLogger(__name__)
    last_commit_message = None
    backup_cwd = get_current_directory_path()
    try:
        os.chdir(path)
        last_commit_message = (
            subprocess.check_output(
                ["git", "log", "-n", "1", "--pretty=format:%H %s"],
                stderr=subprocess.STDOUT,
            )
            .decode()
            .strip()
        )
    except subprocess.CalledProcessError as e:
        _logger.warning(f"Error executing Git: {e}")

    os.chdir(backup_cwd)
    return last_commit_message


def get_dirname_path(input_path):
    """
    It returns the directory name of the path.
    """
    return pathlib.Path(input_path).parent


def get_abs_path(input_path):
    """
    It returns the absolutized version of the path.
    """
    return pathlib.Path(input_path).absolute()


def get_basename_abs_path(input_path):
    """
    It returns the base name of a normalized and absolutized version of the
    path.
    """
    return pathlib.Path(input_path).absolute().name


def get_basename_path(input_path):
    """
    It returns the base name of the path.
    """
    return pathlib.Path(input_path).name


def get_path(*args):
    """
    It returns a new posixPath object.
    """
    return pathlib.Path(*args)


def get_path_size(input_path):
    """
    It returns the size of a path (in bytes)
    """
    return pathlib.Path(input_path).stat().st_size


def build_directory(input_path):
    """
    It creates recursively the directory and its leaf directories.

    Parameters:
        input_path (str): The path to the file
    """

    # Check if the provided path points to a directory
    if is_directory_path(input_path):
        pathlib.Path(input_path).mkdir(parents=True, exist_ok=True)
    else:
        pathlib.Path(input_path).parent.mkdir(parents=True, exist_ok=True)


def change_to_script_dir(input_path):
    """
    Change the current working directory to the directory containing the given
    script.

    Parameters:
        input_path (str): The path to the file.
    """

    # Get the absolutized and normalized path of directory containing the file
    directory_path = pathlib.Path(input_path).absolute().parent

    # Change the current working directory to the script directory
    os.chdir(directory_path)


def get_current_directory_path():
    """
    It returns the current directory path.
    """
    return str(pathlib.Path.cwd())


def is_directory_path(input_path):
    """
    It returns True if the path points to a directory.

    False otherwise.
    """
    return pathlib.Path(input_path).is_dir()


def is_file_path(input_path):
    """
    It returns True if the path points to a file.

    False otherwise.
    """
    return pathlib.Path(input_path).is_file()


def get_relative_path(input_path, start_path="."):
    """
    It returns a relative path to input_path from start_path.

    Default for start_path is the current directory
    """
    return pathlib.Path(input_path).relative_to(start_path)
