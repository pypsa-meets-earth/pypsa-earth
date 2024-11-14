# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import io
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import time
import urllib

import country_converter as coco
import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import requests
import snakemake as sm
import yaml
from fake_useragent import UserAgent
from pypsa.clustering.spatial import _make_consense
from pypsa.components import component_attrs, components
from pypsa.descriptors import Dict
from shapely.geometry import Point
from snakemake.script import Snakemake
from tqdm import tqdm
from vresutils.costdata import annuity

logger = logging.getLogger(__name__)

# list of recognised nan values (NA and na excluded as may be confused with Namibia 2-letter country code)
NA_VALUES = ["NULL", "", "N/A", "NAN", "NaN", "nan", "Nan", "n/a", "null"]

REGION_COLS = ["geometry", "name", "x", "y", "country"]

# filename of the regions definition config file
REGIONS_CONFIG = "regions_definition_config.yaml"


def get_base_dir(file_path):
    return str(pathlib.Path(file_path).parent.parent.absolute())


def get_config_default_path(base_dir_path):
    return str(pathlib.Path(base_dir_path, "config.default.yaml"))


# prefix when running pypsa-earth rules in different directories (if running in pypsa-earth as subworkflow)
BASE_DIR = get_base_dir(__file__)

# absolute path to config.default.yaml
CONFIG_DEFAULT_PATH = get_config_default_path(BASE_DIR)


def check_config_version(config, fp_config=CONFIG_DEFAULT_PATH):
    """
    Check that a version of the local config.yaml matches to the actual config
    version as defined in config.default.yaml.
    """

    # using snakemake capabilities to deal with yanl configs
    with open(fp_config, "r") as f:
        actual_config = yaml.safe_load(f)
    actual_config_version = actual_config.get("version")

    current_config_version = config.get("version")

    if actual_config_version != current_config_version:
        logger.error(
            f"The current version of 'config.yaml' doesn't match to the code version:\n\r"
            f" {current_config_version} provided, {actual_config_version} expected.\n\r"
            f"That can lead to weird errors during execution of the workflow.\n\r"
            f"Please update 'config.yaml' according to 'config.default.yaml' ."
        )


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Customise errors traceback.
    """
    tb = exc_traceback
    while tb.tb_next:
        tb = tb.tb_next
    fl_name = tb.tb_frame.f_globals.get("__file__")
    func_name = tb.tb_frame.f_code.co_name

    if issubclass(exc_type, KeyboardInterrupt):
        logger.error(
            "Manual interruption %r, function %r: %s",
            fl_name,
            func_name,
            exc_value,
        )
    else:
        logger.error(
            "An error happened in module %r, function %r: %s",
            fl_name,
            func_name,
            exc_value,
            exc_info=(exc_type, exc_value, exc_traceback),
        )


def copy_default_files():
    fn = pathlib.Path(BASE_DIR, "config.yaml")
    if not fn.exists():
        fn.write_text(
            "# Write down config entries differing from config.default.yaml\n\nrun: {}"
        )


def create_logger(logger_name, level=logging.INFO):
    """
    Create a logger for a module and adds a handler needed to capture in logs
    traceback from exceptions emerging during the workflow.
    """
    logger_instance = logging.getLogger(logger_name)
    logger_instance.setLevel(level)
    handler = logging.StreamHandler(stream=sys.stdout)
    logger_instance.addHandler(handler)
    sys.excepthook = handle_exception
    return logger_instance


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
        base_folder = pathlib.Path(__file__).parent
        if not pathlib.Path(get_path(base_folder, "configs")).exists():
            base_folder = pathlib.Path(base_folder).parent
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

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = get_path(
            pathlib.Path(__file__).parent, "..", "logs", f"{snakemake.rule}.log"
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

    override_components = None
    override_component_attrs_dict = None

    if custom_components is not None:
        override_components = pypsa.components.components.copy()
        override_component_attrs_dict = Dict(
            {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
        )
        for k, v in custom_components.items():
            override_components.loc[k] = v["component"]
            override_component_attrs_dict[k] = pd.DataFrame(
                columns=["type", "unit", "default", "description", "status"]
            )
            for attr, val in v["attributes"].items():
                override_component_attrs_dict[k].loc[attr] = val

    return pypsa.Network(
        import_name=import_name,
        override_components=override_components,
        override_component_attrs=override_component_attrs_dict,
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
    """
    If extendable carriers (solar/onwind/...) have capacity >= 0, e.g. existing
    assets from the OPSD project are included to the network, the installed
    capacity might exceed the expansion limit.

    Hence, we update the assumptions.
    """
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
    components_dict = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
        Transformer=("s_nom", None),
    )

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components_dict.keys(), skip_empty=False),
        components_dict.values(),
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
    url, file, data=None, headers=None, disable_progress=False, round_to_value=1.0
):
    """
    Function to download data from an url with a progress bar progress in
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
    round_to_value : float
        (default 0) Precision used to report the progress
        e.g. 0.1 stands for 88.1, 10 stands for 90, 80
    """

    pbar = tqdm(total=100, disable=disable_progress)

    def dl_progress(count, block_size, total_size):
        pbar.n = (
            round(count * block_size * 100 / total_size / round_to_value)
            * round_to_value
        )
        pbar.refresh()

    if data is not None:
        data = urllib.parse.urlencode(data).encode()

    if headers:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req) as response:
            with open(file, "wb") as f:
                f.write(response.read())

    else:
        urllib.request.urlretrieve(url, file, reporthook=dl_progress, data=data)


def content_retrieve(url, data=None, headers=None, max_retries=3, backoff_factor=0.3):
    """
    Retrieve the content of a url with improved robustness.

    This function uses a more robust approach to handle permission issues
    and avoid being blocked by the server. It implements exponential backoff
    for retries and rotates user agents.

    Parameters
    ----------
    url : str
        URL to retrieve the content from
    data : dict, optional
        Data for the request, by default None
    headers : dict, optional
        Headers for the request, defaults to a fake user agent
        If no headers are wanted at all, pass an empty dict.
    max_retries : int, optional
        Maximum number of retries, by default 3
    backoff_factor : float, optional
        Factor to apply between attempts, by default 0.3
    """
    if headers is None:
        ua = UserAgent()
        headers = {
            "User-Agent": ua.random,
            "Upgrade-Insecure-Requests": "1",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
            "Accept-Encoding": "gzip, deflate, br",
            "DNT": "1",
            "Connection": "keep-alive",
            "Referer": "https://www.google.com/",
        }

    session = requests.Session()

    for i in range(max_retries):
        try:
            response = session.get(url, headers=headers, data=data)
            response.raise_for_status()
            return io.BytesIO(response.content)
        except (
            requests.exceptions.RequestException,
            requests.exceptions.HTTPError,
        ) as e:
            if i == max_retries - 1:  # last attempt
                raise
            else:
                # Exponential backoff
                wait_time = backoff_factor * (2**i) + np.random.uniform(0, 0.1)
                time.sleep(wait_time)

                # Rotate user agent for next attempt
                headers["User-Agent"] = UserAgent().random

    raise Exception("Max retries exceeded")


def get_aggregation_strategies(aggregation_strategies):
    """
    Default aggregation strategies that cannot be defined in .yaml format must
    be specified within the function, otherwise (when defaults are passed in
    the function's definition) they get lost when custom values are specified
    in the config.
    """

    bus_strategies = dict(country=_make_consense("Bus", "country"))
    bus_strategies.update(aggregation_strategies.get("buses", {}))

    generator_strategies = {"build_year": lambda x: 0, "lifetime": lambda x: np.inf}
    generator_strategies.update(aggregation_strategies.get("generators", {}))

    return bus_strategies, generator_strategies


def mock_snakemake(rule_name, root_dir=None, submodule_dir=None, config_file=None, **wildcards):
    """
    This function is expected to be executed from the "scripts"-directory of "
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards**.

    Parameters
    ----------
    rule_name: str
        name of the rule for which the snakemake object should be generated
    root_dir: str
        path to the root directory
    submodule_dir: str
        path to the submodule directory
    config_file: str
        path to config file to be used in mock_snakemake
    wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """

    script_dir = pathlib.Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = pathlib.Path(root_dir).resolve()

    user_in_script_dir = pathlib.Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(pathlib.Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif pathlib.Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in sm.SNAKEFILE_CHOICES:
            if pathlib.Path(p).exists():
                snakefile = p
                break

        if isinstance(config_file, str):
            with open(config_file, "r") as file:
                config_file = yaml.safe_load(file)

        workflow = sm.Workflow(
            snakefile,
            overwrite_configfiles=[],
            rerun_triggers=[],
            overwrite_config=config_file
        )
        workflow.include(snakefile)
        workflow.global_resources = {}
        try:
            rule = workflow.get_rule(rule_name)
        except Exception as exception:
            print(
                exception,
                f"The {rule_name} might be a conditional rule in the Snakefile.\n"
                f"Did you enable {rule_name} in the config?",
            )
            raise
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i in range(len(io)):
                    io[i] = pathlib.Path(io[i]).absolute()

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

    finally:
        if user_in_script_dir:
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


def two_digits_2_name_country(
    two_code_country, name_string="name_short", no_comma=False, remove_start_words=[]
):
    """
    Convert 2-digit country code to full name country:

    Parameters
    ----------
    two_code_country: str
        2-digit country name
    name_string: str (optional, default name_short)
        When name_short    CD -> DR Congo
        When name_official CD -> Democratic Republic of the Congo
    no_comma: bool (optional, default False)
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
    if remove_start_words is None:
        remove_start_words = list()
    if two_code_country == "SN-GM":
        return f"{two_digits_2_name_country('SN')}-{two_digits_2_name_country('GM')}"

    full_name = coco.convert(two_code_country, to=name_string)

    if no_comma:
        # separate list by delimiter
        splits = full_name.split(", ")

        # reverse the order
        splits.reverse()

        # return the merged string
        full_name = " ".join(splits)

    # when list is non-empty
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


def get_path(*args):
    """
    It returns a new path string.
    """
    return pathlib.Path(*args)


def get_path_size(path):
    """
    It returns the size of a path (in bytes)
    """
    return pathlib.Path(path).stat().st_size


def build_directory(path, just_parent_directory=True):
    """
    It creates recursively the directory and its leaf directories.

    Parameters:
        path (str): The path to the file
        just_parent_directory (Boolean): given a path dir/subdir
            True: it creates just the parent directory dir
            False: it creates the full directory tree dir/subdir
    """

    # Check if the provided path points to a directory
    if just_parent_directory:
        pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    else:
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def change_to_script_dir(path):
    """
    Change the current working directory to the directory containing the given
    script.

    Parameters:
        path (str): The path to the file.
    """

    # Get the absolutized and normalized path of directory containing the file
    directory_path = pathlib.Path(path).absolute().parent

    # Change the current working directory to the script directory
    os.chdir(directory_path)


def get_current_directory_path():
    """
    It returns the current directory path.
    """
    return pathlib.Path.cwd()


def get_relative_path(path, start_path="."):
    """
    It returns a relative path to path from start_path.

    Default for start_path is the current directory
    """
    return pathlib.Path(path).relative_to(start_path)


# PYPSA-EARTH-SEC


def prepare_costs(
    cost_file: str, USD_to_EUR: float, fill_values: dict, Nyears: float | int = 1
):
    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= USD_to_EUR

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = (
        costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    )
    costs = costs.fillna(fill_values)

    def annuity_factor(v):
        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    costs["fixed"] = [
        annuity_factor(v) * v["investment"] * Nyears for i, v in costs.iterrows()
    ]

    return costs


def create_network_topology(n, prefix, connector=" <-> ", bidirectional=True):
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

    # TODO: temporary fix for when underwater_fraction is not found
    if "underwater_fraction" not in n.links.columns:
        if n.links.empty:
            n.links["underwater_fraction"] = None
        else:
            n.links["underwater_fraction"] = 0.0

    candidates = pd.concat(
        [n.lines[ln_attrs], n.links.loc[n.links.carrier == "DC", lk_attrs]]
    ).fillna(0)

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
        topo = pd.concat([topo, topo_reverse])

    return topo


def create_dummy_data(n, sector):
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
    data = (
        np.random.randint(10, 500, size=(len(ind), len(col))) * 1000 * 1
    )  # TODO change 1 with temp. resolution

    return pd.DataFrame(data, index=ind, columns=col)


def cycling_shift(df, steps=1):
    """
    Cyclic shift on index of pd.Series|pd.DataFrame by number of steps.
    """
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
    Dictionary of overridden component attributes.
    """

    attrs = {k: v.copy() for k, v in component_attrs.items()}

    for component, list_name in components.list_name.items():
        fn = f"{directory}/{list_name}.csv"
        if pathlib.Path(fn).is_file():
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs


def get_gadm_filename(country_code, file_prefix="gadm41_"):
    """
    Function to get three digits country code for GADM.
    """
    special_codes_gadm = {
        "XK": "XKO",  # kosovo
        "CP": "XCL",  # clipperton island
        "SX": "MAF",  # saint-martin
        "TF": "ATF",  # french southern territories
        "AX": "ALA",  # aland
        "IO": "IOT",  # british indian ocean territory
        "CC": "CCK",  # cocos island
        "NF": "NFK",  # norfolk
        "PN": "PCN",  # pitcairn islands
        "JE": "JEY",  # jersey
        "XS": "XSP",  # spratly islands
        "GG": "GGY",  # guernsey
        "UM": "UMI",  # United States minor outlying islands
        "SJ": "SJM",  # svalbard
        "CX": "CXR",  # Christmas island
    }

    if country_code in special_codes_gadm:
        return file_prefix + special_codes_gadm[country_code]
    else:
        return file_prefix + two_2_three_digits_country(country_code)


def get_gadm_url(gadm_url_prefix, gadm_filename):
    """
    Function to get the gadm url given a gadm filename.
    """
    return gadm_url_prefix + gadm_filename + ".gpkg"


def download_gadm(
    country_code,
    file_prefix,
    gadm_url_prefix,
    gadm_input_file_args,
    update=False,
    out_logging=False,
):
    """
    Download gpkg file from GADM for a given country code.

    Parameters
    ----------
    country_code : str
        2-digit country name of the downloaded files
    file_prefix : str
        file prefix string
    gadm_url_prefix: str
        gadm url prefix
    gadm_input_file_args: list[str]
        gadm input file arguments list
    update : bool
        Update = true, forces re-download of files
    out_logging : bool
        out_logging = true, enables output logging

    Returns
    -------
    gpkg file per country
    """

    _logger = logging.getLogger(__name__)

    gadm_filename = get_gadm_filename(country_code, file_prefix)
    gadm_url = get_gadm_url(gadm_url_prefix, gadm_filename)
    gadm_input_file = get_path(
        get_current_directory_path(),
        *gadm_input_file_args,
        gadm_filename,
        gadm_filename,
    )

    gadm_input_file_gpkg = get_path(
        str(gadm_input_file) + ".gpkg"
    )  # Input filepath gpkg

    if not pathlib.Path(gadm_input_file_gpkg).exists() or update is True:
        if out_logging:
            _logger.warning(
                f"{gadm_filename} of country {two_digits_2_name_country(country_code)} does not exist, downloading to {gadm_input_file_gpkg}"
            )

        #  create data/osm directory
        build_directory(str(gadm_input_file_gpkg))

        try:
            r = requests.get(gadm_url, stream=True, timeout=300)
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
            raise Exception(
                f"GADM server is down at {gadm_url}. Data needed for building shapes can't be extracted.\n\r"
            )
        except Exception as exception:
            raise Exception(
                f"An error happened when trying to load GADM data by {gadm_url}.\n\r"
                + str(exception)
                + "\n\r"
            )
        else:
            with open(gadm_input_file_gpkg, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return gadm_input_file_gpkg, gadm_filename


def get_gadm_layer_name(file_prefix, layer_id):

    if file_prefix == "gadm41_":
        return "ADM_ADM_" + str(layer_id)
    else:
        raise Exception(
            f"The requested GADM data version {file_prefix} does not exist."
        )


def filter_gadm(
    geo_df,
    layer,
    cc,
    contended_flag,
    output_nonstd_to_csv=False,
):
    # identify non-standard geo_df rows
    geo_df_non_std = geo_df[geo_df["GID_0"] != two_2_three_digits_country(cc)].copy()

    if not geo_df_non_std.empty:
        logger.info(
            f"Contended areas have been found for gadm layer {layer}. They will be treated according to {contended_flag} option"
        )

        # NOTE: in these options GID_0 is not changed because it is modified below
        if contended_flag == "drop":
            geo_df.drop(geo_df_non_std.index, inplace=True)
        elif contended_flag != "set_by_country":
            # "set_by_country" option is the default; if this elif applies, the desired option falls back to the default
            logger.warning(
                f"Value '{contended_flag}' for option contented_flag is not recognized.\n"
                + "Fallback to 'set_by_country'"
            )

    # force GID_0 to be the country code for the relevant countries
    geo_df["GID_0"] = cc

    # country shape should have a single geometry
    if (layer == 0) and (geo_df.shape[0] > 1):
        logger.warning(
            f"Country shape is composed by multiple shapes that are being merged in agreement to contented_flag option '{contended_flag}'"
        )
        # take the first row only to re-define geometry keeping other columns
        geo_df = geo_df.iloc[[0]].set_geometry([geo_df.unary_union])

    # debug output to file
    if output_nonstd_to_csv and not geo_df_non_std.empty:
        geo_df_non_std.to_csv(
            f"resources/non_standard_gadm{layer}_{cc}_raw.csv", index=False
        )

    return geo_df


def get_gadm_layer(
    country_list,
    layer_id,
    geo_crs,
    file_prefix,
    gadm_url_prefix,
    gadm_input_file_args,
    contended_flag,
    update=False,
    out_logging=False,
):
    """
    Function to retrieve a specific layer id of a geopackage for a selection of
    countries.

    Parameters
    ----------
    country_list : str
        List of the countries
    layer_id : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested
    geo_crs: str
        General geographic projection
    file_prefix : str
        file prefix string
    gadm_url_prefix : str
        gadm url prefix
    gadm_input_file_args: list[str]
        gadm input file arguments list
    contended_flag : str
        contended areas
    update : bool
        Update = true, forces re-download of files
    out_logging : bool
        out_logging = true, enables output logging
    """
    # initialization of the list of geo dataframes
    geo_df_list = []

    for country_code in country_list:
        # download file gpkg
        file_gpkg, name_file = download_gadm(
            country_code,
            file_prefix,
            gadm_url_prefix,
            gadm_input_file_args,
            update,
            out_logging,
        )

        # get layers of a geopackage
        list_layers = fiona.listlayers(file_gpkg)

        # get layer name
        if layer_id < 0 or layer_id >= len(list_layers):
            # when layer id is negative or larger than the number of layers, select the last layer
            layer_id = len(list_layers) - 1
        code_layer = np.mod(layer_id, len(list_layers))
        layer_name = get_gadm_layer_name(file_prefix, layer_id)

        # read gpkg file
        geo_df_temp = gpd.read_file(
            file_gpkg, layer=layer_name, engine="pyogrio"
        ).to_crs(geo_crs)

        country_sub_index = ""
        if file_prefix == "gadm41_":
            country_sub_index = f"GID_{layer_id}"
            geo_df_temp = filter_gadm(
                geo_df=geo_df_temp,
                layer=layer_id,
                cc=country_code,
                contended_flag=contended_flag,
                output_nonstd_to_csv=False,
            )
        elif file_prefix == "gadm36_":
            country_sub_index = f"GID_{code_layer}"
            geo_df_temp["GID_0"] = [
                three_2_two_digits_country(twoD_c) for twoD_c in geo_df_temp["GID_0"]
            ]
        else:
            raise Exception(
                f"The requested GADM data version {file_prefix} does not exist."
            )

        geo_df_temp["GADM_ID"] = geo_df_temp[country_sub_index]

        # append geo data frames
        geo_df_list.append(geo_df_temp)

    geo_df_gadm = gpd.GeoDataFrame(pd.concat(geo_df_list, ignore_index=True))
    geo_df_gadm.set_crs(geo_crs, inplace=True)

    return geo_df_gadm


def locate_bus(
    coords,
    co,
    gadm_level,
    geo_crs,
    file_prefix,
    gadm_url_prefix,
    gadm_input_file_args,
    contended_flag,
    col_name="name",
    path_to_gadm=None,
    update=False,
    out_logging=False,
    gadm_clustering=False,
):
    """
    Function to locate the right node for a coordinate set input coords of
    point.

    Parameters
    ----------
    coords: pandas dataseries
        dataseries with 2 rows x & y representing the longitude and latitude
    co: string (code for country where coords are MA Morocco)
        code of the countries where the coordinates are
    gadm_level : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested
    geo_crs : str
        General geographic projection
    file_prefix : str
        file prefix string
    gadm_url_prefix: str
        gadm url prefix
    gadm_input_file_args: list[str]
        gadm input file arguments list
    contended_flag : str
        contended areas
    path_to_gadm : str
        path to gadm
    update : bool
        Update = true, forces re-download of files
    out_logging : bool
        out_logging = true, enables output logging
    gadm_clustering : bool
        gadm_cluster = true, to enable clustering
    col_name: str
        column to use to filter the GeoDataFrame
    """
    if not gadm_clustering:
        gdf = gpd.read_file(path_to_gadm)
    else:
        if path_to_gadm:
            gdf = gpd.read_file(path_to_gadm)
            if "GADM_ID" in gdf.columns:
                col = "GADM_ID"

                if gdf[col][0][
                    :3
                ].isalpha():  # TODO clean later by changing all codes to 2 letters
                    gdf[col] = gdf[col].apply(
                        lambda name: three_2_two_digits_country(name[:3]) + name[3:]
                    )
        else:
            gdf = get_gadm_layer(
                co,
                gadm_level,
                geo_crs,
                file_prefix,
                gadm_url_prefix,
                gadm_input_file_args,
                contended_flag,
                update,
                out_logging,
            )
            col_name = "GID_{}".format(gadm_level)

        # gdf.set_index("GADM_ID", inplace=True)
    gdf_co = gdf[
        gdf[col_name].str.contains(co)
    ]  # geodataframe of entire continent - output of prev function {} are placeholders
    # in strings - conditional formatting
    # insert any variable into that place using .format - extract string and filter for those containing co (MA)
    point = Point(coords["x"], coords["y"])  # point object

    try:
        return gdf_co[gdf_co.contains(point)][
            col_name
        ].item()  # filter gdf_co which contains point and returns the bus

    except ValueError:
        return gdf_co[gdf_co.geometry == min(gdf_co.geometry, key=(point.distance))][
            col_name
        ].item()  # looks for closest one shape=node


def get_conv_factors(sector):
    """
    Create a dictionary with all the conversion factors for the standard net calorific value
    from Tera Joule per Kilo Metric-ton to Tera Watt-hour based on
    https://unstats.un.org/unsd/energy/balance/2014/05.pdf.

    Considering that 1 Watt-hour = 3600 Joule, one obtains the values below dividing
    the standard net calorific values from the pdf by 3600.

    For example, the value "hard coal": 0.007167 is given by 25.8 / 3600, where 25.8 is the standard
    net calorific value.
    """

    conversion_factors_dict = {
        "additives and oxygenates": 0.008333,
        "anthracite": 0.005,
        "aviation gasoline": 0.01230,
        "bagasse": 0.002144,
        "biodiesel": 0.01022,
        "biogasoline": 0.007444,
        "bio jet kerosene": 0.011111,
        "bitumen": 0.01117,
        "brown coal": 0.003889,
        "brown coal briquettes": 0.00575,
        "charcoal": 0.00819,
        "coal tar": 0.007778,
        "coke-oven coke": 0.0078334,
        "coke-oven gas": 0.000277,
        "coking coal": 0.007833,
        "conventional crude oil": 0.01175,
        "crude petroleum": 0.011750,
        "ethane": 0.01289,
        "fuel oil": 0.01122,
        "fuelwood": 0.00254,
        "gas coke": 0.007326,
        "gas oil/ diesel oil": 0.01194,
        "gasoline-type jet fuel": 0.01230,
        "hard coal": 0.007167,
        "kerosene-type jet fuel": 0.01225,
        "lignite": 0.003889,
        "liquefied petroleum gas (lpg)": 0.01313,
        "lubricants": 0.011166,
        "motor gasoline": 0.01230,
        "naphtha": 0.01236,
        "natural gas": 0.00025,
        "natural gas liquids": 0.01228,
        "oil shale": 0.00247,
        "other bituminous coal": 0.005556,
        "other kerosene": 0.01216,
        "paraffin waxes": 0.01117,
        "patent fuel": 0.00575,
        "peat": 0.00271,
        "peat products": 0.00271,
        "petroleum coke": 0.009028,
        "refinery gas": 0.01375,
        "sub-bituminous coal": 0.005555,
    }

    if sector == "industry":
        return conversion_factors_dict
    else:
        logger.info(f"No conversion factors available for sector {sector}")
        return np.nan


def aggregate_fuels(sector):
    gas_fuels = [
        "blast furnace gas",
        "natural gas (including lng)",
    ]

    oil_fuels = [
        "additives and oxygenates",
        "aviation gasoline",
        "bitumen",
        "conventional crude oil",
        "crude petroleum",
        "ethane",
        "fuel oil",
        "gas oil/ diesel oil",
        "gasoline-type jet fuel",
        "kerosene-type jet fuel",
        "liquefied petroleum gas (lpg)",
        "lubricants",
        "motor gasoline",
        "naphtha",
        "natural gas liquids",
        "other kerosene",
        "paraffin waxes",
        "petroleum coke",
        "refinery gas",
    ]

    coal_fuels = [
        "anthracite",
        "blast furnace gas",
        "brown coal",
        "brown coal briquettes",
        "coal coke",
        "coal tar",
        "coke-oven coke",
        "coke-oven gas",
        "coking coal",
        "gas coke",
        "gasworks gas",
        "hard coal",
        "lignite",
        "oil shale",
        "other bituminous coal",
        "patent fuel",
        "peat",
        "peat products",
        "recovered gases",
        "sub-bituminous coal",
    ]

    biomass_fuels = [
        "bagasse",
        "fuelwood",
        "biogases",
        "biogasoline",
        "biodiesel",
        "charcoal",
        "black liquor",
        "bio jet kerosene",
        "animal waste",
        "industrial waste",
        "municipal wastes",
        "vegetal waste",
    ]

    electricity = ["electricity"]

    heat = ["heat", "direct use of geothermal heat", "direct use of solar thermal heat"]

    if sector == "industry":
        return gas_fuels, oil_fuels, biomass_fuels, coal_fuels, heat, electricity
    else:
        logger.info(f"No fuels available for sector {sector}")
        return np.nan


def modify_commodity(commodity):
    if commodity.strip() == "Hrad coal":
        commodity = "Hard coal"
    elif commodity.strip().casefold() == "coke oven gas":
        commodity = "Coke-oven gas"
    elif commodity.strip().casefold() == "coke oven coke":
        commodity = "Coke-oven coke"
    elif commodity.strip() == "Liquified Petroleum Gas (LPG)":
        commodity = "Liquefied Petroleum Gas (LPG)"
    elif commodity.strip() == "Gas Oil/Diesel Oil":
        commodity = "Gas Oil/ Diesel Oil"
    elif commodity.strip() == "Lignite brown coal- recoverable resources":
        commodity = "Lignite brown coal - recoverable resources"
    return commodity.strip().casefold()


def safe_divide(numerator, denominator, default_value=np.nan):
    """
    Safe division function that returns NaN when the denominator is zero.
    """
    if denominator != 0.0:
        return numerator / denominator
    else:
        logging.warning(
            f"Division by zero: {numerator} / {denominator}, returning NaN."
        )
        return default_value


def normed(x):
    return (x / x.sum()).fillna(0.0)
