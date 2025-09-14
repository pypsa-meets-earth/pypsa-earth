# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import calendar
import io
import logging
import os
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import yaml
from currency_converter import CurrencyConverter
from fake_useragent import UserAgent
from pypsa.components import component_attrs, components

logger = logging.getLogger(__name__)

currency_converter = CurrencyConverter(
    fallback_on_missing_rate=True,
    fallback_on_wrong_date=True,
)

# list of recognised nan values (NA and na excluded as may be confused with Namibia 2-letter country code)
NA_VALUES = ["NULL", "", "N/A", "NAN", "NaN", "nan", "Nan", "n/a", "null"]

REGION_COLS = ["geometry", "name", "x", "y", "country"]

# filename of the regions definition config file
REGIONS_CONFIG = "regions_definition_config.yaml"

# prefix when running pypsa-earth rules in different directories (if running in pypsa-earth as subworkflow)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# absolute path to config.default.yaml
CONFIG_DEFAULT_PATH = os.path.join(BASE_DIR, "config.default.yaml")


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
            f"Please update 'config.yaml' according to 'config.default.yaml.'\n\r"
            "If issues persist, consider to update the environment to the latest version."
        )


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


def copy_default_files():
    fn = Path(os.path.join(BASE_DIR, "config.yaml"))
    if not fn.exists():
        fn.write_text(
            "# Write down config entries differing from config.default.yaml\n\nrun: {}"
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
        base_folder = os.path.dirname(__file__)
        if not os.path.exists(os.path.join(base_folder, "configs")):
            base_folder = os.path.dirname(base_folder)
    else:
        base_folder = os.getcwd()
    osm_config_path = os.path.join(base_folder, "configs", REGIONS_CONFIG)
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
    import logging

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
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

    try:
        from pypsa.descriptors import Dict
    except:
        from pypsa.definitions.structures import Dict  # from pypsa version v0.31

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
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req) as response:
            with open(file, "wb") as f:
                f.write(response.read())

    else:
        urllib.request.urlretrieve(url, file, reporthook=dlProgress, data=data)


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


def mock_snakemake(
    rulename, root_dir=None, submodule_dir=None, configfile=None, **wildcards
):
    """
    This function is expected to be executed from the "scripts"-directory of "
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards**.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    configfile: str
        path to config file to be used in mock_snakemake
    wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm

    try:
        from pypsa.descriptors import Dict
    except:
        from pypsa.definitions.structures import Dict  # from pypsa version v0.31
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    user_in_script_dir = Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in sm.SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break

        if isinstance(configfile, str):
            with open(configfile, "r") as file:
                configfile = yaml.safe_load(file)

        workflow = sm.Workflow(
            snakefile,
            overwrite_configfiles=[],
            rerun_triggers=[],
            overwrite_config=configfile,
        )
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
        snakemake.benchmark = job.benchmark

        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

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

    if os.stat(file).st_size > 0:
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
    if os.path.exists(fn):
        os.unlink(fn)  # remove file if it exists

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
    if os.path.getsize(fn) > 0:
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
    backup_cwd = os.getcwd()
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


def update_config_dictionary(
    config_dict,
    parameter_key_to_fill="lines",
    dict_to_use={"geometry": "first", "bounds": "first"},
):
    config_dict.setdefault(parameter_key_to_fill, {})
    config_dict[parameter_key_to_fill].update(dict_to_use)
    return config_dict


# PYPSA-EARTH-SEC
def annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and.

    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
    """

    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


# Simple cache to avoid repeated computations and logging for same (currency, output_currency, year)
_currency_conversion_cache = {}


# Simple cache to avoid repeated computations and logging for same (currency, output_currency, year)
def get_yearly_currency_exchange_rate(
    initial_currency: str,
    output_currency: str,
    year: int,
    default_exchange_rate: float = None,
    _currency_conversion_cache: dict = None,
    future_exchange_rate_strategy: str = "latest",
    custom_future_exchange_rate: float = None,
):
    """
    Returns the average EUR-to-output currency exchange rate and the currency year.

    Uses cached values if available; otherwise computes the average from daily rates.
    Falls back to a default exchange rate if provided and no data is available.
    """

    if _currency_conversion_cache is None:
        _currency_conversion_cache = {}  # Use empty cache if not provided

    key = (initial_currency, output_currency, year)
    if key in _currency_conversion_cache:
        return _currency_conversion_cache[key]

    successful_years = []
    default_years = []

    # Helper inner function to handle one year at a time
    def _average_for_year(y):
        if initial_currency == "EUR":
            # EUR has no direct rates, use USD dates as reference
            available_dates = sorted(currency_converter._rates["USD"].keys())
        else:
            if initial_currency not in currency_converter._rates:
                if default_exchange_rate is not None:
                    default_years.append(y)
                    return default_exchange_rate
                raise RuntimeError(
                    f"No data for currency {initial_currency} and no default rate provided."
                )
            available_dates = sorted(currency_converter._rates[initial_currency].keys())

        max_date = available_dates[-1]
        # If year is beyond available data and strategy is "custom", return custom value
        if y > max_date.year:
            if future_exchange_rate_strategy == "custom":
                if custom_future_exchange_rate is not None:
                    logger.info(
                        f"Using custom future exchange rate ({custom_future_exchange_rate}) for {initial_currency}->{output_currency} in {y}."
                    )
                    return custom_future_exchange_rate
                else:
                    raise RuntimeError(
                        "Custom future exchange rate strategy selected, but no value was provided."
                    )
            # fallback to latest available year
            effective_year = max_date.year
            logger.info(
                f"Using latest available year ({effective_year}) for future exchange rate of {initial_currency}->{output_currency} in {y}."
            )
        else:
            effective_year = y
        dates_to_use = [d for d in available_dates if d.year == effective_year]

        rates = []
        for date in dates_to_use:
            try:
                rate = currency_converter.convert(
                    1, initial_currency, output_currency, date
                )
                rates.append(rate)
            except Exception:
                continue

        if rates:
            successful_years.append(effective_year)
            return sum(rates) / len(rates)

        if default_exchange_rate is not None:
            default_years.append(effective_year)
            return default_exchange_rate

        raise RuntimeError(
            f"No exchange rate data found for {initial_currency}->{output_currency} in {effective_year}, and no default rate provided."
        )

    avg_rate = _average_for_year(year)

    # Log only once per call, avoiding multiple repeated messages
    if successful_years:
        logger.info(f"Currency conversion succeeded for years: {successful_years}")
    if default_years:
        logger.warning(f"Using default exchange rate for years: {default_years}")

    # Save computed rate to cache
    _currency_conversion_cache[key] = avg_rate
    return avg_rate


def build_currency_conversion_cache(
    df,
    output_currency,
    default_exchange_rate=None,
    future_exchange_rate_strategy: str = "latest",
    custom_future_exchange_rate: float = None,
):
    """
    Builds a cache of exchange rates for all unique (output_currency, year) pairs in the dataset.

    Rates are computed once and stored for reuse to improve performance.
    """
    currency_list = currency_converter.currencies

    unique_keys = {
        (x["unit"][:3], output_currency, int(x["currency_year"]))
        for _, x in df.iterrows()
        if x["unit"][:3] in currency_list
    }

    _currency_conversion_cache = {}
    for key in unique_keys:
        initial_currency, _, year = key
        try:
            rate = get_yearly_currency_exchange_rate(
                initial_currency,
                output_currency,
                year,
                default_exchange_rate,
                _currency_conversion_cache=_currency_conversion_cache,
                future_exchange_rate_strategy=future_exchange_rate_strategy,
                custom_future_exchange_rate=custom_future_exchange_rate,
            )
            _currency_conversion_cache[key] = rate
        except Exception as e:
            logger.warning(f"Failed to get rate for {key}: {e}")
            continue

    return _currency_conversion_cache


def apply_currency_conversion(cost_dataframe, output_currency, cache):
    """
    Applies exchange rates from the cache to convert all cost values and units.

    Converts only rows with monetary units that start with a known currency symbol and contain '/'.
    """
    currency_list = currency_converter.currencies

    def convert_row(x):
        unit = x["unit"]
        value = x["value"]

        if not isinstance(unit, str) or "/" not in unit:
            return pd.Series([value, unit])

        currency = unit[:3]
        year = x.get("currency_year")

        if currency not in currency_list:
            return pd.Series([value, unit])

        if pd.isnull(year):
            logger.warning(
                f"Missing currency_year for row with unit '{unit}' and value '{value}'. Skipping currency conversion."
            )
            return pd.Series([value, unit])

        try:
            key = (currency, output_currency, int(year))
            rate = cache.get(key)
            if rate is not None:
                new_value = value * rate
                new_unit = unit.replace(currency, output_currency)
                return pd.Series([new_value, new_unit])
            else:
                logger.warning(f"Missing exchange rate for {key}. Skipping conversion.")
        except Exception as e:
            logger.warning(f"Failed to convert row {x.name}: {e}")

        return pd.Series([value, unit])

    cost_dataframe[["value", "unit"]] = cost_dataframe.apply(convert_row, axis=1)
    return cost_dataframe


def prepare_costs(
    cost_file: str,
    config: dict,
    output_currency: str,
    fill_values: dict,
    Nyears: float | int = 1,
    default_exchange_rate: float = None,
    future_exchange_rate_strategy: str = "latest",
    custom_future_exchange_rate: float = None,
):
    """
    Loads and processes cost data, converting units and currency to a common format.

    Applies currency conversion, fills missing values, and computes fixed annualized costs.
    """
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3

    # apply filter on financial_case and scenario, if they are contained in the cost dataframe
    wished_cost_scenario = config["cost_scenario"]
    wished_financial_case = config["financial_case"]
    for col in ["scenario", "financial_case"]:
        if col in costs.columns:
            costs[col] = costs[col].replace("", pd.NA)

    if "scenario" in costs.columns:
        costs = costs[
            (costs["scenario"].str.casefold() == wished_cost_scenario.casefold())
            | (costs["scenario"].isnull())
        ]

    if "financial_case" in costs.columns:
        costs = costs[
            (costs["financial_case"].str.casefold() == wished_financial_case.casefold())
            | (costs["financial_case"].isnull())
        ]

    if costs["currency_year"].isnull().any():
        logger.warning(
            "Some rows are missing 'currency_year' and will be skipped in currency conversion."
        )

    # Create a shared cache for exchange rates
    _currency_conversion_cache = build_currency_conversion_cache(
        costs,
        output_currency,
        default_exchange_rate,
        future_exchange_rate_strategy,
        custom_future_exchange_rate,
    )

    modified_costs = apply_currency_conversion(
        costs, output_currency, _currency_conversion_cache
    )

    modified_costs = (
        modified_costs.loc[:, "value"]
        .unstack(level=1)
        .groupby("technology")
        .sum(min_count=1)
    )
    modified_costs = modified_costs.fillna(fill_values)

    def annuity_factor(v):
        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    modified_costs["fixed"] = [
        annuity_factor(v) * v["investment"] * Nyears
        for _, v in modified_costs.iterrows()
    ]

    return modified_costs


def create_network_topology(
    n, prefix, like="ac", connector=" <-> ", bidirectional=True
):
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
    data = (
        np.random.randint(10, 500, size=(len(ind), len(col))) * 1000 * 1
    )  # TODO change 1 with temp. resolution

    return pd.DataFrame(data, index=ind, columns=col)


# def create_transport_data_dummy(pop_layout,
#                                 transport_data,
#                                 cars=4000000,
#                                 average_fuel_efficiency=0.7):

#     for country in pop_layout.ct.unique():

#         country_data = pd.DataFrame(
#             data=[[cars, average_fuel_efficiency]],
#             columns=transport_data.columns,
#             index=[country],
#         )
#         transport_data = pd.concat([transport_data, country_data], axis=0)

#     transport_data_dummy = transport_data

#     return transport_data_dummy

# def create_temperature_dummy(pop_layout, temperature):

#     temperature_dummy = pd.DataFrame(index=temperature.index)

#     for index in pop_layout.index:
#         temperature_dummy[index] = temperature["ES0 0"]

#     return temperature_dummy

# def create_energy_totals_dummy(pop_layout, energy_totals):
#     """
#     Function to add additional countries specified in pop_layout.index to energy_totals, these countries take the same values as Spain
#     """
#     # All countries in pop_layout get the same values as Spain
#     for country in pop_layout.ct.unique():
#         energy_totals.loc[country] = energy_totals.loc["ES"]

#     return energy_totals


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
        if os.path.isfile(fn):
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs


def get_country(target, **keys):
    """
    Function to convert country codes using pycountry.

    Parameters
    ----------
    target: str
        Desired type of country code.
        Examples:
        - 'alpha_3' for 3-digit
        - 'alpha_2' for 2-digit
        - 'name' for full country name
    keys: dict
        Specification of the country name and reference system.
        Examples:
        - alpha_3="ZAF" for 3-digit
        - alpha_2="ZA" for 2-digit
        - name="South Africa" for full country name

    Returns
    -------
    country code as requested in keys or np.nan, when country code is not recognized

    Example of usage
    -------
    - Convert 2-digit code to 3-digit codes: get_country('alpha_3', alpha_2="ZA")
    - Convert 3-digit code to 2-digit codes: get_country('alpha_2', alpha_3="ZAF")
    - Convert 2-digit code to full name: get_country('name', alpha_2="ZA")
    """
    import pycountry as pyc

    assert len(keys) == 1
    try:
        return getattr(pyc.countries.get(**keys), target)
    except (KeyError, AttributeError):
        return np.nan


def download_GADM(country_code, update=False, out_logging=False):
    """
    Download gpkg file from GADM for a given country code.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files
    update : bool
        Update = true, forces re-download of files

    Returns
    -------
    gpkg file per country
    """

    GADM_filename = f"gadm36_{two_2_three_digits_country(country_code)}"
    GADM_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/{GADM_filename}_gpkg.zip"
    _logger = logging.getLogger(__name__)
    GADM_inputfile_zip = os.path.join(
        os.getcwd(),
        "data",
        "raw",
        "gadm",
        GADM_filename,
        GADM_filename + ".zip",
    )  # Input filepath zip

    GADM_inputfile_gpkg = os.path.join(
        os.getcwd(),
        "data",
        "raw",
        "gadm",
        GADM_filename,
        GADM_filename + ".gpkg",
    )  # Input filepath gpkg

    if not os.path.exists(GADM_inputfile_gpkg) or update is True:
        if out_logging:
            _logger.warning(
                f"Stage 4/4: {GADM_filename} of country {two_digits_2_name_country(country_code)} does not exist, downloading to {GADM_inputfile_zip}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_zip), exist_ok=True)

        with requests.get(GADM_url, stream=True) as r:
            with open(GADM_inputfile_zip, "wb") as f:
                shutil.copyfileobj(r.raw, f)

        with zipfile.ZipFile(GADM_inputfile_zip, "r") as zip_ref:
            zip_ref.extractall(os.path.dirname(GADM_inputfile_zip))

    return GADM_inputfile_gpkg, GADM_filename


def _get_shape_col_gdf(path_to_gadm, co, gadm_layer_id, gadm_clustering):
    """
    Parameters
    ----------
    country_list : str
        List of the countries
    layer_id : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested
    """
    from build_shapes import get_GADM_layer

    col = "name"
    if not gadm_clustering:
        gdf_shapes = gpd.read_file(path_to_gadm)
    else:
        if path_to_gadm:
            gdf_shapes = gpd.read_file(path_to_gadm)
            if "GADM_ID" in gdf_shapes.columns:
                col = "GADM_ID"

                if gdf_shapes[col][0][
                    :3
                ].isalpha():  # TODO clean later by changing all codes to 2 letters
                    gdf_shapes[col] = gdf_shapes[col].apply(
                        lambda name: three_2_two_digits_country(name[:3]) + name[3:]
                    )
        else:
            gdf_shapes = get_GADM_layer(co, gadm_layer_id)
            col = "GID_{}".format(gadm_layer_id)
    gdf_shapes = gdf_shapes[gdf_shapes[col].str.contains(co)]
    return gdf_shapes, col


def locate_bus(
    df,
    countries,
    gadm_level,
    path_to_gadm=None,
    gadm_clustering=False,
    dropnull=True,
    col_out=None,
):
    """
    Function to locate the points of the dataframe df into the GADM shapefile.

    Parameters
    ----------
    df: pd.Dataframe
        Dataframe with mandatory x, y and country columns
    countries: list
        List of target countries
    gadm_level: int
        GADM level to be used
    path_to_gadm: str (default None)
        Path to the GADM shapefile
    gadm_clustering: bool (default False)
        True if gadm clustering is adopted
    dropnull: bool (default True)
        True if the rows with null values should be dropped
    col_out: str (default gadm_{gadm_level})
        Name of the output column
    """
    if col_out is None:
        col_out = "gadm_{}".format(gadm_level)
    df = df[df.country.isin(countries)]
    df[col_out] = None
    for co in countries:
        gdf_shape, col = _get_shape_col_gdf(
            path_to_gadm, co, gadm_level, gadm_clustering
        )
        sub_df = df.loc[df.country == co, ["x", "y", "country"]]
        gdf = gpd.GeoDataFrame(
            sub_df,
            geometry=gpd.points_from_xy(sub_df.x, sub_df.y),
            crs="EPSG:4326",
        )

        gdf_merged = gpd.sjoin_nearest(gdf, gdf_shape, how="inner", rsuffix="right")

        df.loc[gdf_merged.index, col_out] = gdf_merged[col]

    if dropnull:
        df = df[df[col_out].notnull()]

    return df


def get_conv_factors(sector):
    # Create a dictionary with all the conversion factors from ktons or m3 to TWh based on https://unstats.un.org/unsd/energy/balance/2014/05.pdf
    if sector == "industry":
        fuels_conv_toTWh = {
            "Gas Oil/ Diesel Oil": 0.01194,
            "Motor Gasoline": 0.01230,
            "Kerosene-type Jet Fuel": 0.01225,
            "Aviation gasoline": 0.01230,
            "Biodiesel": 0.01022,
            "Natural gas liquids": 0.01228,
            "Biogasoline": 0.007444,
            "Bitumen": 0.01117,
            "Fuel oil": 0.01122,
            "Liquefied petroleum gas (LPG)": 0.01313,
            "Liquified Petroleum Gas (LPG)": 0.01313,
            "Lubricants": 0.01117,
            "Naphtha": 0.01236,
            "Fuelwood": 0.00254,
            "Charcoal": 0.00819,
            "Patent fuel": 0.00575,
            "Brown coal briquettes": 0.00575,
            "Hard coal": 0.007167,
            "Hrad coal": 0.007167,
            "Other bituminous coal": 0.005556,
            "Anthracite": 0.005,
            "Peat": 0.00271,
            "Peat products": 0.00271,
            "Lignite": 0.003889,
            "Brown coal": 0.003889,
            "Sub-bituminous coal": 0.005555,
            "Coke-oven coke": 0.0078334,
            "Coke oven coke": 0.0078334,
            "Coke Oven Coke": 0.0078334,
            "Gasoline-type jet fuel": 0.01230,
            "Conventional crude oil": 0.01175,
            "Brown Coal Briquettes": 0.00575,
            "Refinery Gas": 0.01375,
            "Petroleum coke": 0.009028,
            "Coking coal": 0.007833,
            "Peat Products": 0.00271,
            "Petroleum Coke": 0.009028,
            "Additives and Oxygenates": 0.008333,
            "Bagasse": 0.002144,
            "Bio jet kerosene": 0.011111,
            "Crude petroleum": 0.011750,
            "Gas coke": 0.007326,
            "Gas Coke": 0.007326,
            "Refinery gas": 0.01375,
            "Coal Tar": 0.007778,
            "Paraffin waxes": 0.01117,
            "Ethane": 0.01289,
            "Oil shale": 0.00247,
            "Other kerosene": 0.01216,
        }
    return fuels_conv_toTWh


def aggregate_fuels(sector):
    gas_fuels = [
        "Natural gas (including LNG)",  #
        "Natural Gas (including LNG)",  #
    ]

    oil_fuels = [
        "Motor Gasoline",  ##
        "Liquefied petroleum gas (LPG)",  ##
        "Liquified Petroleum Gas (LPG)",  ##
        "Fuel oil",  ##
        "Kerosene-type Jet Fuel",  ##
        "Conventional crude oil",  #
        "Crude petroleum",  ##
        "Lubricants",
        "Naphtha",  ##
        "Gas Oil/ Diesel Oil",  ##
        "Petroleum coke",  ##
        "Petroleum Coke",  ##
        "Ethane",  ##
        "Bitumen",  ##
        "Refinery gas",  ##
        "Additives and Oxygenates",  #
        "Refinery Gas",  ##
        "Aviation gasoline",  ##
        "Gasoline-type jet fuel",  ##
        "Paraffin waxes",  ##
        "Natural gas liquids",  #
        "Other kerosene",
    ]

    biomass_fuels = [
        "Bagasse",  #
        "Fuelwood",  #
        "Biogases",
        "Biogasoline",  #
        "Biodiesel",  #
        "Charcoal",  #
        "Black Liquor",  #
        "Bio jet kerosene",  #
        "Animal waste",  #
        "Industrial Waste",  #
        "Industrial waste",
        "Municipal Wastes",  #
        "Vegetal waste",
    ]

    coal_fuels = [
        "Anthracite",
        "Brown coal",  #
        "Brown coal briquettes",  #
        "Coke oven coke",
        "Coke-oven coke",
        "Coke Oven Coke",
        "Coking coal",
        "Hard coal",  #
        "Hrad coal",  #
        "Other bituminous coal",
        "Sub-bituminous coal",
        "Coking coal",
        "Coke Oven Gas",  ##
        "Gas Coke",
        "Gasworks Gas",  ##
        "Lignite",  #
        "Peat",  #
        "Peat products",
        "Coal Tar",  ##
        "Brown Coal Briquettes",  ##
        "Gas coke",
        "Peat Products",
        "Oil shale",  #
        "Oil Shale",  #
        "Coal coke",  ##
        "Patent fuel",  ##
        "Blast Furnace Gas",  ##
        "Recovered gases",  ##
    ]

    electricity = ["Electricity"]

    heat = ["Heat", "Direct use of geothermal heat", "Direct use of solar thermal heat"]

    return gas_fuels, oil_fuels, biomass_fuels, coal_fuels, heat, electricity


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
        return pd.DataFrame(np.nan, index=numerator.index, columns=numerator.columns)


def lossy_bidirectional_links(n, carrier):
    """
    Split bidirectional links of type carrier into two unidirectional links to include transmission losses.
    """

    # identify all links of type carrier
    carrier_i = n.links.query("carrier == @carrier").index

    if carrier_i.empty:
        return

    logger.info(f"Splitting bidirectional links with the carrier {carrier}")

    # set original links to be unidirectional
    n.links.loc[carrier_i, "p_min_pu"] = 0

    # add a new links that mirror the original links, but represent the reversed flow direction
    # the new links have a cost and length of 0 to not distort the overall cost and network length
    rev_links = (
        n.links.loc[carrier_i].copy().rename({"bus0": "bus1", "bus1": "bus0"}, axis=1)
    )
    rev_links["length_original"] = rev_links[
        "length"
    ]  # tracker for the length of the original links length
    rev_links["capital_cost"] = 0
    rev_links["length"] = 0
    rev_links["reversed"] = True  # tracker for easy identification of reversed links
    rev_links.index = rev_links.index.map(lambda x: x + "-reversed")

    # add the new reversed links to the network and fill the newly created trackers with default values for the other links
    n.links = pd.concat([n.links, rev_links], sort=False)
    n.links["reversed"] = n.links["reversed"].fillna(False).infer_objects(copy=False)
    n.links["length_original"] = n.links["length_original"].fillna(n.links.length)


def set_length_based_efficiency(n, carrier, bus_suffix, transmission_efficiency):
    """
    Set the efficiency of all links of type carrier in network n based on their length and the values specified in the config.
    Additionally add the length based electricity demand required for compression (if applicable).
    The bus_suffix refers to the suffix that differentiates the links bus0 from the corresponding electricity bus, i.e. " H2".
    Important:
    Call this function AFTER lossy_bidirectional_links when creating links that are both bidirectional and lossy,
    and have a length based electricity demand for compression. Otherwise the compression will not consistently take place at
    the inflow bus and instead vary between the inflow and the outflow bus.
    """

    # get the links length based efficiency and required compression
    if carrier not in transmission_efficiency:
        raise KeyError(
            f"An error occurred when setting the length based efficiency for the Links of type {carrier}."
            f"The Link type {carrier} was not found in the config under config['sector']['transmission_efficiency']."
        )
    efficiencies = transmission_efficiency[carrier]
    efficiency_static = efficiencies.get("efficiency_static", 1)
    efficiency_per_1000km = efficiencies.get("efficiency_per_1000km", 1)
    compression_per_1000km = efficiencies.get("compression_per_1000km", 0)

    # indetify all links of type carrier
    carrier_i = n.links.loc[n.links.carrier == carrier].index

    # identify the lengths of all links of type carrier
    # use "length_original" for lossy bidirectional links and "length" for any other link
    if ("reversed" in n.links.columns) and any(n.links.loc[carrier_i, "reversed"]):
        lengths = n.links.loc[carrier_i, "length_original"]
    else:
        lengths = n.links.loc[carrier_i, "length"]

    # set the links' length based efficiency
    n.links.loc[carrier_i, "efficiency"] = (
        efficiency_static * efficiency_per_1000km ** (lengths / 1e3)
    )

    # set the links's electricity demand for compression
    if compression_per_1000km > 0:
        # connect the links to their corresponding electricity buses
        n.links.loc[carrier_i, "bus2"] = n.links.loc[
            carrier_i, "bus0"
        ].str.removesuffix(bus_suffix)
        # TODO: use these lines to set bus 2 instead, once n.buses.location is functional and remove bus_suffix.
        """
        n.links.loc[carrier_i, "bus2"] = n.links.loc[carrier_i, "bus0"].map(
            n.buses.location
        )  # electricity
        """
        # set the required compression demand
        n.links.loc[carrier_i, "efficiency2"] = -compression_per_1000km * lengths / 1e3


def branch(condition, then, otherwise=None):
    """
    This is a placeholder function that exists in Snakemake versions > 8.3.0.
    It can be removed once Snakemake is updated to a compatible version.
    """
    if condition:
        return then

    if otherwise is None:
        if isinstance(then, dict):
            return {}
        elif isinstance(then, str):
            return []
        else:
            return None

    return otherwise


def rename_techs(label):
    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch",
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        "H2": "hydrogen storage",
        "battery": "battery storage",
        "CCS": "CCS",
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines",
    }

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label


def add_missing_carriers(n, carriers):
    """
    Function to add missing carriers to the network without raising errors.
    """
    valid_carriers = {c for c in carriers if isinstance(c, str) and c.strip() != ""}
    missing_carriers = valid_carriers - set(n.carriers.index)
    if len(missing_carriers) > 0:
        for carrier in missing_carriers:
            n.add("Carrier", carrier)


def sanitize_carriers(n, config):
    """
    Sanitize the carrier information in a PyPSA Network object.

    The function ensures that all unique carrier names are present in the network's
    carriers attribute, and adds nice names and colors for each carrier according
    to the provided configuration dictionary.

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA Network object that represents an electrical power system.
    config : dict
        A dictionary containing configuration information, specifically the
        "plotting" key with "nice_names" and "tech_colors" keys for carriers.

    Returns
    -------
    None
        The function modifies the 'n' PyPSA Network object in-place, updating the
        carriers attribute with nice names and colors.

    Warnings
    --------
    Raises a warning if any carrier's "tech_colors" are not defined in the config dictionary.
    """

    for c in n.iterate_components():
        if "carrier" in c.df:
            add_missing_carriers(n, c.df.carrier)

    carrier_i = n.carriers.index
    nice_names = (
        pd.Series(config["plotting"]["nice_names"])
        .reindex(carrier_i)
        .fillna(carrier_i.to_series())
    )
    n.carriers["nice_name"] = n.carriers.nice_name.where(
        n.carriers.nice_name != "", nice_names
    )

    tech_colors = config["plotting"]["tech_colors"]
    colors = pd.Series(tech_colors).reindex(carrier_i)
    # try to fill missing colors with tech_colors after renaming
    missing_colors_i = colors[colors.isna()].index
    colors[missing_colors_i] = missing_colors_i.map(rename_techs).map(tech_colors)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(f"tech_colors for carriers {missing_i} not defined in config.")
    n.carriers["color"] = n.carriers.color.where(n.carriers.color != "", colors)


def sanitize_locations(n):
    if "location" in n.buses.columns:
        n.buses["x"] = n.buses.x.where(n.buses.x != 0, n.buses.location.map(n.buses.x))
        n.buses["y"] = n.buses.y.where(n.buses.y != 0, n.buses.location.map(n.buses.y))
        n.buses["country"] = n.buses.country.where(
            n.buses.country.ne("") & n.buses.country.notnull(),
            n.buses.location.map(n.buses.country),
        )
