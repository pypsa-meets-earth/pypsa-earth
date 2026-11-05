# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pandas as pd
import xarray as xr
from _helpers import configure_logging, create_logger
from add_electricity import load_powerplants

logger = create_logger(__name__)


def extract_inflow_df(
    ppl_df: pd.DataFrame,
    glofas_xr: xr.Dataset,
    # TODO Implement normalisation
    k: int = 1,
) -> pd.DataFrame:
    """
    Extract inflow for locations of hydropowerplants
    """
    glofas_copy_xr = glofas_xr.copy(deep=True)

    # TODO Account for the case when there is no hydro generation
    # NB 'technology' contains data on 'Reservoir' and 'Run-Of-River'
    ppl_hydro_df = ppl_df.query("carrier=='hydro'")

    ppl_hydro_lat = xr.DataArray(
        ppl_hydro_df["lat"].to_numpy(),
        dims="plant",
        coords={"plant": ppl_hydro_df.index},
    )

    ppl_hydro_lon = xr.DataArray(
        ppl_hydro_df["lon"].to_numpy(),
        dims="plant",
        coords={"plant": ppl_hydro_df.index},
    )
    # TODO Average by a few cells instead taking only the nearest one
    ppl_hydro_inflow_xr = glofas_copy_xr["dis24"].sel(
        latitude=ppl_hydro_lat,
        longitude=ppl_hydro_lon,
        method="nearest",
    )

    ppl_hydro_daily_inflow_df = ppl_hydro_inflow_xr.to_pandas()
    ppl_hydro_daily_inflow_df.index.name = "time"
    ppl_hydro_daily_inflow_df.columns.name = "plant"

    ppl_hydro_daily_inflow_df.index = pd.to_datetime(ppl_hydro_daily_inflow_df.index)
    ppl_hydro_inflow_df = ppl_hydro_daily_inflow_df.resample("1h").interpolate(
        method="linear"
    )

    return ppl_hydro_inflow_df


def transform_to_xr(inflow_df: pd.DataFrame) -> pd.DataFrame:
    """
    Transform dataframe into xarray dataset with structure
    of hydro renewable profile
    """

    tmp_df = inflow_df

    hydro_xr = xr.Dataset(
        data_vars={"inflow": (("plant", "time"), tmp_df.to_numpy().T)},
        coords={
            "plant": tmp_df.columns.to_numpy(),
            "time": tmp_df.index.to_numpy(),
        },
    )

    return hydro_xr


if __name__ == "__main__":

    # TODO Avoid excessive import
    from _helpers import mock_snakemake

    # if "snakemake" not in globals():
    #     from _helpers import mock_snakemake

    snakemake = mock_snakemake("build_glofas_profile")
    configure_logging(snakemake)

    ppls = load_powerplants(snakemake.input.powerplants)
    glofas_xr = xr.open_dataset(snakemake.input.glofas)

    inflow_ppl_df = extract_inflow_df(ppl_df=ppls, glofas_xr=glofas_xr)

    inflow_xr = transform_to_xr(inflow_ppl_df)

    inflow_xr.to_netcdf(snakemake.output.profile)
