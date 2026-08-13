# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Inspect cutout databundles and record their spatial and temporal coverage.

The script reads ``configs/bundle_config.yaml``, selects Zenodo-backed cutout
bundles, downloads each archive into a temporary workspace, extracts only
NetCDF cutout files, records their lon/lat extent and time range in one CSV file,
plots the spatial extents on a world map, and removes temporary data.

Examples
--------
Run all cutout bundles:

    python scripts/non_workflow/extract_cutout_bundle_info.py

Run only tutorial cutouts:

    python scripts/non_workflow/extract_cutout_bundle_info.py --tutorial-only

Create only the CSV:

    python scripts/non_workflow/extract_cutout_bundle_info.py --skip-plots

Create only the plots from an existing CSV:

    python scripts/non_workflow/extract_cutout_bundle_info.py --skip-csv
"""

import argparse
import logging
import os
import shutil
import tempfile
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlretrieve
from zipfile import ZipFile

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "matplotlib"),
)

import matplotlib

matplotlib.use("Agg")

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator

logger = logging.getLogger(__name__)

GROUPS = ("tutorial", "default")
X_COORDS = ("x", "lon", "longitude")
Y_COORDS = ("y", "lat", "latitude")
TIME_COORDS = ("time",)

CSV_COLUMNS = [
    "bundle",
    "group",
    "countries",
    "zenodo_url",
    "cutout_member",
    "time_start",
    "time_end",
    "x_min",
    "x_max",
    "y_min",
    "y_max",
]
PLOT_COLUMNS = ["bundle", "group", "x_min", "x_max", "y_min", "y_max"]
FAILURE_COLUMNS = ["bundle", "group", "zenodo_url", "error"]


@dataclass(frozen=True)
class CutoutBundle:
    """Zenodo-backed cutout bundle entry from ``bundle_config.yaml``."""

    name: str
    group: str
    countries: tuple[str, ...]
    url: str


def parse_args() -> argparse.Namespace:
    """Parse command-line options and validate temp-file handling."""

    parser = argparse.ArgumentParser(
        description=(
            "Download cutout bundles, collect metadata, and plot spatial coverage."
        )
    )
    parser.add_argument(
        "--config",
        default="configs/bundle_config.yaml",
        type=Path,
        help="Databundle configuration file.",
    )
    parser.add_argument(
        "--output-dir",
        default="resources/cutout_bundle_extents",
        type=Path,
        help="Directory where the metadata CSV and plots are written.",
    )
    parser.add_argument(
        "--csv-path",
        type=Path,
        help="CSV input/output path. Defaults to <output-dir>/cutout_bundle_info.csv.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        help="Optional temporary workspace. Defaults to a cleaned system temp dir.",
    )
    parser.add_argument(
        "--disable-progress",
        action="store_true",
        help="Hide download progress bars.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep downloaded zips and extracted NetCDFs in --work-dir.",
    )
    parser.add_argument(
        "--tutorial-only",
        action="store_true",
        help="Only process tutorial cutout bundles.",
    )
    parser.add_argument(
        "--skip-csv",
        action="store_true",
        help="Skip CSV creation and create plots from the existing CSV.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Create the CSV without creating plots.",
    )

    args = parser.parse_args()
    if args.skip_csv and args.skip_plots:
        parser.error("--skip-csv and --skip-plots cannot be used together.")
    if args.keep_temp and args.work_dir is None:
        parser.error("--keep-temp requires --work-dir.")

    return args


def load_cutout_bundles(config_path: Path) -> tuple[list[CutoutBundle], list[dict]]:
    """Return cutout bundles with Zenodo URLs and note unavailable ones."""

    with open(config_path) as f:
        config = yaml.safe_load(f)

    bundles: list[CutoutBundle] = []
    skipped: list[dict] = []

    for name, bundle_config in config.get("databundles", {}).items():
        if bundle_config.get("category") != "cutouts":
            continue

        group = "tutorial" if bundle_config.get("tutorial", False) else "default"
        url = bundle_config.get("urls", {}).get("zenodo")
        if not url:
            skipped.append(
                {
                    "bundle": name,
                    "group": group,
                    "reason": "cutout bundle has no Zenodo URL",
                }
            )
            continue

        bundles.append(
            CutoutBundle(
                name=name,
                group=group,
                countries=tuple(
                    str(country) for country in bundle_config.get("countries", [])
                ),
                url=url,
            )
        )

    return bundles, skipped


def selected_groups(tutorial_only: bool) -> tuple[str, ...]:
    """Return the group names requested by the CLI options."""

    return ("tutorial",) if tutorial_only else GROUPS


def filter_groups(
    bundles: list[CutoutBundle],
    skipped: list[dict],
    groups: tuple[str, ...],
) -> tuple[list[CutoutBundle], list[dict]]:
    """Filter bundle and skipped records to the requested group names."""

    return (
        [bundle for bundle in bundles if bundle.group in groups],
        [record for record in skipped if record["group"] in groups],
    )


def download(url: str, destination: Path, disable_progress: bool) -> None:
    """Download ``url`` to ``destination`` with an optional byte progress bar."""

    destination.parent.mkdir(parents=True, exist_ok=True)

    if disable_progress:
        urlretrieve(url, destination)
        return

    from tqdm import tqdm

    with tqdm(unit="B", unit_scale=True, unit_divisor=1024) as progress:
        downloaded = 0

        def report(blocks: int, block_size: int, total_size: int) -> None:
            nonlocal downloaded
            if total_size > 0:
                progress.total = total_size

            new_downloaded = blocks * block_size
            progress.update(new_downloaded - downloaded)
            downloaded = new_downloaded

        urlretrieve(url, destination, reporthook=report)


def zip_netcdfs(zip_path: Path) -> list[str]:
    """List NetCDF members in a downloaded bundle archive."""

    with ZipFile(zip_path) as archive:
        return [
            member.filename
            for member in archive.infolist()
            if not member.is_dir() and member.filename.endswith(".nc")
        ]


def extract_zip_member(zip_path: Path, member: str, output_dir: Path) -> Path:
    """Extract one archive member into ``output_dir`` using a safe flat name."""

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / member.replace("/", "__")

    with ZipFile(zip_path) as archive:
        with archive.open(member) as source, open(output_path, "wb") as target:
            shutil.copyfileobj(source, target)

    return output_path


def coordinate_name(dataset: xr.Dataset, candidates: tuple[str, ...]) -> str:
    """Find the first supported coordinate name in a NetCDF dataset."""

    for name in candidates:
        if name in dataset.coords or name in dataset.variables:
            return name

    raise ValueError(f"None of the coordinate names {candidates} were found.")


def coordinate_bounds(dataset: xr.Dataset, coord_name: str) -> tuple[float, float]:
    """Return finite min/max bounds for one coordinate."""

    values = np.asarray(dataset[coord_name].values, dtype=float)

    return float(np.nanmin(values)), float(np.nanmax(values))


def format_date(value: object) -> str:
    """Format a decoded NetCDF time value as DD-MM-YYYY."""

    if isinstance(value, np.datetime64):
        value = pd.Timestamp(value)
    if hasattr(value, "strftime"):
        return value.strftime("%d-%m-%Y")
    raise TypeError(f"Unsupported time value: {value!r}")


def time_bounds(dataset: xr.Dataset) -> tuple[str, str]:
    """Return date-formatted bounds for a cutout time coordinate."""

    time_coord = coordinate_name(dataset, TIME_COORDS)
    values = np.asarray(dataset[time_coord].values).reshape(-1)
    values = [value for value in values if not pd.isna(value)]
    if not values:
        raise ValueError("The time coordinate does not contain any valid values.")

    return format_date(min(values)), format_date(max(values))


def cutout_metadata(cutout_path: Path) -> dict:
    """Read spatial and temporal bounds from a NetCDF cutout file."""

    with xr.open_dataset(cutout_path) as dataset:
        x_coord = coordinate_name(dataset, X_COORDS)
        y_coord = coordinate_name(dataset, Y_COORDS)
        x_min, x_max = coordinate_bounds(dataset, x_coord)
        y_min, y_max = coordinate_bounds(dataset, y_coord)
        time_start, time_end = time_bounds(dataset)

    return {
        "time_start": time_start,
        "time_end": time_end,
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
    }


def process_bundle(
    bundle: CutoutBundle,
    work_dir: Path,
    disable_progress: bool,
    keep_temp: bool,
) -> list[dict]:
    """Download one bundle and return metadata records for its NetCDF cutouts."""

    bundle_dir = work_dir / bundle.name
    zip_path = bundle_dir / f"{bundle.name}.zip"
    extracted_dir = bundle_dir / "extracted"
    records: list[dict] = []

    try:
        logger.info("Downloading %s from %s", bundle.name, bundle.url)
        download(bundle.url, zip_path, disable_progress)

        members = zip_netcdfs(zip_path)
        if not members:
            raise ValueError("No NetCDF cutouts found in downloaded zip file.")

        for member in members:
            cutout_path = extract_zip_member(zip_path, member, extracted_dir)
            try:
                records.append(
                    {
                        "bundle": bundle.name,
                        "group": bundle.group,
                        "countries": ",".join(bundle.countries),
                        "zenodo_url": bundle.url,
                        "cutout_member": member,
                        **cutout_metadata(cutout_path),
                    }
                )
            finally:
                if not keep_temp:
                    with suppress(FileNotFoundError):
                        cutout_path.unlink()

        return records
    finally:
        if not keep_temp:
            shutil.rmtree(bundle_dir, ignore_errors=True)


def collect_metadata(
    bundles: list[CutoutBundle],
    work_dir: Path,
    disable_progress: bool,
    keep_temp: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process all selected bundles and collect successes and failures."""

    records: list[dict] = []
    failures: list[dict] = []

    for bundle in bundles:
        try:
            records.extend(
                process_bundle(bundle, work_dir, disable_progress, keep_temp)
            )
        except Exception as exc:
            logger.exception("Failed to process %s", bundle.name)
            failures.append(
                {
                    "bundle": bundle.name,
                    "group": bundle.group,
                    "zenodo_url": bundle.url,
                    "error": str(exc),
                }
            )

    return (
        pd.DataFrame(records, columns=CSV_COLUMNS),
        pd.DataFrame(failures, columns=FAILURE_COLUMNS),
    )


def bundle_label(bundle_name: str) -> str:
    """Convert a bundle identifier into a compact plot label."""

    return bundle_name


def add_world_map(ax) -> ccrs.PlateCarree:
    """Configure a full-world Equal Earth map with visible country borders."""

    data_crs = ccrs.PlateCarree()
    countries = cfeature.NaturalEarthFeature(
        category="cultural",
        name="admin_0_countries",
        scale="50m",
        facecolor="#f7f4ec",
        edgecolor="#7a7a7a",
        linewidth=0.35,
    )

    ax.set_global()
    ax.set_facecolor("#d8ecf5")
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#d8ecf5", zorder=0)
    ax.add_feature(countries, zorder=1)
    ax.coastlines(resolution="50m", color="#545454", linewidth=0.45, zorder=2)

    gridlines = ax.gridlines(
        crs=data_crs,
        draw_labels=False,
        linewidth=0.45,
        color="#9b9b9b",
        alpha=0.55,
        linestyle="--",
        zorder=3,
    )
    gridlines.xlocator = FixedLocator(np.arange(-180, 181, 60))
    gridlines.ylocator = FixedLocator(np.arange(-60, 61, 30))

    return data_crs


def plot_group_extents(extents: pd.DataFrame, group: str, output_path: Path) -> None:
    """Plot one group of cutout extents on a full-world Equal Earth map."""

    group_extents = extents[extents["group"] == group].copy()
    if group_extents.empty:
        logger.warning("No extents available for %s bundles; skipping plot.", group)
        return

    map_crs = ccrs.EqualEarth()
    fig, ax = plt.subplots(figsize=(15, 8.5), subplot_kw={"projection": map_crs})
    data_crs = add_world_map(ax)
    colors = plt.get_cmap("tab20", len(group_extents))

    for color_index, row in enumerate(group_extents.itertuples(index=False)):
        width = row.x_max - row.x_min
        height = row.y_max - row.y_min
        label = bundle_label(row.bundle)
        color = colors(color_index)

        ax.add_patch(
            Rectangle(
                (row.x_min, row.y_min),
                width,
                height,
                facecolor=color,
                edgecolor=color,
                alpha=0.28,
                linewidth=2.0,
                label=label,
                transform=data_crs,
                zorder=4,
            )
        )
        ax.text(
            row.x_min + width / 2,
            row.y_min + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=8,
            color="black",
            transform=data_crs,
            zorder=5,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
                "pad": 1.6,
            },
        )

    ax.set_title(f"{group.title()} cutout bundle extents")
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def write_metadata_csv(metadata: pd.DataFrame, csv_path: Path) -> None:
    """Write cutout metadata without replacing an existing CSV."""

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(csv_path, index=False, mode="x")


def create_plots_from_csv(
    csv_path: Path,
    output_dir: Path,
    groups: tuple[str, ...],
) -> None:
    """Read extent records from CSV and create one plot per selected group."""

    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Cannot create plots because the extent CSV does not exist: {csv_path}"
        )

    extents = pd.read_csv(csv_path)
    missing_columns = set(PLOT_COLUMNS).difference(extents.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Extent CSV is missing required columns: {missing}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for group in groups:
        plot_group_extents(
            extents,
            group,
            output_dir / f"cutout_bundle_extents_{group}.png",
        )
    logger.info("Wrote extent plots to %s", output_dir)


def collect_with_workspace(
    bundles: list[CutoutBundle],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collect metadata in either a user-provided or auto-cleaned workspace."""

    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=True)
        return collect_metadata(
            bundles,
            args.work_dir.resolve(),
            args.disable_progress,
            args.keep_temp,
        )

    with tempfile.TemporaryDirectory(prefix="cutout_bundle_extents_") as tmpdir:
        return collect_metadata(
            bundles,
            Path(tmpdir),
            args.disable_progress,
            args.keep_temp,
        )


def create_metadata_csv(
    args: argparse.Namespace,
    groups: tuple[str, ...],
    csv_path: Path,
) -> None:
    """Download selected bundles and write their metadata to CSV."""

    if csv_path.exists():
        raise FileExistsError(
            "Output file already exists. Refusing to overwrite existing metadata "
            f"CSV: {csv_path}"
        )

    bundles, skipped = load_cutout_bundles(args.config.resolve())
    bundles, skipped = filter_groups(bundles, skipped, groups)

    logger.info(
        "Found %s Zenodo cutout bundles for %s (%s tutorial, %s default).",
        len(bundles),
        ", ".join(groups),
        sum(bundle.group == "tutorial" for bundle in bundles),
        sum(bundle.group == "default" for bundle in bundles),
    )

    metadata, failures = collect_with_workspace(bundles, args)
    write_metadata_csv(metadata, csv_path)
    logger.info("Wrote cutout metadata CSV to %s", csv_path)

    if not failures.empty:
        logger.warning(
            "%s bundle(s) failed: %s",
            len(failures),
            ", ".join(failures["bundle"]),
        )
    if skipped:
        logger.info(
            "%s cutout bundle(s) skipped because they do not have Zenodo URLs: %s",
            len(skipped),
            ", ".join(record["bundle"] for record in skipped),
        )


def run(args: argparse.Namespace) -> None:
    """Create the metadata CSV and plots unless either step is skipped."""

    groups = selected_groups(args.tutorial_only)
    output_dir = args.output_dir.resolve()
    csv_path = (
        args.csv_path.resolve()
        if args.csv_path
        else output_dir / "cutout_bundle_info.csv"
    )

    if not args.skip_csv:
        create_metadata_csv(args, groups, csv_path)

    if not args.skip_plots:
        create_plots_from_csv(csv_path, output_dir, groups)


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(name)s:%(message)s",
    )
    run(args)


if __name__ == "__main__":
    main()
