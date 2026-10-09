from pathlib import Path

import pandas as pd

# ============================================================
# Paths
# ============================================================

TABLE_DIR = Path("results/scenarios/thesis_tables")

INPUT = TABLE_DIR / "table_04_parameter_source_assumption_provenance.csv"

OUTPUT_CSV = TABLE_DIR / "table_04_parameter_provenance_FINAL.csv"

OUTPUT_TEX = TABLE_DIR / "table_04_parameter_provenance_FINAL.tex"

SOURCE_CSV = TABLE_DIR / "table_04_source_bibliography_seed.csv"

STATUS_CSV = TABLE_DIR / "table_04_final_status_summary.csv"


if not INPUT.exists():
    raise FileNotFoundError(INPUT)


df = pd.read_csv(INPUT)


# ============================================================
# Final classification
#
# Exactly four categories:
#
# Verified source
# Derived
# Thesis assumption
# Model configuration
# ============================================================

STATUS = {
    "Study region": "Model configuration",
    "Model year": "Thesis assumption",
    "Meteorological year": "Model configuration",
    "Temporal resolution": "Model configuration",
    "Spatial resolution": "Model configuration",
    "2045 annual electricity demand": "Verified source",
    "Demand scaling factor": "Derived",
    "Existing conventional/hydro fleet": "Verified source",
    "Existing coal capacity": "Verified source",
    "Existing CCGT capacity": "Verified source",
    "Existing OCGT capacity": "Verified source",
    "Reservoir hydro power": "Verified source",
    "Reservoir duration": "Thesis assumption",
    "Implied reservoir energy capacity": "Derived",
    "Solar and onshore-wind expansion": "Model configuration",
    "Transmission treatment": "Model configuration",
    "Strict electricity-sector CO2 limit": "Thesis assumption",
    "Common connection cluster": "Thesis assumption",
    "Mining electrical capacity": "Thesis assumption",
    "Mining flexibility": "Thesis assumption",
    "BTC gross electricity value": "Derived",
    "PEM electrical capacity": "Thesis assumption",
    "Annual H2 production target": "Thesis assumption",
    "PEM LHV efficiency": "Derived",
    "Specific electricity consumption": "Derived",
    "Annual PEM electricity requirement": "Derived",
    "PEM capacity factor": "Derived",
    "PEM hourly operating range": "Thesis assumption",
}


# ============================================================
# Citation keys
#
# These are stable suggested BibTeX keys.
# Existing Zotero keys can be substituted later without
# changing the scientific classification.
# ============================================================

CITATIONS = {
    "Study region": "",
    "Model year": "",
    "Meteorological year": "CopernicusERA5",
    "Temporal resolution": "PyPSAEarth",
    "Spatial resolution": "PyPSAEarth",
    "2045 annual electricity demand": "KAZENERGY2023NPP; KZMinEnergy2035",
    "Demand scaling factor": "KAZENERGY2023NPP; Parzen2022GlobalDemand",
    "Existing conventional/hydro fleet": "PyPSAKZData2024",
    "Existing coal capacity": "PyPSAKZData2024",
    "Existing CCGT capacity": "PyPSAKZData2024",
    "Existing OCGT capacity": "PyPSAKZData2024",
    "Reservoir hydro power": "PyPSAKZData2024",
    "Reservoir duration": "",
    "Implied reservoir energy capacity": "PyPSAKZData2024",
    "Solar and onshore-wind expansion": "PyPSAEarth; IRENASTAT2023",
    "Transmission treatment": "PyPSAEarth; OpenStreetMap",
    "Strict electricity-sector CO2 limit": "",
    "Common connection cluster": "",
    "Mining electrical capacity": "KEGOC2021",
    "Mining flexibility": "Stoll2019; CambridgeMining2025",
    "BTC gross electricity value": (
        "BITMAIN2024S21XP; "
        "HashrateIndex2026Aug10; "
        "ECB2026FX; "
        "EurostatHICP; "
        "EC2026Forecast; "
        "Stoll2019"
    ),
    "PEM electrical capacity": "HyrasiaOne2022",
    "Annual H2 production target": "HyrasiaOne2022",
    "PEM LHV efficiency": "DanishEnergyAgencyRenewableFuels",
    "Specific electricity consumption": "DanishEnergyAgencyRenewableFuels",
    "Annual PEM electricity requirement": "DanishEnergyAgencyRenewableFuels",
    "PEM capacity factor": "DanishEnergyAgencyRenewableFuels",
    "PEM hourly operating range": "Tremel2018",
}


# ============================================================
# Exact primary URL / DOI
# ============================================================

URLS = {
    "Study region": "",
    "Model year": "",
    "Meteorological year": "https://doi.org/10.24381/cds.adbb2d47",
    "Temporal resolution": "https://pypsa-earth.readthedocs.io/en/v0.8.0/",
    "Spatial resolution": "https://pypsa-earth.readthedocs.io/en/v0.8.0/",
    "2045 annual electricity demand": (
        "https://www.kazenergy.com/upload/document/"
        "operation/forum/2023/"
        "Status%20and%20progress%20of%20the%20NPP%20"
        "programme%20in%20Kazakhstan.pdf"
    ),
    "Demand scaling factor": (
        "https://doi.org/10.5281/zenodo.6569890 | "
        "https://www.kazenergy.com/upload/document/"
        "operation/forum/2023/"
        "Status%20and%20progress%20of%20the%20NPP%20"
        "programme%20in%20Kazakhstan.pdf"
    ),
    "Existing conventional/hydro fleet": (
        "https://github.com/" "pypsa-meets-earth/pypsa-kz-data"
    ),
    "Existing coal capacity": ("https://github.com/" "pypsa-meets-earth/pypsa-kz-data"),
    "Existing CCGT capacity": ("https://github.com/" "pypsa-meets-earth/pypsa-kz-data"),
    "Existing OCGT capacity": ("https://github.com/" "pypsa-meets-earth/pypsa-kz-data"),
    "Reservoir hydro power": ("https://github.com/" "pypsa-meets-earth/pypsa-kz-data"),
    "Reservoir duration": "",
    "Implied reservoir energy capacity": (
        "https://github.com/" "pypsa-meets-earth/pypsa-kz-data"
    ),
    "Solar and onshore-wind expansion": (
        "https://pypsa-earth.readthedocs.io/en/v0.8.0/ | " "https://pxweb.irena.org/"
    ),
    "Transmission treatment": (
        "https://pypsa-earth.readthedocs.io/en/v0.8.0/ | "
        "https://www.openstreetmap.org/"
    ),
    "Strict electricity-sector CO2 limit": "",
    "Common connection cluster": "",
    "Mining electrical capacity": (
        "https://ar2021.kegoc.kz/pdf/" "AR2021_KEGOC_eng.pdf"
    ),
    "Mining flexibility": (
        "https://doi.org/10.1016/j.joule.2019.05.012 | "
        "https://www.jbs.cam.ac.uk/wp-content/uploads/"
        "2025/04/"
        "2025-04-cambridge-digital-mining-industry-report.pdf"
    ),
    "BTC gross electricity value": (
        "https://file12.bitmain.com/shop-product-s3/"
        "firmware/68414f17-e491-4879-a0df-6619a994dbeb/"
        "2024/07/29/16/"
        "S21%20XP%20Product%20Manual%20V1.0.8.pdf | "
        "https://hashrateindex.com/blog/"
        "hashrate-index-roundup-august-10-2026/ | "
        "https://www.ecb.europa.eu/stats/"
        "policy_and_exchange_rates/"
        "euro_reference_exchange_rates/html/"
        "eurofxref-graph-usd.en.html | "
        "https://ec.europa.eu/eurostat/ | "
        "https://economy-finance.ec.europa.eu/"
        "economic-forecast-and-surveys/economic-forecasts/"
        "spring-2026-economic-forecast-slowdown-growth-"
        "energy-shock-drives-inflation_en"
    ),
    "PEM electrical capacity": (
        "https://hyrasia.one/wp-content/uploads/"
        "2024/04/"
        "221027_Pressemitteilung-IA_HYRASIA-ONE.pdf"
    ),
    "Annual H2 production target": (
        "https://hyrasia.one/wp-content/uploads/"
        "2024/04/"
        "221027_Pressemitteilung-IA_HYRASIA-ONE.pdf"
    ),
    "PEM LHV efficiency": (
        "https://ens.dk/en/analyses-and-statistics/" "technology-data-renewable-fuels"
    ),
    "Specific electricity consumption": (
        "https://ens.dk/en/analyses-and-statistics/" "technology-data-renewable-fuels"
    ),
    "Annual PEM electricity requirement": (
        "https://ens.dk/en/analyses-and-statistics/" "technology-data-renewable-fuels"
    ),
    "PEM capacity factor": (
        "https://ens.dk/en/analyses-and-statistics/" "technology-data-renewable-fuels"
    ),
    "PEM hourly operating range": "https://doi.org/10.1007/978-3-319-72459-1",
}


# ============================================================
# Verification / methodological notes
# ============================================================

NOTES = {
    "Study region": (
        "Kazakhstan is the explicitly selected thesis " "system boundary."
    ),
    "Model year": (
        "2045 is the selected long-term thesis planning "
        "horizon and is not an official national "
        "decarbonization-target year."
    ),
    "Meteorological year": (
        "ERA5 is the verified meteorological data source; "
        "selection of 2013 is a model configuration choice. "
        "A single weather year remains a limitation."
    ),
    "Temporal resolution": ("Final scenarios use 8760 hourly snapshots."),
    "Spatial resolution": (
        "Electricity network aggregated to ten principal "
        "regions. Additional carrier buses mean the complete "
        "PyPSA network contains more than ten buses."
    ),
    "2045 annual electricity demand": (
        "VERIFIED: KAZENERGY-hosted 3 Oct 2023 presentation, "
        "slide 3, gives 187.0 billion kWh for 2045. "
        "The slide states 'Forecast balances up to 2035, "
        "vision up to 2050'. Treat 187 TWh as a planning/"
        "vision value, not a statutory demand target."
    ),
    "Demand scaling factor": (
        "Derived as 187.0 TWh divided by the original "
        "130.421511 TWh interpolated PyPSA-Earth/GEGIS "
        "2045 profile. Hourly shape is preserved."
    ),
    "Existing conventional/hydro fleet": (
        "Final workflow uses the supervisor-provided "
        "Kazakhstan-specific custom_powerplants.csv with "
        "custom_powerplants.method=replace."
    ),
    "Existing coal capacity": (
        "Technology total obtained from the Kazakhstan-"
        "specific custom powerplant dataset used by the "
        "final workflow."
    ),
    "Existing CCGT capacity": (
        "Technology total obtained from the Kazakhstan-"
        "specific custom powerplant dataset used by the "
        "final workflow."
    ),
    "Existing OCGT capacity": (
        "Technology total obtained from the Kazakhstan-"
        "specific custom powerplant dataset used by the "
        "final workflow."
    ),
    "Reservoir hydro power": (
        "Power capacity comes from the Kazakhstan-specific "
        "custom powerplant dataset."
    ),
    "Reservoir duration": (
        "Uniform 72 h duration is an explicit thesis "
        "simplification. It is not claimed to be measured "
        "physical storage duration for every reservoir."
    ),
    "Implied reservoir energy capacity": (
        "Derived as 2.66319 GW x 72 h = 191.74968 GWh."
    ),
    "Solar and onshore-wind expansion": (
        "Solar and onshore wind are endogenous expansion "
        "options. Existing/minimum renewable capacity is "
        "handled using PyPSA-Earth IRENA statistics."
    ),
    "Transmission treatment": (
        "OSM-based topology and copt treatment are workflow "
        "configuration choices retained consistently across "
        "the final scenarios."
    ),
    "Strict electricity-sector CO2 limit": (
        "0 MtCO2/a is an explicit thesis scenario assumption. "
        "It means zero direct fossil CO2 in the modeled "
        "electricity system, not economy-wide lifecycle "
        "carbon neutrality."
    ),
    "Common connection cluster": (
        "KZ0 0 was selected from the ten-cluster diagnostic "
        "because it had the highest mean wind+solar capacity "
        "factor. BTC and PEM use the same cluster to control "
        "for location. It is not claimed as a uniquely "
        "optimal real-world project site."
    ),
    "Mining electrical capacity": (
        "1 GW is a controlled thesis scenario scale. "
        "KEGOC reported digital-miner consumption above "
        "1,000 MW in 2021 excluding shadow mining, providing "
        "Kazakhstan-specific plausibility, not a 2045 forecast."
    ),
    "Mining flexibility": (
        "0-100% is a simplified controllable-load "
        "representation. Curtailed mining is lost computation "
        "and creates no deferred annual-energy obligation."
    ),
    "BTC gross electricity value": (
        "Derived from a rounded 32 USD/(PH/s day) Aug-2026 "
        "hashprice anchor, S21 XP 13.5 J/TH, PUE 1.05, "
        "USD/EUR conversion and conversion to EUR2020. "
        "62.81 EUR2020/MWh is a gross modeled electricity "
        "value/willingness-to-pay, not an electricity tariff "
        "and not a 2045 BTC-price forecast."
    ),
    "PEM electrical capacity": (
        "1 GW is a controlled system-scale thesis assumption "
        "chosen for comparison with the 1 GW BTC case. "
        "Hyrasia provides Kazakhstan-specific scale "
        "plausibility only."
    ),
    "Annual H2 production target": (
        "100 ktH2/a is a thesis scenario assumption. "
        "Hyrasia's planned 20 GW electrolysis / up to "
        "2 MtH2/a provides a Kazakhstan-specific scale "
        "benchmark, not a direct parameter transfer."
    ),
    "PEM LHV efficiency": (
        "2045 value is DERIVED by linear interpolation "
        "between the DEA PEMEC values used in the thesis "
        "provenance: 0.616 in 2040 and 0.664 in 2050, "
        "giving 0.640 in 2045."
    ),
    "Specific electricity consumption": (
        "Derived on an LHV basis: "
        "33.33 kWh_H2/kg / 0.640 = "
        "52.078125 kWh_el/kgH2."
    ),
    "Annual PEM electricity requirement": (
        "Derived from 100 million kgH2/a and "
        "52.078125 kWh_el/kgH2 = "
        "5.2078125 TWh_el/a."
    ),
    "PEM capacity factor": (
        "Derived from 5.2078125 TWh/a divided by "
        "1 GW x 8760 h = approximately 59.45%."
    ),
    "PEM hourly operating range": (
        "0-100% hourly dispatch is a deliberate simplified "
        "PEM-flexibility assumption. Detailed minimum-load, "
        "startup, degradation and sub-hourly dynamics are "
        "outside the principal PyPSA-Earth scenarios."
    ),
}


# ============================================================
# Strict checks
# ============================================================

parameters = set(df["Parameter"])

expected = set(STATUS)

missing_in_mapping = parameters - expected

missing_in_table = expected - parameters

if missing_in_mapping:
    raise RuntimeError(
        "Parameters found in table but not final mapping:\n"
        + "\n".join(sorted(missing_in_mapping))
    )

if missing_in_table:
    raise RuntimeError(
        "Mapped parameters missing from table:\n" + "\n".join(sorted(missing_in_table))
    )


# ============================================================
# Apply final provenance
# ============================================================

df["Final status"] = df["Parameter"].map(STATUS)

df["Citation key(s)"] = df["Parameter"].map(CITATIONS)

df["Primary URL / DOI"] = df["Parameter"].map(URLS)

df["Verification note"] = df["Parameter"].map(NOTES)


# Keep the previous Classification column for audit trail,
# but rename it so it cannot be confused with Final status.

if "Classification" in df.columns:
    df = df.rename(columns={"Classification": "Previous classification"})


# ============================================================
# Save full final CSV
# ============================================================

df.to_csv(
    OUTPUT_CSV,
    index=False,
)


# ============================================================
# Compact thesis LaTeX table
# ============================================================

compact = df[
    [
        "Category",
        "Parameter",
        "Final value",
        "Unit",
        "Final status",
        "Citation key(s)",
    ]
].copy()

compact.to_latex(
    OUTPUT_TEX,
    index=False,
    escape=True,
    longtable=True,
    caption=(
        "Final provenance classification of the principal "
        "PyPSA-Earth model parameters."
    ),
    label="tab:pypsa_parameter_provenance_final",
)


# ============================================================
# Source bibliography seed
# ============================================================

sources = pd.DataFrame(
    [
        {
            "Citation key": "KAZENERGY2023NPP",
            "Author / organisation": "KAZENERGY-hosted Kazakhstan Energy Week presentation",
            "Year": 2023,
            "Title": "Status and Progress of the NPP Programme in Kazakhstan",
            "Primary URL / DOI": (
                "https://www.kazenergy.com/upload/document/"
                "operation/forum/2023/"
                "Status%20and%20progress%20of%20the%20NPP%20"
                "programme%20in%20Kazakhstan.pdf"
            ),
            "Used for": "2045 electricity-demand planning/vision value",
            "Important location": "Slide 3",
        },
        {
            "Citation key": "KZMinEnergy2035",
            "Author / organisation": "Ministry of Energy of the Republic of Kazakhstan",
            "Year": 2022,
            "Title": "Energy balance of the Republic of Kazakhstan to 2035",
            "Primary URL / DOI": (
                "https://www.gov.kz/memleket/entities/"
                "energo/press/news/details/338994?lang=ru"
            ),
            "Used for": "Independent official cross-check of 2035 demand",
            "Important location": "152.9 billion kWh by 2035",
        },
        {
            "Citation key": "Parzen2022GlobalDemand",
            "Author / organisation": "Parzen, Franken and Fioriti",
            "Year": 2022,
            "Title": (
                "Global demand data for PyPSA-Earth: "
                "An Open Optimisation Model of the Earth Energy System"
            ),
            "Primary URL / DOI": "https://doi.org/10.5281/zenodo.6569890",
            "Used for": "PyPSA-Earth/GEGIS demand time-series basis",
            "Important location": "Zenodo version 0.0.3",
        },
        {
            "Citation key": "CopernicusERA5",
            "Author / organisation": "Copernicus Climate Change Service / ECMWF",
            "Year": "",
            "Title": "ERA5 hourly data on single levels",
            "Primary URL / DOI": "https://doi.org/10.24381/cds.adbb2d47",
            "Used for": "2013 weather data",
            "Important location": "",
        },
        {
            "Citation key": "PyPSAKZData2024",
            "Author / organisation": "PyPSA meets Earth / Agora Energiewende / OET",
            "Year": 2024,
            "Title": "pypsa-kz-data",
            "Primary URL / DOI": (
                "https://github.com/" "pypsa-meets-earth/pypsa-kz-data"
            ),
            "Used for": "Kazakhstan-specific powerplant dataset",
            "Important location": "data/custom_powerplants.csv",
        },
        {
            "Citation key": "PyPSAEarth",
            "Author / organisation": "PyPSA meets Earth",
            "Year": "",
            "Title": "PyPSA-Earth documentation",
            "Primary URL / DOI": "https://pypsa-earth.readthedocs.io/en/v0.8.0/",
            "Used for": "Workflow/model configuration",
            "Important location": "",
        },
        {
            "Citation key": "IRENASTAT2023",
            "Author / organisation": "International Renewable Energy Agency",
            "Year": 2023,
            "Title": "IRENASTAT renewable capacity statistics",
            "Primary URL / DOI": "https://pxweb.irena.org/",
            "Used for": "Existing renewable-capacity baseline",
            "Important location": "",
        },
        {
            "Citation key": "OpenStreetMap",
            "Author / organisation": "OpenStreetMap contributors",
            "Year": "",
            "Title": "OpenStreetMap",
            "Primary URL / DOI": "https://www.openstreetmap.org/",
            "Used for": "Transmission topology",
            "Important location": "",
        },
        {
            "Citation key": "KEGOC2021",
            "Author / organisation": "KEGOC JSC",
            "Year": 2021,
            "Title": "Annual Report 2021",
            "Primary URL / DOI": "https://ar2021.kegoc.kz/pdf/AR2021_KEGOC_eng.pdf",
            "Used for": "Kazakhstan BTC-mining scale plausibility",
            "Important location": "Digital miners >1000 MW excluding shadow mining",
        },
        {
            "Citation key": "Stoll2019",
            "Author / organisation": "Stoll, Klaassen and Gallersdorfer",
            "Year": 2019,
            "Title": "The Carbon Footprint of Bitcoin",
            "Primary URL / DOI": "https://doi.org/10.1016/j.joule.2019.05.012",
            "Used for": "Mining energy/PUE methodological support",
            "Important location": "",
        },
        {
            "Citation key": "CambridgeMining2025",
            "Author / organisation": "Cambridge Centre for Alternative Finance",
            "Year": 2025,
            "Title": "Cambridge Digital Mining Industry Report",
            "Primary URL / DOI": (
                "https://www.jbs.cam.ac.uk/wp-content/uploads/"
                "2025/04/"
                "2025-04-cambridge-digital-mining-industry-report.pdf"
            ),
            "Used for": "Modern Bitcoin-mining operating context",
            "Important location": "",
        },
        {
            "Citation key": "BITMAIN2024S21XP",
            "Author / organisation": "BITMAIN",
            "Year": 2024,
            "Title": "S21 XP Product Manual",
            "Primary URL / DOI": (
                "https://file12.bitmain.com/shop-product-s3/"
                "firmware/68414f17-e491-4879-a0df-6619a994dbeb/"
                "2024/07/29/16/"
                "S21%20XP%20Product%20Manual%20V1.0.8.pdf"
            ),
            "Used for": "270 TH/s, 3645 W, 13.5 J/TH",
            "Important location": "Specification table",
        },
        {
            "Citation key": "HashrateIndex2026Aug10",
            "Author / organisation": "Hashrate Index / Kaan Farahani",
            "Year": 2026,
            "Title": "Hashrate Index Roundup (August 10, 2026)",
            "Primary URL / DOI": (
                "https://hashrateindex.com/blog/"
                "hashrate-index-roundup-august-10-2026/"
            ),
            "Used for": "Dated BTC hashprice anchor",
            "Important location": (
                "Spot 31.73; 7-day avg 32.31; " "30-day avg 32.11 USD/PH/s/day"
            ),
        },
        {
            "Citation key": "ECB2026FX",
            "Author / organisation": "European Central Bank",
            "Year": 2026,
            "Title": "Euro foreign exchange reference rates",
            "Primary URL / DOI": (
                "https://www.ecb.europa.eu/stats/"
                "policy_and_exchange_rates/"
                "euro_reference_exchange_rates/html/"
                "eurofxref-graph-usd.en.html"
            ),
            "Used for": "EUR/USD conversion",
            "Important location": "10 August 2026",
        },
        {
            "Citation key": "EurostatHICP",
            "Author / organisation": "Eurostat",
            "Year": 2026,
            "Title": "HICP annual average inflation rates",
            "Primary URL / DOI": "https://ec.europa.eu/eurostat/",
            "Used for": "2021-2025 EUR inflation chain",
            "Important location": "Dataset prc_hicp_aind",
        },
        {
            "Citation key": "EC2026Forecast",
            "Author / organisation": "European Commission, DG ECFIN",
            "Year": 2026,
            "Title": "European Economic Forecast, Spring 2026",
            "Primary URL / DOI": "https://doi.org/10.2765/0071034",
            "Used for": "2026 EU inflation assumption",
            "Important location": "EU inflation 2026 = 3.1%",
        },
        {
            "Citation key": "HyrasiaOne2022",
            "Author / organisation": "HYRASIA ONE / Svevind Energy",
            "Year": 2022,
            "Title": "HYRASIA ONE Investment Agreement press release",
            "Primary URL / DOI": (
                "https://hyrasia.one/wp-content/uploads/"
                "2024/04/"
                "221027_Pressemitteilung-IA_HYRASIA-ONE.pdf"
            ),
            "Used for": "Kazakhstan H2 system-scale plausibility",
            "Important location": "40 GW RES; 20 GW electrolysis; up to 2 MtH2/a",
        },
        {
            "Citation key": "DanishEnergyAgencyRenewableFuels",
            "Author / organisation": "Danish Energy Agency",
            "Year": "",
            "Title": "Technology Data for Renewable Fuels",
            "Primary URL / DOI": (
                "https://ens.dk/en/analyses-and-statistics/"
                "technology-data-renewable-fuels"
            ),
            "Used for": "PEMEC technology parameter provenance",
            "Important location": "Hydrogen production via electrolysis",
        },
        {
            "Citation key": "Tremel2018",
            "Author / organisation": "Tremel",
            "Year": 2018,
            "Title": "Electricity-based Fuels",
            "Primary URL / DOI": "https://doi.org/10.1007/978-3-319-72459-1",
            "Used for": "PEM flexibility / system-integration support",
            "Important location": "",
        },
    ]
)


sources.to_csv(
    SOURCE_CSV,
    index=False,
)


# ============================================================
# Status summary
# ============================================================

summary = (
    df.groupby("Final status")
    .size()
    .rename("Number of parameters")
    .sort_values(ascending=False)
)

summary.to_csv(STATUS_CSV)


# ============================================================
# Print verification
# ============================================================

print()
print("=" * 120)
print("FINAL PYPSA-EARTH PARAMETER PROVENANCE")
print("=" * 120)

print(
    df[
        [
            "Category",
            "Parameter",
            "Final value",
            "Unit",
            "Final status",
            "Citation key(s)",
        ]
    ].to_string(index=False)
)

print()
print("=" * 120)
print("FINAL STATUS SUMMARY")
print("=" * 120)

print(summary.to_string())

print()
print("=" * 120)
print("187 TWh DEMAND SOURCE CHECK")
print("=" * 120)

demand = df.loc[df["Parameter"] == "2045 annual electricity demand"].iloc[0]

print(f"Value        : " f"{demand['Final value']} " f"{demand['Unit']}")

print(f"Status       : " f"{demand['Final status']}")

print(f"Citation     : " f"{demand['Citation key(s)']}")

print(f"Source       : " f"{demand['Primary URL / DOI']}")

print(f"Verification : " f"{demand['Verification note']}")

print()
print("Written:")
print(OUTPUT_CSV)
print(OUTPUT_TEX)
print(SOURCE_CSV)
print(STATUS_CSV)

print("=" * 120)
