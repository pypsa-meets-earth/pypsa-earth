from pathlib import Path
import pandas as pd


OUT = Path("results/scenarios/thesis_tables")
OUT.mkdir(parents=True, exist_ok=True)


rows = [

    # ========================================================
    # SYSTEM DEFINITION
    # ========================================================

    {
        "Category": "System",
        "Parameter": "Study region",
        "Final value": "Kazakhstan",
        "Unit": "-",
        "Classification": "Scope definition",
        "Evidence / source": "Bachelor-thesis study scope",
        "Methodological role": (
            "Defines the geographical system boundary."
        ),
        "Important caveat": (
            "International power-system interactions are not "
            "represented as a full regional interconnected model."
        ),
    },

    {
        "Category": "System",
        "Parameter": "Model year",
        "Final value": "2045",
        "Unit": "year",
        "Classification": "Scenario definition",
        "Evidence / source": "Bachelor-thesis scenario design",
        "Methodological role": (
            "Future-year basis for demand, technology costs "
            "and system optimization."
        ),
        "Important caveat": (
            "2045 is the thesis study year. The zero-direct-CO2 "
            "electricity-system scenario is a thesis counterfactual "
            "and is not presented as Kazakhstan's statutory national target."
        ),
    },

    {
        "Category": "System",
        "Parameter": "Meteorological year",
        "Final value": "2013",
        "Unit": "year",
        "Classification": "Source/model input",
        "Evidence / source": "ERA5 / PyPSA-Earth weather dataset",
        "Methodological role": (
            "Provides hourly wind, solar and hydro-related "
            "meteorological conditions."
        ),
        "Important caveat": (
            "The 2045 system is optimized for one historical "
            "weather year."
        ),
    },

    {
        "Category": "System",
        "Parameter": "Temporal resolution",
        "Final value": "1",
        "Unit": "h",
        "Classification": "Model choice",
        "Evidence / source": "Final PyPSA-Earth configuration",
        "Methodological role": (
            "8760 hourly snapshots resolve renewable variability "
            "and flexible-load operation."
        ),
        "Important caveat": (
            "Sub-hourly flexibility and transient operation are "
            "outside the model scope."
        ),
    },

    {
        "Category": "System",
        "Parameter": "Spatial resolution",
        "Final value": "10",
        "Unit": "clusters",
        "Classification": "Model choice",
        "Evidence / source": "Final PyPSA-Earth configuration",
        "Methodological role": (
            "Balances regional representation and computational effort."
        ),
        "Important caveat": (
            "Local substation- and plant-level network constraints "
            "are aggregated."
        ),
    },

    # ========================================================
    # DEMAND
    # ========================================================

    {
        "Category": "Demand",
        "Parameter": "2045 annual electricity demand",
        "Final value": "187.0",
        "Unit": "TWh/a",
        "Classification": "Source-supported input",
        "Evidence / source": (
            "Kazakhstan long-term electricity-demand planning "
            "trajectory used in thesis source register"
        ),
        "Methodological role": (
            "Annual demand target; PyPSA-Earth hourly profile "
            "was normalized to this value."
        ),
        "Important caveat": (
            "Planning projection rather than a statutory demand target."
        ),
    },

    {
        "Category": "Demand",
        "Parameter": "Demand scaling factor",
        "Final value": "1.4338125546",
        "Unit": "-",
        "Classification": "Derived",
        "Evidence / source": (
            "187 TWh target divided by original "
            "130.421511 TWh PyPSA-Earth 2045 profile"
        ),
        "Methodological role": (
            "Preserves the hourly profile shape while matching "
            "the 2045 annual demand target."
        ),
        "Important caveat": (
            "No additional peak-shape modification was applied."
        ),
    },

    # ========================================================
    # GENERATION FLEET
    # ========================================================

    {
        "Category": "Generation fleet",
        "Parameter": "Existing conventional/hydro fleet",
        "Final value": "Kazakhstan-specific custom fleet",
        "Unit": "-",
        "Classification": "Source-supported dataset",
        "Evidence / source": (
            "Supervisor-provided pypsa-kz-data repository; "
            "custom_powerplants.csv"
        ),
        "Methodological role": (
            "Replaces the generic powerplantmatching fleet."
        ),
        "Important caveat": (
            "Used because the generic fleet contained material "
            "duplication/misclassification for Kazakhstan."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Existing coal capacity",
        "Final value": "12.967",
        "Unit": "GW",
        "Classification": "Source-supported dataset",
        "Evidence / source": "Kazakhstan-specific custom powerplant dataset",
        "Methodological role": "Existing fixed coal capacity.",
        "Important caveat": (
            "Capacity can remain installed in zero-direct-CO2 scenarios "
            "even when fossil dispatch is forced to zero."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Existing CCGT capacity",
        "Final value": "3.480",
        "Unit": "GW",
        "Classification": "Source-supported dataset",
        "Evidence / source": "Kazakhstan-specific custom powerplant dataset",
        "Methodological role": "Existing fixed CCGT capacity.",
        "Important caveat": (
            "Installed capacity is distinct from dispatch."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Existing OCGT capacity",
        "Final value": "1.6254",
        "Unit": "GW",
        "Classification": "Source-supported dataset",
        "Evidence / source": "Kazakhstan-specific custom powerplant dataset",
        "Methodological role": (
            "Existing gas-turbine capacity; expansion remains "
            "possible according to final model configuration."
        ),
        "Important caveat": (
            "The zero-direct-CO2 constraint prevents direct "
            "fossil-CO2-emitting generation in S1/S4/S5/S6."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Reservoir hydro power",
        "Final value": "2.66319",
        "Unit": "GW",
        "Classification": "Source-supported dataset",
        "Evidence / source": "Kazakhstan-specific custom powerplant dataset",
        "Methodological role": "Existing reservoir hydro fleet.",
        "Important caveat": (
            "Reservoir energy duration is separately simplified."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Reservoir duration",
        "Final value": "72",
        "Unit": "h",
        "Classification": "Explicit thesis assumption",
        "Evidence / source": (
            "Uniform modelling assumption applied to reservoir plants"
        ),
        "Methodological role": (
            "Represents simplified reservoir-energy availability."
        ),
        "Important caveat": (
            "Not a measured physical storage duration for each plant."
        ),
    },

    {
        "Category": "Generation fleet",
        "Parameter": "Implied reservoir energy capacity",
        "Final value": "191.74968",
        "Unit": "GWh",
        "Classification": "Derived",
        "Evidence / source": (
            "2.66319 GW reservoir capacity x 72 h"
        ),
        "Methodological role": (
            "Energy-equivalent representation of reservoir storage."
        ),
        "Important caveat": (
            "Derived from the uniform-duration assumption."
        ),
    },

    # ========================================================
    # RENEWABLES / GRID
    # ========================================================

    {
        "Category": "Renewables",
        "Parameter": "Solar and onshore-wind expansion",
        "Final value": "Endogenous",
        "Unit": "-",
        "Classification": "Model assumption",
        "Evidence / source": "PyPSA-Earth capacity-expansion formulation",
        "Methodological role": (
            "Optimizer determines economically optimal new capacity."
        ),
        "Important caveat": (
            "Available renewable potentials are very large and "
            "non-binding in the principal scenarios."
        ),
    },

    {
        "Category": "Transmission",
        "Parameter": "Transmission treatment",
        "Final value": "OSM topology; copt",
        "Unit": "-",
        "Classification": "PyPSA-Earth methodology",
        "Evidence / source": "Final PyPSA-Earth configuration",
        "Methodological role": (
            "Represents spatial transmission constraints and "
            "optimized line-capacity expansion."
        ),
        "Important caveat": (
            "AC line volume is a capacity-distance indicator, "
            "not transported electrical energy."
        ),
    },

    # ========================================================
    # ZERO DIRECT CO2
    # ========================================================

    {
        "Category": "Scenario",
        "Parameter": "Strict electricity-sector CO2 limit",
        "Final value": "0",
        "Unit": "MtCO2/a",
        "Classification": "Thesis scenario assumption",
        "Evidence / source": "S1/S4/S5/S6 scenario definition",
        "Methodological role": (
            "Prevents direct fossil-CO2-emitting generation."
        ),
        "Important caveat": (
            "Represents zero direct fossil CO2 in the modeled "
            "electricity system. It does not imply economy-wide or lifecycle "
            "carbon neutrality and is not Kazakhstan's official 2045 policy target."
        ),
    },

    # ========================================================
    # FLEXIBLE-CONSUMER LOCATION
    # ========================================================

    {
        "Category": "Flexible consumers",
        "Parameter": "Common connection cluster",
        "Final value": "KZ0 0",
        "Unit": "-",
        "Classification": "Model choice",
        "Evidence / source": (
            "10-cluster diagnostic: highest mean wind+solar "
            "capacity factor among candidate AC buses"
        ),
        "Methodological role": (
            "Common node for BTC and H2 prevents location from "
            "biasing their comparison."
        ),
        "Important caveat": (
            "Resource-oriented model cluster; not claimed as the "
            "uniquely optimal real-world project site."
        ),
    },

    # ========================================================
    # BITCOIN
    # ========================================================

    {
        "Category": "Bitcoin",
        "Parameter": "Mining electrical capacity",
        "Final value": "1.0",
        "Unit": "GW",
        "Classification": "Thesis scenario assumption",
        "Evidence / source": (
            "Controlled system-scale scenario; Kazakhstan mining "
            "scale used as plausibility context"
        ),
        "Methodological role": (
            "Maximum meter-side flexible mining demand."
        ),
        "Important caveat": (
            "Not a forecast of Kazakhstan mining capacity in 2045."
        ),
    },

    {
        "Category": "Bitcoin",
        "Parameter": "Mining flexibility",
        "Final value": "0-100",
        "Unit": "% of capacity",
        "Classification": "Source-supported model assumption",
        "Evidence / source": (
            "Bitcoin-mining demand-response literature and thesis "
            "source register"
        ),
        "Methodological role": (
            "Mining can be curtailed without a later electricity "
            "make-up requirement."
        ),
        "Important caveat": (
            "Operational detail such as startup delays is neglected "
            "at hourly resolution."
        ),
    },

    {
        "Category": "Bitcoin",
        "Parameter": "BTC gross electricity value",
        "Final value": "62.81",
        "Unit": "EUR2020/MWh",
        "Classification": "Derived market anchor",
        "Evidence / source": (
            "2026 hashprice anchor + ASIC efficiency + PUE + "
            "USD/EUR and inflation conversion"
        ),
        "Methodological role": (
            "Economic willingness-to-pay used for flexible BTC dispatch."
        ),
        "Important caveat": (
            "Dated 2026 gross-revenue anchor, not a 2045 BTC-price "
            "or profitability forecast; full mining CAPEX/OPEX is "
            "outside this PyPSA-Earth cost proxy."
        ),
    },

    # ========================================================
    # HYDROGEN
    # ========================================================

    {
        "Category": "Hydrogen",
        "Parameter": "PEM electrical capacity",
        "Final value": "1.0",
        "Unit": "GW_el",
        "Classification": "Thesis scenario assumption",
        "Evidence / source": (
            "Controlled comparison with 1 GW Bitcoin load; "
            "Kazakhstan project scale used as plausibility context"
        ),
        "Methodological role": (
            "Maximum electrical input of the PEM plant."
        ),
        "Important caveat": (
            "Not a forecast of installed Kazakhstan PEM capacity "
            "in 2045."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "Annual H2 production target",
        "Final value": "100",
        "Unit": "ktH2/a",
        "Classification": "Thesis scenario assumption",
        "Evidence / source": (
            "Scenario scaling informed by large Kazakhstan "
            "green-hydrogen project concepts"
        ),
        "Methodological role": (
            "Annual product requirement for S3/S5/S6."
        ),
        "Important caveat": (
            "Used as a controlled system-scale case, not a national "
            "production forecast."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "PEM LHV efficiency",
        "Final value": "0.640",
        "Unit": "-",
        "Classification": "Source-supported / derived",
        "Evidence / source": (
            "2045 PEM parameter provenance in thesis source register"
        ),
        "Methodological role": (
            "Converts AC electricity to hydrogen chemical energy."
        ),
        "Important caveat": (
            "Constant efficiency; part-load efficiency variation "
            "is neglected."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "Specific electricity consumption",
        "Final value": "52.078",
        "Unit": "kWh_el/kgH2",
        "Classification": "Derived",
        "Evidence / source": (
            "33.33 kWh_LHV/kgH2 divided by 0.640 efficiency"
        ),
        "Methodological role": (
            "Electrical energy required per kilogram of hydrogen."
        ),
        "Important caveat": (
            "Consistent with constant-LHV-efficiency representation."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "Annual PEM electricity requirement",
        "Final value": "5.207813",
        "Unit": "TWh_el/a",
        "Classification": "Derived",
        "Evidence / source": (
            "100 ktH2/a, 33.33 kWh_LHV/kg and eta_LHV=0.640"
        ),
        "Methodological role": (
            "Annual electricity requirement imposed indirectly "
            "through the H2 product target."
        ),
        "Important caveat": (
            "Hourly timing remains endogenous."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "PEM capacity factor",
        "Final value": "59.45",
        "Unit": "%",
        "Classification": "Derived",
        "Evidence / source": (
            "5.207813 TWh divided by 1 GW x 8760 h"
        ),
        "Methodological role": (
            "Annual utilization implied by capacity and production target."
        ),
        "Important caveat": (
            "The optimizer determines which specific hours are used."
        ),
    },

    {
        "Category": "Hydrogen",
        "Parameter": "PEM hourly operating range",
        "Final value": "0-100",
        "Unit": "% of capacity",
        "Classification": "Simplified model assumption",
        "Evidence / source": (
            "PEM flexibility literature interpreted at hourly resolution"
        ),
        "Methodological role": (
            "Allows production to shift toward low-marginal-value hours."
        ),
        "Important caveat": (
            "Minimum stable load, startup dynamics and degradation "
            "are not explicitly resolved."
        ),
    },

]


df = pd.DataFrame(rows)


# ============================================================
# Save CSV
# ============================================================

csv_path = (
    OUT
    / "table_04_parameter_source_assumption_provenance.csv"
)

df.to_csv(
    csv_path,
    index=False,
)


# ============================================================
# Save LaTeX longtable
# ============================================================

tex_path = (
    OUT
    / "table_04_parameter_source_assumption_provenance.tex"
)

latex = df.to_latex(
    index=False,
    escape=True,
    longtable=True,
    caption=(
        "Provenance and classification of the principal "
        "PyPSA-Earth model parameters used in the thesis."
    ),
    label="tab:model_parameter_provenance",
)

tex_path.write_text(
    latex,
    encoding="utf-8",
)


# ============================================================
# Classification summary
# ============================================================

summary = (
    df.groupby(
        "Classification"
    )
    .size()
    .rename("Number of parameters")
    .sort_values(
        ascending=False
    )
)

summary_path = (
    OUT
    / "table_04a_parameter_classification_summary.csv"
)

summary.to_csv(
    summary_path
)


# ============================================================
# Print
# ============================================================

print()
print("=" * 110)
print("MODEL PARAMETER PROVENANCE TABLE CREATED")
print("=" * 110)

print(
    df[
        [
            "Category",
            "Parameter",
            "Final value",
            "Unit",
            "Classification",
        ]
    ].to_string(
        index=False
    )
)

print()
print("Classification summary:")
print(summary.to_string())

print()
print("Written:")
print(csv_path)
print(tex_path)
print(summary_path)

print("=" * 110)
