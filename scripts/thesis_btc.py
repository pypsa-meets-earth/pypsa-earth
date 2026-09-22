from __future__ import annotations


# =============================================================================
# Thesis BTC component names
# =============================================================================

BTC_LINK_NAME = "THESIS BTC mining link"
BTC_BUS_NAME = "THESIS BTC service bus"
BTC_STORE_NAME = "THESIS BTC service accumulator"

OLD_BTC_GENERATOR_NAME = "THESIS BTC flexible mining"


# =============================================================================
# Thesis-central BTC assumptions
# =============================================================================

DEFAULT_HASHPRICE_USD2026_PER_PH_DAY = 32.0

DEFAULT_ASIC_EFFICIENCY_J_PER_TH = 13.5

DEFAULT_PUE = 1.05

# ECB reference rate:
# 10 Aug 2026: 1 EUR = 1.1555 USD
DEFAULT_USD_PER_EUR = 1.1555

# Cumulative EUR2020 -> EUR2026 monetary factor
# derived from the thesis inflation chain
DEFAULT_EUR2026_TO_EUR2020_FACTOR = 1.29631

DEFAULT_NON_ELECTRIC_OPEX_EUR2020_PER_MWH = 0.0


# =============================================================================
# Economic conversion
# =============================================================================

def btc_hashprice_eur2020_per_th_day(
    hashprice_usd_per_ph_day: float
    = DEFAULT_HASHPRICE_USD2026_PER_PH_DAY,
    usd_per_eur: float = DEFAULT_USD_PER_EUR,
    eur2026_to_eur2020_factor: float
    = DEFAULT_EUR2026_TO_EUR2020_FACTOR,
) -> float:
    """
    Convert BTC hashprice from

        USD2026 / (PH/s) / day

    to

        EUR2020 / (TH/s) / day.

    Conversion chain
    ----------------
    1 PH/s = 1000 TH/s

    USD -> EUR:
        divide by USD per EUR

    EUR2026 -> EUR2020:
        divide by cumulative monetary factor
    """

    if hashprice_usd_per_ph_day <= 0:
        raise ValueError(
            "BTC hashprice must be positive."
        )

    if usd_per_eur <= 0:
        raise ValueError(
            "USD per EUR exchange rate must be positive."
        )

    if eur2026_to_eur2020_factor <= 0:
        raise ValueError(
            "EUR2026-to-EUR2020 factor must be positive."
        )

    return (
        hashprice_usd_per_ph_day
        / 1000.0
        / usd_per_eur
        / eur2026_to_eur2020_factor
    )


def btc_electricity_value_eur2020_per_mwh(
    hashprice_usd_per_ph_day: float
    = DEFAULT_HASHPRICE_USD2026_PER_PH_DAY,
    asic_efficiency_j_per_th: float
    = DEFAULT_ASIC_EFFICIENCY_J_PER_TH,
    pue: float = DEFAULT_PUE,
    usd_per_eur: float = DEFAULT_USD_PER_EUR,
    eur2026_to_eur2020_factor: float
    = DEFAULT_EUR2026_TO_EUR2020_FACTOR,
) -> float:
    """
    Calculate gross BTC mining revenue per MWh of
    facility-meter electricity.

    Returns
    -------
    float
        EUR2020 / MWh_facility

    Notes
    -----
    ASIC efficiency is expressed in J/TH.

    PUE converts ASIC electrical demand to total
    facility-meter electrical demand.

    The resulting value is a gross marginal revenue
    value, not an electricity price, LCOE, or complete
    mining profit.
    """

    if asic_efficiency_j_per_th <= 0:
        raise ValueError(
            "ASIC efficiency must be positive."
        )

    if pue < 1.0:
        raise ValueError(
            "PUE must be greater than or equal to 1."
        )

    hashprice_eur2020_per_th_day = (
        btc_hashprice_eur2020_per_th_day(
            hashprice_usd_per_ph_day=
                hashprice_usd_per_ph_day,
            usd_per_eur=usd_per_eur,
            eur2026_to_eur2020_factor=
                eur2026_to_eur2020_factor,
        )
    )

    # 1 MW = 1,000,000 J/s.
    #
    # ASIC efficiency:
    # J / TH
    #
    # Including PUE gives the hashrate supported
    # by 1 MW of total facility electricity.
    hashrate_th_per_s_per_mw = (
        1e6
        / (
            asic_efficiency_j_per_th
            * pue
        )
    )

    # Hashprice is daily revenue per TH/s.
    # Divide by 24 to obtain revenue per MWh
    # for a continuously operating 1 MW facility.
    value_eur2020_per_mwh = (
        hashrate_th_per_s_per_mw
        * hashprice_eur2020_per_th_day
        / 24.0
    )

    return value_eur2020_per_mwh


# =============================================================================
# PyPSA component
# =============================================================================

def add_btc_mining_link(
    n,
    electricity_bus: str = "KZ0 0",
    p_nom_mw: float = 1000.0,
    hashprice_usd_per_ph_day: float
    = DEFAULT_HASHPRICE_USD2026_PER_PH_DAY,
    asic_efficiency_j_per_th: float
    = DEFAULT_ASIC_EFFICIENCY_J_PER_TH,
    pue: float = DEFAULT_PUE,
    usd_per_eur: float = DEFAULT_USD_PER_EUR,
    eur2026_to_eur2020_factor: float
    = DEFAULT_EUR2026_TO_EUR2020_FACTOR,
    non_electric_opex_eur2020_per_mwh: float
    = DEFAULT_NON_ELECTRIC_OPEX_EUR2020_PER_MWH,
):
    """
    Add flexible Bitcoin mining as a PyPSA Link.

    Physical interpretation
    -----------------------
    Electricity is withdrawn from the Kazakhstan
    electricity bus and converted into a bookkeeping
    BTC-service carrier.

    The BTC-service carrier is not a physical energy
    carrier. It is only used to provide a balanced PyPSA
    representation of flexible electricity consumption.

    Dispatch
    --------
    BTC mining is fully flexible between 0 and p_nom.

    There is:
      - no annual BTC electricity target,
      - no deferred electricity demand,
      - no minimum utilization constraint.

    Economics
    ---------
    The Link marginal cost is negative because consuming
    electricity creates BTC mining revenue.

    Full mining CAPEX, ASIC replacement cost, and other
    fixed mining costs remain outside this national
    electricity-system optimization.
    """

    if electricity_bus not in n.buses.index:
        raise KeyError(
            f"BTC connection bus '{electricity_bus}' "
            "does not exist in the network."
        )

    if p_nom_mw <= 0:
        raise ValueError(
            "BTC electrical capacity must be positive."
        )

    if non_electric_opex_eur2020_per_mwh < 0:
        raise ValueError(
            "BTC non-electric marginal OPEX cannot "
            "be negative."
        )

    # -------------------------------------------------------------------------
    # Remove legacy BTC representation if present
    # -------------------------------------------------------------------------

    if OLD_BTC_GENERATOR_NAME in n.generators.index:
        n.remove(
            "Generator",
            OLD_BTC_GENERATOR_NAME,
        )

    # Allow safe rebuilding of the new representation.
    if BTC_STORE_NAME in n.stores.index:
        n.remove(
            "Store",
            BTC_STORE_NAME,
        )

    if BTC_LINK_NAME in n.links.index:
        n.remove(
            "Link",
            BTC_LINK_NAME,
        )

    # -------------------------------------------------------------------------
    # Economic calculation
    # -------------------------------------------------------------------------

    hashprice_eur2020_per_th_day = (
        btc_hashprice_eur2020_per_th_day(
            hashprice_usd_per_ph_day=
                hashprice_usd_per_ph_day,
            usd_per_eur=usd_per_eur,
            eur2026_to_eur2020_factor=
                eur2026_to_eur2020_factor,
        )
    )

    gross_value = (
        btc_electricity_value_eur2020_per_mwh(
            hashprice_usd_per_ph_day=
                hashprice_usd_per_ph_day,
            asic_efficiency_j_per_th=
                asic_efficiency_j_per_th,
            pue=pue,
            usd_per_eur=usd_per_eur,
            eur2026_to_eur2020_factor=
                eur2026_to_eur2020_factor,
        )
    )

    net_value = (
        gross_value
        - non_electric_opex_eur2020_per_mwh
    )

    if net_value <= 0:
        raise ValueError(
            "BTC net electricity value must be positive. "
            f"Calculated value: {net_value:.6f} "
            "EUR2020/MWh."
        )

    # -------------------------------------------------------------------------
    # Carriers
    # -------------------------------------------------------------------------

    if "bitcoin_mining" not in n.carriers.index:
        n.add(
            "Carrier",
            "bitcoin_mining",
        )

    if "bitcoin_service" not in n.carriers.index:
        n.add(
            "Carrier",
            "bitcoin_service",
        )

    # -------------------------------------------------------------------------
    # BTC bookkeeping bus
    # -------------------------------------------------------------------------

    if BTC_BUS_NAME not in n.buses.index:
        n.add(
            "Bus",
            BTC_BUS_NAME,
            carrier="bitcoin_service",
        )

    # -------------------------------------------------------------------------
    # BTC Link
    # -------------------------------------------------------------------------

    n.add(
        "Link",
        BTC_LINK_NAME,

        bus0=electricity_bus,
        bus1=BTC_BUS_NAME,

        carrier="bitcoin_mining",

        # Fixed facility-meter electrical capacity
        p_nom=p_nom_mw,
        p_nom_extendable=False,

        # Fully flexible load
        p_min_pu=0.0,
        p_max_pu=1.0,

        # Bookkeeping conversion only:
        # 1 MWh electricity ->
        # 1 MWh-equivalent BTC service
        efficiency=1.0,

        # Link p0 is positive electricity consumption.
        # Negative cost represents marginal BTC value.
        marginal_cost=-net_value,
    )

    # -------------------------------------------------------------------------
    # BTC-service accumulator
    # -------------------------------------------------------------------------
    #
    # This is NOT physical storage.
    #
    # It only absorbs Link output so that the BTC service
    # bus remains energy-balanced.
    #
    # The energy capacity equals the theoretical maximum
    # BTC service production over the modeled period.
    # -------------------------------------------------------------------------

    snapshot_hours = float(
        n.snapshot_weightings.stores.sum()
    )

    max_service_mwh = (
        p_nom_mw
        * snapshot_hours
    )

    n.add(
        "Store",
        BTC_STORE_NAME,

        bus=BTC_BUS_NAME,
        carrier="bitcoin_service",

        e_nom=max_service_mwh,
        e_nom_extendable=False,

        e_min_pu=0.0,
        e_max_pu=1.0,

        e_initial=0.0,
        e_cyclic=False,
        e_cyclic_per_period=False,

        standing_loss=0.0,

        marginal_cost=0.0,
        capital_cost=0.0,
    )

    # -------------------------------------------------------------------------
    # Reporting
    # -------------------------------------------------------------------------

    hashrate_th_per_s_per_mw = (
        1e6
        / (
            asic_efficiency_j_per_th
            * pue
        )
    )

    print(
        "\nBTC mining Link added"
        f"\n  Connection bus       : {electricity_bus}"
        f"\n  Electrical capacity  : {p_nom_mw:.3f} MW"
        f"\n  Hashprice             : "
        f"{hashprice_usd_per_ph_day:.3f} "
        "USD2026/(PH/s)/day"
        f"\n  Converted hashprice   : "
        f"{hashprice_eur2020_per_th_day:.9f} "
        "EUR2020/(TH/s)/day"
        f"\n  ASIC efficiency       : "
        f"{asic_efficiency_j_per_th:.3f} J/TH"
        f"\n  PUE                   : {pue:.4f}"
        f"\n  Hashrate per MW       : "
        f"{hashrate_th_per_s_per_mw:.3f} TH/s"
        f"\n  Gross electricity val.: "
        f"{gross_value:.6f} EUR2020/MWh"
        f"\n  Non-electric var OPEX : "
        f"{non_electric_opex_eur2020_per_mwh:.6f} "
        "EUR2020/MWh"
        f"\n  Net dispatch value    : "
        f"{net_value:.6f} EUR2020/MWh"
        f"\n  BTC service e_nom     : "
        f"{max_service_mwh / 1e6:.6f} TWh-equivalent\n"
    )

    return {
        "hashprice_usd2026_per_ph_day":
            hashprice_usd_per_ph_day,
        "hashprice_eur2020_per_th_day":
            hashprice_eur2020_per_th_day,
        "asic_efficiency_j_per_th":
            asic_efficiency_j_per_th,
        "pue":
            pue,
        "usd_per_eur":
            usd_per_eur,
        "eur2026_to_eur2020_factor":
            eur2026_to_eur2020_factor,
        "gross_value_eur2020_per_mwh":
            gross_value,
        "non_electric_opex_eur2020_per_mwh":
            non_electric_opex_eur2020_per_mwh,
        "net_value_eur2020_per_mwh":
            net_value,
        "hashrate_th_per_s_per_mw":
            hashrate_th_per_s_per_mw,
        "max_service_mwh":
            max_service_mwh,
    }
