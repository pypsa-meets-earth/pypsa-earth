# ASEAN Grid Infrastructure

This section details the modeling of the ASEAN electricity grid infrastructure within PyPSA-ASEAN, focusing on how national power development plans and the broader ASEAN Power Grid initiatives are integrated.

## National Power Development Plans and Prebuilt Networks

The ASEAN region is actively pursuing significant grid expansions extending into the 2040s. To accurately reflect future scenarios, these planned transmission lines must be incorporated into the model. While the standardization of national transmission plans across ASEAN countries varies, PyPSA-ASEAN addresses this by integrating available data.

**Key Data Sources for National Plans:**

The current methodology for incorporating national power development plans uses a piecemeal approach, drawing on various national sources for each country. (Further details on specific country sources would be listed here).

> **Note**: Current implementations assume that all identified transmission plans are already in place from the initial years of the simulation. This simplification is necessary because comprehensive data on the exact commissioning year for all planned transmissions is often unavailable. Additionally, managing disconnected grids and nodes can complicate the clustering process, making it challenging to simplify large-scale models efficiently.

**Leveraging Prebuilt Networks:**

A distinguishing feature in PyPSA-ASEAN (and PyPSA-Eur, unlike PyPSA-Earth) is the utilization of prebuilt networks. Generating, processing, and merging data from OpenStreetMap (OSM) is a time-consuming process. To streamline this, PyPSA-ASEAN provides:

*   `data/osm-prebuilt`: Contains OpenStreetMap (OSM) data that has already been retrieved and processed. This serves as a foundational network without incorporating specific national expansion plans.
*   `data/osm-plus-prebuilt`: Includes the processed OSM data along with the integrated National Power Development Plans. This option provides a more comprehensive base network reflecting planned expansions.

You can select which prebuilt network to use by configuring the `electricity` section in your `config.yaml`:

```yaml
electricity:
  base_network: osm-plus-prebuilt # Choose either 'osm-prebuilt' or 'osm-plus-prebuilt'
  base_network_version: 0.1.1 # Specify the version of the prebuilt network to use
```

The selection and integration of these base networks are managed through the Snakemake workflow, specifically within the `base_network` rule. This rule utilizes the `OSMDIR` variable, which is set based on the `base_network` and `base_network_version` configurations, to correctly access the prebuilt network files.

## ASEAN Power Grid (APG) Interconnections

The ASEAN Power Grid represents a crucial infrastructure for realizing an ASEAN-wide integrated energy system. PyPSA-ASEAN actively models the implementation of these interconnections, which are summarized within the `data/transmission_projects` directory.

The system categorizes transmission projects as follows:

*   **New Lines**: Refers to newly planned AC interconnections.
*   **New Links**: Refers to newly planned DC interconnections.
*   **Upgraded Lines/Links**: Indicates existing interconnections where capacity is planned to be increased.

PyPSA-ASEAN adopts and reframes a script originally developed for PyPSA-Eur (used for ENTSOE Ten Year Network Development Plan and Germany's National Development Plan) to apply to the context of the ASEAN Interconnection Masterplan (AIMS) and the Indonesia Supergrid (ID_Supergrid).

**Configuring Transmission Projects:**

Unlike the prebuilt networks, the integration of specific transmission projects (new lines, new links, and upgrades) is highly configurable via the `transmission_projects` section in your `config.yaml`. This section is frequently adjusted to model various scenarios and policy interventions.

Here's a breakdown of the key configuration options and their impact on the workflow, particularly on the `scripts/build_transmission_projects.py` and `scripts/add_transmission_projects.py` scripts:

```yaml
transmission_projects:
  enable: true
  include:
    ID_SuperGrid: true
    AIMS: true
  skip: []
  status:
  - in_planning
  - under_construction
  - existing
  - under_consideration
  new_link_capacity: keep #keep or zero
  distance_upper_bound: 0.3
  set_by_build_year: true
  delay_construction: 0
```

*   **`enable`**: (Boolean) When `true`, this global flag activates the processing and integration of all transmission projects defined within this section. If `false`, no transmission projects will be added or modified.

*   **`include`**: (Dictionary) This allows you to selectively include or exclude specific transmission project datasets. For instance, setting `ID_SuperGrid: true` and `AIMS: true` (as in the example) means the Indonesia Supergrid and ASEAN Interconnection Masterplan projects will be considered.

*   **`skip`**: (List of Strings) This option allows you to exclude specific individual files within the included project directories. For example, if you include `AIMS` but want to exclude a specific file named `problematic_lines.csv` within the `data/transmission_projects/AIMS/` directory, you would add `problematic_lines` to this list.

*   **`status`**: (List of Strings) This is a critical filter that determines which projects, based on their `project_status` attribute (e.g., `in_planning`, `under_construction`, `existing`), are considered for inclusion.

*   **`new_link_capacity`**: (String, `keep` or `zero`) This setting dictates how the capacity (`p_nom`) of newly added DC links (`new_links.csv`) is handled. If set to `keep`, the `p_nom` value specified in the `new_links.csv` file will be used. If set to `zero`, the capacity of new links will be initialized to 0, allowing the optimization model to determine their optimal capacity.

*   **`distance_upper_bound`**: (Float) This value (in degrees, approximately equivalent to 0.1 degrees ≈ 15 km) defines the maximum distance for matching new line/link endpoints to existing buses in the network. If a new line/link endpoint is further than this bound from any existing bus, the system will attempt to create a new bus for it, ensuring connectivity.

*   **`set_by_build_year`**: (Boolean) If `true`, the `build_year` attribute of projects will be considered. This is particularly relevant for myopic foresight models where infrastructure is added sequentially over time.

*   **`delay_construction`**: (Integer) This parameter allows you to simulate delays in the construction of new transmission projects. A value of `0` means no delay. If set to `X`, projects will be considered available `X` years after their planned `build_year`. This setting affects how projects are integrated into the network at different planning horizons.

These configuration options provide extensive control over how the ASEAN Power Grid infrastructure is modeled, enabling researchers to explore a wide range of future scenarios and policy impacts. By adjusting these parameters, you can fine-tune the representation of grid expansion and interconnection within your PyPSA-ASEAN simulations.
