# osm_data_config
# configuration paramaters necessary for osm_data_ level scripts
# OFFICIAL TWO DIGIT COUNTRY CODES WITH OSM/GEOFABRIK STRING NAMES
AFRICA_CC = {
    "DZ": "algeria",
    "AO": "angola",
    "BJ": "benin",
    "BW": "botswana",
    "BF": "burkina-faso",
    "BI": "burundi",
    "CM": "cameroon",
    # canary-islands    # Island
    # "CV": "cape-verde", # Island
    "CF": "central-african-republic",
    "TD": "chad",
    # "KM": "comores", # Island
    "CG": "congo-brazzaville",
    "CD": "congo-democratic-republic",
    "DJ": "djibouti",
    "EG": "egypt",
    "GQ": "equatorial-guinea",
    "ER": "eritrea",
    "ET": "ethiopia",
    "GA": "gabon",
    "GH": "ghana",
    # "GW": "guinea-bissau", # No Data
    "GN": "guinea",
    "CI": "ivory-coast",
    "KE": "kenya",
    "LS": "lesotho",
    "LR": "liberia",
    "LY": "libya",
    "MG": "madagascar",
    "MW": "malawi",
    "ML": "mali",
    "MR": "mauritania",
    "MU": "mauritius",
    "MA": "morocco",
    "MZ": "mozambique",
    "NA": "namibia",
    "NE": "niger",
    "NG": "nigeria",
    "RW": "rwanda",
    # saint-helena-ascension-and-tristan-da-cunha
    # "ST": "sao-tome-and-principe", #Island
    # "SN-GM": "senegal-and-gambia",  # Merged Country Code, See Map, See https://github.com/pypsa-meets-africa/pypsa-africa/issues/40#issuecomment-885424094
    # "SC": "seychelles", #Island
    "SL": "sierra-leone",
    # "SO": "somalia", # No Data
    # south-africa-and-lesotho
    "ZA": "south-africa",
    "SS": "south-sudan",
    "SD": "sudan",
    "SZ": "swaziland",
    "TZ": "tanzania",
    "TG": "togo",
    "TN": "tunisia",
    "UG": "uganda",
    "ZM": "zambia",
    "ZW": "zimbabwe",
}

# The feature category represents the final representation of the feature
# For node features: ways are converted to nodes
# For way features: only ways are used

feature_category = {
    "substation": "node",
    "generator": "node",
    "line": "way",
    "tower": "node",
    "cable": "way",
}

COMP_CC = {
    "DZ": "algeria",
    "AO": "angola",
    "CM": "cameroon",
    "CF": "central-african-republic",
    "TD": "chad",
    "CG": "congo-brazzaville",
    "CD": "congo-democratic-republic",
    "EG": "egypt",
    "GQ": "equatorial-guinea",
    "ER": "eritrea",
    "ET": "ethiopia",
    "CI": "ivory-coast",
    "LS": "lesotho",
    "LR": "liberia",
    "LY": "libya",
    "MR": "mauritania",
    "MU": "mauritius",
    "MA": "morocco",
    "MZ": "mozambique",
    "RW": "rwanda",
    "TZ": "tanzania",
    "TG": "togo",
    "TN": "tunisia",
    "ZM": "zambia",
    "ZW": "zimbabwe",
}


# ===============================
# OSM FEATURE COLUMNS
# ===============================
# These configurations are used to specify which OSM tags are kept as columns in DataFrame.
# Follows the OSM Wiki: https://wiki.openstreetmap.org/wiki/Power


# "Length" is added for way features
# "Area" is added for node features

# ========================
# BASIC INFO TAGS
# ========================
# A list of tags that are relevant for most OSM keys

columns_basic = [
    "id",
    "lonlat",
    "tags.power",
    "Type",
    "Country",
    # "refs"
]

# ========================
# SUBSTATION TAGS
# ========================

# Default tags to keep as columns with substation
# Based on: https://wiki.openstreetmap.org/wiki/Key:substation
columns_substation = [
    "Area",
    "tags.substation",
    "tags.voltage",

    # Other tags which are not kept by default
    # =====================================
    # "TODO:ADD Tags not kept here",
]

# ========================
# GENERATOR TAGS
# ========================

# Default tags to keep as columns with generator
# Based on: https://wiki.openstreetmap.org/wiki/Key:generator

columns_generator = [
    "Area",
    "tags.generator:type",
    "tags.generator:method",
    "tags.generator:source",
    "tags.generator:output:electricity",
    # Other tags which are not kept by default
    # =====================================
    # "TODO:ADD Tags not kept here",
]

# ========================
# LINE TAGS
# ========================

# Default tags to keep as columns with line
# Based on: https://wiki.openstreetmap.org/wiki/Key:line

columns_line = [
    "Length",
    "tags.cables",
    "tags.voltage",
    "tags.circuits",
    "tags.frequency",

    # Other tags which are not kept by default
    # =====================================
    # "TODO:ADD Tags not kept here",
]

# ========================
# CABLE TAGS
# ========================

# Default tags to keep as columns with substation
# Based on: https://wiki.openstreetmap.org/wiki/Key:cable

columns_cable = [
    "Length",
    "tags.cables",
    "tags.voltage",
    "tags.circuits",
    "tags.frequency",
    "tags.location",

    # Other tags which are not kept by default
    # =====================================
    # "TODO:ADD Tags not kept here",
]

# ========================
# TOWER TAGS
# ========================

# Default tags to keep as columns with tower
# Based on: https://wiki.openstreetmap.org/wiki/Key:tower

columns_tower = [
    "Area",
    "tags.tower",
    "tags.material",
    "tags.structure",
    "tags.operator",
    "tags.line_attachment",
    "tags.line_management",
    "tags.ref",
    "tags.height",

    # Other tags which are not kept by default
    # =====================================
    # "TODO:ADD Tags not kept here",
]

# FINAL DICTIONARY

feature_columns = {
    "substation": columns_basic + columns_substation,
    "generator": columns_basic + columns_generator,
    "line": columns_basic + columns_line,
    "cable": columns_basic + columns_cable,
    "tower": columns_basic + columns_tower,
}
