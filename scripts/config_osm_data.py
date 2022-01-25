# SPDX-FileCopyrightText: : 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# ===============================
# OSM FEATURE COLUMNS
# ===============================
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
    "tags.name",
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

# Python dictionary of ISO 3166-1-alpha-2 codes, as per publicly
# available data on official ISO site in July 2015.
#
# Available under MIT license
# Dimitris Karagkasidis, https://github.com/pageflt

continents = {
    "LA": "latin_america",
    "SA": "south_america",
    "CA": "central_america",
    "AS": "asia",
    "OC": "australia",
    "AF": "africa",
    "EU": "europe",
    # "AN": "antarctica"
}

world_iso = {
    "africa": {
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
        "GW": "guinea-bissau",  # No Data
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
        # "MU": "mauritius", Island
        "MA": "morocco",
        "MZ": "mozambique",
        "NA": "namibia",
        "NE": "niger",
        "NG": "nigeria",
        "RW": "rwanda",
        # saint-helena-ascension-and-tristan-da-cunha
        # "ST": "sao-tome-and-principe", #Island
        "SN": "senegal",
        "GM": "gambia",
        # "SC": "seychelles", #Island
        "SL": "sierra-leone",
        "SO": "somalia",  # No Data
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
        "EH": "western-sahara",
    },
    "asia": {
        "AF": "afghanistan",
        "AM": "armenia",
        "AZ": "azerbaijan",
        "BH": "bahrain",
        "BD": "bangladesh",
        "BT": "bhutan",
        # "IO": "british indian ocean territory",
        # "BN": "brunei darussalam", # merged with MY
        "KH": "cambodia",
        "CN": "china",
        # "CX": "christmas island",
        # "CC": "cocos (keeling) islands",
        "CY": "cyprus",
        "GE": "georgia",
        "HK": "hong kong",
        "IN": "india",
        "ID": "indonesia",
        "IR": "iran",
        "IQ": "iraq",
        "IL": "israel",
        "JP": "japan",
        "JO": "jordan",
        "KZ": "kazakhstan",
        "KP": "north-korea",
        "KR": "south-korea",
        "KW": "kuwait",
        "KG": "kyrgyzstan",
        "LA": "lao-people's-democratic-republic",
        "LB": "lebanon",
        "MO": "macao",
        "MY": "malaysia",
        "SG": "singapore",
        "BN": "brunei",
        # "MV": "maldives",  # Island
        "MN": "mongolia",
        "MM": "myanmar",
        "NP": "nepal",
        "OM": "oman",
        "PK": "pakistan",
        "PS": "palestine",
        "PH": "philippines",
        "QA": "qatar",
        "SA": "saudi-arabia",
        "SG": "singapore",  # merged with MY
        "LK": "sri-lanka",
        "SY": "syria",
        "TW": "taiwan",
        "TJ": "tajikistan",
        "TH": "thailand",
        "TM": "turkmenistan",
        "AE": "united-arab-emirates",
        "UZ": "uzbekistan",
        "VN": "vietnam",
        "YE": "yemen",
    },
    "australia": {
        # "AS": "american-oceania",  # Island
        "AU": "australia",
        # "CK": "cook islands",  # Island
        # "FJ": "fiji",  # Island
        # "PF": "french-polynesia",  # Island
        # "GU": "guam",  # Island
        # "KI": "kiribati",  # Island
        # "MH": "marshall islands",  # Island
        # "FM": "micronesia",  # Island
        # "NR": "nauru",  # Island
        "NC": "new-caledonia",
        "NZ": "new-zealand",
        # "NU": "niue",  # Island
        # "NF": "norfolk island",  # Island
        # "MP": "northern mariana islands",
        # "PW": "palau",  # Island
        "PG": "papua-new-guinea",
        # "WS": "samoa",  # Island
        # "SB": "solomon islands",
        # "TK": "tokelau",  # Island
        # "TO": "tonga",  # Island
        # "TV": "tuvalu",  # Island
        # "VU": "vanuatu",  # Island
        # "WF": "wallis-et-futuna",  # Island
    },
    "europe": {
        "AL": "albania",
        "AD": "andorra",
        "AT": "austria",
        "BY": "belarus",
        "BE": "belgium",
        "BA": "bosnia-herzegovina",
        "BG": "bulgaria",
        "HR": "croatia",
        "CZ": "czech-republic",
        "DK": "denmark",
        "EE": "estonia",
        # "FO": "faroe islands",
        "FI": "finland",
        "FR": "france",
        "DE": "germany",
        # "GI": "gibraltar", Island ?
        "GR": "greece",
        # "GG": "guernsey", Island
        "HU": "hungary",
        "IS": "iceland",
        "IE": "ireland-and-northern-ireland",
        # "IM": "isle of man",
        "IT": "italy",
        # "JE": "jersey",
        "LV": "latvia",
        "LI": "liechtenstein",
        "LT": "lithuania",
        "LU": "luxembourg",
        "MK": "macedonia",
        "MT": "malta",
        "MD": "moldova",
        "MC": "monaco",
        "ME": "montenegro",
        "NL": "netherlands",
        "NO": "norway",
        "PL": "poland",
        "PT": "portugal",
        "RO": "romania",
        "RU": "russia",
        # "SM": "san-marino",
        "RS": "serbia",
        "SK": "slovakia",
        "SI": "slovenia",
        "ES": "spain",
        # "SJ": "svalbard-and-jan-mayen",
        "SE": "sweden",
        "CH": "switzerland",
        "UA": "ukraine",
        "GB": "great-britain",
        "TR": "turkey",
    },
    "north_america": {
        "CA": "canada",
        "GL": "greenland",
        "MX": "mexico",
        "US": "united states",
    },
    "latin_america": {
        "AR": "argentina",
        "BO": "bolivia",
        "BR": "brazil",
        "CL": "chile",
        "CO": "colombia",
        "EC": "ecuador",
        "GF": "french-guyane",
        "GY": "guyane",
        "PE": "peru",
        "PY": "paraguay",
        "SR": "suriname",
        "UY": "uruguay",
        "VE": "venezuela",
    },
    "central_america": {
        "BZ": "belize",
        "CR": "costa-rica",
        "HN": "honduras",
        "GT": "guatemala",
        "NI": "nicaragua",
        "PA": "panama",
        "SV": "el-salvador",
    },
}

continent_regions = {
    # Based on: https://waml.org/waml-information-bulletin/46-3/index-to-lc-g-schedule/1-world/
    # Eurpean regions
    "SCR": ["DK", "NO", "SE", "FI", "IS"],  # SCANDANAVIAN REGION
    # EASTREN EUROPIAN REGION
    "EER": ["BY", "PL", "BU", "CZ", "RU", "SK", "UA", "LT", "LV", "EE", "SI"],
    # CENTRAL EUROPIAN REGION
    "CER": ["AT", "CH", "CZ", "DE", "HU", "PL", "SK"],
    # BALKAN PENISULAN REGION
    "BPR": [
        "AL",
        "AN",
        "BA",
        "BG",
        "GR",
        "HR",
        "MD",
        "MT",
        "RO",
        "SL",
        "RS",
        "ME",
        "MK",
    ],
    # WESTREN EUROPE
    "WER": ["FR", "BE", "GB", "IE", "LU", "MC", "NL", "AD", "LI"],
    "SER": ["ES", "IT", "PT"],  # SOUTHERN EUROPAIN REGION
    # African regions
    "NAR": ["EG", "DZ", "LY", "MA", "SD", "SS"],  # NORTHERN AFRICAN REGION
    # WESTREN AFRICAN REGION
    "WAR": [
        "MR",
        "ML",
        "NE",
        "NG",
        "BJ",
        "BF",
        "TG",
        "GH",
        "CI",
        "LR",
        "SL",
        "GN",
        "GM",
        "SL",
    ],
    # CENTRAL AFRICAN REGION
    "CAR": ["TD", "CF", "CM", "GQ", "GA", "CD", "CG", "AO"],
    # EASTREN AFRICAN REGION
    "EAR": ["ET", "UG", "KE", "RW", "BI", "TZ", "MZ", "DJ", "MG"],
    # SOUTHERN AFRICAN REGION
    "SAR": ["MW", "ZM", "ZW", "BW", "NA", "SZ", "LS", "ZA"],
    # Asian regions
    "WAS": [
        "TR",  # TURKEY
        "AM",  # ARMENIA
        "AZ",  # AZERBAIJAN
        "BH",  # BAHREIN
        "CY",  # CYPRUS (NORTH AND SOUTH)
        "GE",  # GEORGIA
        "IQ",  # IRAQ
        "IL",  # ISRAEL 
        "JO",  # JORDAN
        "KW",  # KUWAIT
        "LB",  # LEBANON
        "OM",  # OMAN
        "PS",  # PALESTINE
        "QA",  # QATAR
        "SA",  # SAUDI ARABIA
        "SY",  # SYRIA
        "AE",  # UNITED ARAB EMIRATES
        "YE",
        # XNC is missing for nothern cyprus
    ],  # YEMEN
    # FAR EASTREN AISIAN REGION
    "FEAR": ["JP", "KP", "KR", "CN", "TW", "CN", "MN", "HK", "MO"],
    # SOUTHEASTREN AISIAN REGION
    "SEAR": ["LA", "TH", "KH", "VN", "PH", "MY", "SG", "BN", "ID"],
    "SG":
    "MY",  # Singapore -> Malaysia
    "BN":
    "MY",  # Brunei -> Malaysia
    "CAR": ["KZ", "KG", "UZ", "TM", "TJ"],  # CENTRAL AISIAN REGION
    # SOUTHERN AISIAN REGION
    "SAR": ["MM", "BD", "BT", "NP", "IN", "LK", "PK", "AF"],
    # MIDDLE EASTREN ASIAN REGION
    "MEAR": ["TR", "SY", "IQ", "IR", "JO", "IL", "AE", "YE"],
    # American continent regions
    "NACR": ["CA", "GL", "MX", "US"],  # NORTHERN AMERCAN CONTINENT REGION
    # SOUTHERN LATIN AMERICAN REGION
    "LACR": ["AR", "BO", "BR", "CL", "CO", "EC", "PE", "SR", "UY", "VE"],
    "CACR": ["BZ", "GT", "SV", "HN", "NI", "CR"],  # CENTRAL AMERICAN REGION
    # Customized test set
    "TEST": ["NG", "NE", "SL", "MA"],
}

# Geofabrik and iso norm deviate for some countries and domains

# dictionary of correspondance between iso country codes and geofabrik codes containing those information
# This dictionary instructs the script download_osm_data about how to successfully download data
# from countries that are aggregated into osm.
# For example, Senegal (SN) and Gambia (GM) cannot be downloaded from OSM separately, but only jointly as SNGM
#   That's the reason why in this dictionary they can be found the following entries:
#       "SN": "SNGM"
#       "GM": "SNGM"
#   This instruct the workflow that when the country "SN" is requested, then it shall download the "SNGM" file
iso_to_geofk_dict = {
    "EH": "MA",  # Western Sahara -> Morocco
    "SN": "SNGM",  # Senegal -> Senegal-Gambia
    "GM": "SNGM",  # Gambia -> Senegal-Gambia
    "HK": "CN",  # Hong Kong  -> China
    "MO": "CN",  # Macao  -> China
    "SG": "MY-SG-BN",  # Singapore -> Malaysia-Singapore-Brunei
    "BN": "MY-SG-BN",  # Brunei -> Malaysia-Singapore-Brunei
    "MY": "MY-SG-BN",  # Malaysia -> Malaysia-Singapore-Brunei
    "SA": "GCC",  # Saudi Arabia -> Gulf Cooperation Council
    "KW": "GCC",  # Kuwait -> Gulf Cooperation Council
    "BH": "GCC",  # Bahrain -> Gulf Cooperation Council
    "QA": "GCC",  # Qatar -> Gulf Cooperation Council
    "AE": "GCC",  # United Arab Emirates -> Gulf Cooperation Council
    "OM": "GCC",  # Oman -> Gulf Cooperation Council
    "PS": "IL-PL",
    "IL": "IL-PL",
}

# Cyprus and Georgia -> European domain
# Russia -> a separate domain

# data for some islands seem to be merged with some other areas data
# "FO": "faroe islands"
# "NF": "norfolk island",
# "PF": "french-polynesia"
# "GU": "guam"

# "latin_america" -> "south-america"

world_geofk = {
    "africa": {
        "DZ": "algeria",
        "AO": "angola",
        "BJ": "benin",
        "BW": "botswana",
        "BF": "burkina-faso",
        "BI": "burundi",
        "CM": "cameroon",
        # canary-islands, # Island
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
        "GW": "guinea-bissau",  # No Data
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
        # "MU": "mauritius", # Island
        "MA": "morocco",
        "MZ": "mozambique",
        "NA": "namibia",
        "NE": "niger",
        "NG": "nigeria",
        "RW": "rwanda",
        # saint-helena-ascension-and-tristan-da-cunha # Islands
        # "ST": "sao-tome-and-principe", # Island
        "SNGM": "senegal-and-gambia",  # Geofk shortcurt
        # "SC": "seychelles", # Island
        "SL": "sierra-leone",
        "SO": "somalia",  # No Data
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
    },
    "asia": {
        "AF": "afghanistan",
        "AM": "armenia",
        "AZ": "azerbaijan",
        "BD": "bangladesh",
        "BT": "bhutan",
        "KH": "cambodia",
        "CN": "china",
        "GCC": "gcc-states",  # Geofk shortcurt for SA, KW, BH, QA, AE, OM
        "IN": "india",
        "ID": "indonesia",
        "IR": "iran",
        "IQ": "iraq",
        "IL-PL": "israel-and-palestine",
        "JP": "japan",
        "JO": "jordan",
        "KZ": "kazakhstan",
        "KP": "north-korea",
        "KR": "south-korea",
        "KG": "kyrgyzstan",
        "LA": "laos",
        "LB": "lebanon",
        "MY-SG-BN": "malaysia-singapore-brunei",  # Geofk shortcurt
        "MN": "mongolia",
        "MM": "myanmar",
        "NP": "nepal",
        "PK": "pakistan",
        "PH": "philippines",
        "LK": "sri-lanka",
        "SY": "syria",
        "TW": "taiwan",
        "TJ": "tajikistan",
        "TH": "thailand",
        "TM": "turkmenistan",
        "UZ": "uzbekistan",
        "VN": "vietnam",
        "YE": "yemen",
    },
    "australia-oceania": {
        "AU": "australia",
        # "CK": "cook islands",  # Island
        # "FJ": "fiji",  # Islands
        # "PF": "french-polynesia",  # Islands
        # "GU": "guam",  # Island
        # "KI": "kiribati",  # Islands
        # "MH": "marshall islands",  # Islands
        # "FM": "micronesia",  # Islands
        # "NR": "nauru",  # Islands
        "NC": "new-caledonia",
        "NZ": "new-zealand",
        # "NU": "niue",  # Island
        # "NF": "norfolk island",  # Island
        # "MP": "northern mariana islands",  # Islands
        # "PW": "palau",  # Islands
        "PG": "papua-new-guinea",
        # "WS": "samoa",  # Islands
        # "SB": "solomon islands",  # Islands
        # "TK": "tokelau",  # Islands
        # "TO": "tonga",  # Islands
        # "TV": "tuvalu",  # Islands
        # "VU": "vanuatu",  # Islands
        # "WF": "wallis-et-futuna",  # Islands
    },
    "europe": {
        "AL": "albania",
        "AD": "andorra",
        "AT": "austria",
        "BY": "belarus",
        "BE": "belgium",
        "BA": "bosnia-herzegovina",
        "BG": "bulgaria",
        "HR": "croatia",
        "CZ": "czech-republic",
        "CY": "cyprus",
        "DK": "denmark",
        "EE": "estonia",
        # "FO": "faroe islands", # Islands
        "FI": "finland",
        "FR": "france",
        "GE": "georgia",
        "DE": "germany",
        # "GI": "gibraltar", # Peninsula; Isolated PS?
        "GR": "greece",
        # "GG": "guernsey", # Island
        "HU": "hungary",
        "IS": "iceland",
        "IE": "ireland-and-northern-ireland",
        # "IM": "isle of man", # Island
        "IT": "italy",
        # "JE": "jersey", # Island
        "LV": "latvia",
        "LI": "liechtenstein",
        "LT": "lithuania",
        "LU": "luxembourg",
        "MK": "macedonia",
        "MT": "malta",
        "MD": "moldova",
        "MC": "monaco",
        "ME": "montenegro",
        "NL": "netherlands",
        "NO": "norway",
        "PL": "poland",
        "PT": "portugal",
        "RO": "romania",
        "SM": "san-marino",
        "RS": "serbia",
        "SK": "slovakia",
        "SI": "slovenia",
        "ES": "spain",
        # "SJ": "svalbard-and-jan-mayen", # Islands
        "SE": "sweden",
        "CH": "switzerland",
        "UA": "ukraine",
        "GB": "great-britain",
        "TR": "turkey",
    },
    "russia": {
        "CEFD": "central-fed-district",
        "FEFD": "far-eastern-fed-district",
        "NCDF": "north-caucasus-fed-district",
        "NWDF": "northwestern-fed-district",
        "SBFD": "siberian-fed-district",
        "SOFD": "south-fed-district",
        "URDF": "ural-fed-district",
        "VOFD": "volga-fed-district",
        "RU": "russia",
    },
    "north-america": {
        "CA": "canada",
        "GL": "greenland",
        "MX": "mexico",
        "US": "us",
    },
    "south-america": {
        "AR": "argentina",
        "BO": "bolivia",
        "BR": "brazil",
        "CL": "chile",
        "CO": "colombia",
        "EC": "ecuador",
        "PE": "peru",
        "SR": "suriname",
        "PY": "paraguay",
        "UY": "uruguay",
        "VE": "venezuela",
    },
    "central-america": {
        "BZ": "belize",
        "GT": "guatemala",
        "SV": "el-salvador",
        "HN": "honduras",
        "NI": "nicaragua",
        "CR": "costa-rica",
    },
}

world_countries = {
    country_2D: country_name
    for d in world_geofk.values() for (country_2D, country_name) in d.items()
}
