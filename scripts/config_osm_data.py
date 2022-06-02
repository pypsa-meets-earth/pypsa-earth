# -*- coding: utf-8 -*-
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
    "LA": "NorthAmerica",
    "SA": "SouthAmerica",
    "AS": "Asia",
    "OC": "Oceania",
    "AF": "Africa",
    "EU": "Europe",
    # "AN": "antarctica"
}

world_iso = {
    "Africa": {
        "DZ": "algeria",
        "AO": "angola",
        "BJ": "benin",
        "BW": "botswana",
        # "IO": "british-indian-ocean-territory", # Island
        "BF": "burkina-faso",
        "BI": "burundi",
        "CM": "cameroon",
        # "IC": "canary-islands"    # Island
        # "CV": "cape-verde", # Island
        "CF": "central-african-republic",
        "TD": "chad",
        # "KM": "comoros", # Island
        "CG": "congo-brazzaville",
        "CD": "congo-democratic-republic",
        "DJ": "djibouti",
        "EG": "egypt",
        "GQ": "equatorial-guinea",
        "ER": "eritrea",
        "ET": "ethiopia",
        # "TF": "french-southern-territories",  # Island
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
        # "YT": "mayotte",  # Island
        "MA": "morocco",
        "MZ": "mozambique",
        "NA": "namibia",
        "NE": "niger",
        "NG": "nigeria",
        # "RE": "reunion",  # Island
        "RW": "rwanda",
        # saint-helena-ascension-and-tristan-da-cunha # Islands
        # "ST": "sao-tome-and-principe", # Island
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
    "Asia": {
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
        # "CX": "christmas island", # Island
        # "CC": "cocos (keeling) islands", # Island
        "CY": "cyprus",
        # "EG": "egypt",  # leads to bug -> missing ssp file when executing ["Africa"]
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
        "RU": "russian-federation",
        "SA": "saudi-arabia",
        "SG": "singapore",  # merged with MY
        # "XS": "spratly-islands", #Island
        "LK": "sri-lanka",
        "SY": "syria",
        "TW": "taiwan",
        "TJ": "tajikistan",
        "TH": "thailand",
        "TL": "timor-leste",
        "TR": "turkey",
        "TM": "turkmenistan",
        "AE": "united-arab-emirates",
        # "XD": "united-nations-neutral-zone"
        "UZ": "uzbekistan",
        "VN": "vietnam",
        "YE": "yemen",
    },
    "Oceania": {
        # "AS": "american-samoa",  # Island
        "AU": "australia",
        # "CK": "cook islands",  # Island
        # "FJ": "fiji",  # Island
        # "PF": "french-polynesia",  # Island
        # "GU": "guam",  # Island
        # "KI": "kiribati",  # Island
        # "MH": "marshall-islands",  # Island
        # "FM": "micronesia",  # Island
        # "NR": "nauru",  # Island
        "NC": "new-caledonia",
        "NZ": "new-zealand",
        # "NU": "niue",  # Island
        # "NF": "norfolk-island",  # Island
        # "MP": "northern-mariana-islands",
        # "PW": "palau",  # Island
        # "PN": "pitcairn-islands", # Islands
        # "PW": "palau",  # Island
        # "WS": "samoa",  # Island
        # "SB": "solomon-islands",
        # "TK": "tokelau",  # Island
        # "TO": "tonga",  # Island
        # "TV": "tuvalu",  # Island
        # "UM": "united-states-minor-outlying-islands", #Islands
        # "VU": "vanuatu",  # Island
        # "WF": "wallis-and-futuna",  # Island
    },
    "Europe": {
        # "AX":"aland-islands", # Island
        "AL": "albania",
        "AD": "andorra",
        "AM": "armenia",
        "AT": "austria",
        "AZ": "Azerbaijan",
        "BY": "belarus",
        "BE": "belgium",
        "BA": "bosnia-herzegovina",
        "BG": "bulgaria",
        "HR": "croatia",
        "CY": "cyprus",
        "CZ": "czech-republic",
        "DK": "denmark",
        "EE": "estonia",
        # "FO": "faroe islands", # Islands
        "FI": "finland",
        "FR": "france",
        "GE": "georgia",
        "DE": "germany",
        # "GI": "gibraltar", # Island ?
        "GR": "greece",
        # "GG": "guernsey", # Island
        "HU": "hungary",
        "IS": "iceland",
        "IE": "ireland-and-northern-ireland",
        # "IM": "isle of man", # Island
        "IT": "italy",
        # "JE": "jersey", # Island
        "KZ": "kazakhstan",
        "XK": "kosovo",
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
        "SM": "san-marino",
        "RS": "serbia",
        "SK": "slovakia",
        "SI": "slovenia",
        "ES": "spain",
        # "SJ": "svalbard-and-jan-mayen",  # Islands
        "SE": "sweden",
        "CH": "switzerland",
        "UA": "ukraine",
        "GB": "great-britain",
        "TR": "turkey",
        "VA": "vatican",
    },
    "NorthAmerica": {
        "AI": "anguilla",
        # "AG": "antigua-and-barbuda", # Islands
        # "AW": "aruba", # Islands
        # "BS": "bahamas", # Islands
        # "BB": "barbados", # Islands
        # "BM": "bermuda", # Islands
        # "BQ": "bonaire", # Islands
        # "VG": "british-virgin-islands", # Islands
        "CA": "canada",
        # "KY": "cayman-islands", # Islands
        # "CU": "cuba", # Islands
        # "CW": "curacao", # Islands
        # "DM": "dominica", # Islands
        "DO": "dominican-republic",
        "GL": "greenland",
        # "GD": "grenada", # Islands
        # "GP": "guadeloupe", # Islands
        "HT": "haiti",
        # "JM": "jamaica", # Islands
        # "MQ": "martinique", # Islands
        "MX": "mexico",
        # "MS": "montserrat", # Islands
        "US": "united-states-of-america",
        # "PR": "puerto-rico", # Islands
        # "BL": "saint-barthelemy", # Islands
        # "KN": "saint-kitts-and-nevis", # Islands
        # "LC": "saint-lucia", # Islands
        # "MF": "saint-martin", # Islands
        # "PM": "saint-pierre-and-miquelon", # Islands
        # "VC": "saint-vincent-and-the-grenadines", # Islands
        # "SX": "saint-marteen", # Islands
        # "TT": "trinidad-and-tobago", # Islands
        # "TC": "turks-and-caicos", # Islands
        # "UM": "united-states-minor-outlying-islands", #Islands
        # "VI": "united-states-virgin-islands", #Islands
        "BZ": "belize",
        "CR": "costa-rica",
        "HN": "honduras",
        "GT": "guatemala",
        "NI": "nicaragua",
        "PA": "panama",
        "SV": "el-salvador",
    },
    "SouthAmerica": {
        "AR": "argentina",
        "BO": "bolivia",
        "BR": "brazil",
        "CL": "chile",
        "CO": "colombia",
        "EC": "ecuador",
        # "FK": "falkland-islands", #Islands
        "GF": "french-guiana",
        # "GY": "guyana", # No Data
        "PE": "peru",
        "PY": "paraguay",
        "SR": "suriname",
        "UY": "uruguay",
        "VE": "venezuela",
    },
    # "Antarctica": {
    #   "AQ": "antarctica",
    #   "BV": "bouvet-island",
    #   "HM": "heard-island-and-mcdonald-island",
    #   "GS": "south-georgia-and-the-south-sandwich-islands",
    # },
}

# Based on: https://waml.org/waml-information-bulletin/46-3/index-to-lc-g-schedule/1-world/
# Australasia region includes New Caledonia and Papua New Guinea
continent_regions = {
    # European regions
    "SCR": ["DK", "NO", "SE", "FI", "IS"],  # SCANDANAVIAN REGION
    # EASTREN EUROPIAN REGION
    "EER": ["BY", "PL", "CZ", "RU", "SK", "UA", "LT", "LV", "EE", "FI", "MD"],
    # CENTRAL EUROPIAN REGION
    "CER": ["AT", "CH", "CZ", "DE", "HU", "PL", "SK", "LI"],
    # BALKAN PENISULAN REGION
    "BPR": ["AL", "BA", "BG", "GR", "HR", "ME", "RO", "SI", "RS", "ME", "MK"],
    # WESTREN EUROPE
    "WER": ["FR", "BE", "GB", "IE", "LU", "MC", "NL", "AD"],
    # SOUTHERN EUROPAIN REGION
    "SER": ["ES", "AD", "IT", "PT", "SM", "MT"],
    # African regions
    # NORTHERN AFRICAN REGION
    "NAR": ["EG", "LY", "TN", "DZ", "MA", "EH", "SD", "SS"],
    # WESTREN AFRICAN REGION
    # Guinea-Bissau ["GW"] belongs to the region but power data are NA in OSM)
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
        "SN",
        "GM",
    ],
    # CENTRAL AFRICAN REGION
    "CAR": ["TD", "CF", "CM", "GQ", "GA", "CD", "CG", "AO"],
    # EASTREN AFRICAN REGION
    # Somalia ["SO"] belongs to the region but power data are NA in OSM)
    "EAR": ["ER", "ET", "UG", "KE", "RW", "BI", "TZ", "MZ", "DJ", "MG"],
    # SOUTHERN AFRICAN REGION
    "SAR": ["MW", "ZM", "ZW", "BW", "NA", "SZ", "LS", "ZA"],
    # Asian regions
    "KVR": ["AZ", "GE", "AM"],
    "WAS": [
        "TR",
        "AM",
        "AZ",
        "BH",
        "CY",
        "GE",
        "IQ",
        "IL",
        "JO",
        "KW",
        "LB",
        "OM",
        "PS",
        "QA",
        "SA",
        "SY",
        "AE",
        "YE",
    ],
    # FAR EASTREN AISIAN REGION
    "FEAR": ["JP", "KP", "KR", "CN", "TW", "MN", "HK", "MO"],
    # SOUTHEASTREN AISIAN REGION
    "SEAR": ["LA", "TH", "KH", "VN", "PH", "MY", "SG", "BN", "ID"],
    # CENTRAL AISIAN REGION
    "CASR": ["KZ", "KG", "UZ", "TM", "TJ"],
    # SOUTHERN AISIAN REGION
    "SASR": ["MM", "BD", "BT", "NP", "IN", "LK", "PK", "AF"],
    # MIDDLE EASTREN ASIAN REGION
    "MEAR": [
        "TR",
        "SY",
        "LB",
        "CY",
        "IQ",
        "IR",
        "JO",
        "IL",
        "PS",
        "AE",
        "YE",
        "KW",
        "BH",
        "QA",
        "SA",
        "OM",
    ],
    # American continent regions
    "NACR": ["CA", "GL", "MX", "US"],  # NORTHERN AMERCAN CONTINENT REGION
    # SOUTHERN LATIN AMERICAN REGION
    "LACR": ["AR", "BO", "BR", "CL", "CO", "EC", "GF", "PE", "PY", "SR", "UY", "VE"],
    # CENTRAL AMERICAN REGION
    "CACR": ["BZ", "GT", "SV", "HN", "NI", "CR", "PA"],
    # Australasia
    "AUO": ["AU", "NC", "NZ", "PG"],
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
    "PS": "IL-PL",  # Israel and Palestine are merged in OSM
    "IL": "IL-PL",  # Israel and Palestine are merged in OSM
    "SM": "IT",  # San-Marino is merged to Italy
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
        "PA": "panama",
    },
    "europe/france": {
        "GF": "guyane",
    },
}

world_countries = {
    country_2D: country_name
    for d in world_geofk.values()
    for (country_2D, country_name) in d.items()
}
