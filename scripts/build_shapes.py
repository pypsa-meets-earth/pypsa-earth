# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021 PyPSA-Africa authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
import logging
import multiprocessing as mp
import os
import shutil
import zipfile
from itertools import takewhile
from operator import attrgetter

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import requests
import rioxarray as rx
import xarray as xr
from _helpers import (
    configure_logging,
    sets_path_to_root,
    three_2_two_digits_country,
    two_2_three_digits_country,
    two_digits_2_name_country,
)
from rasterio.mask import mask
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from shapely.validation import make_valid
from tqdm import tqdm

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)

sets_path_to_root("pypsa-earth")


def download_GADM(country_code, update=False, out_logging=False):
    """
    Download gpkg file from GADM for a given country code

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

    GADM_filename = f"gadm41_{two_2_three_digits_country(country_code)}"
    GADM_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/{GADM_filename}.gpkg"

    GADM_inputfile_gpkg = os.path.join(
        os.getcwd(),
        "data",
        "gadm",
        GADM_filename,
        GADM_filename + ".gpkg",
    )  # Input filepath gpkg

    if not os.path.exists(GADM_inputfile_gpkg) or update is True:
        if out_logging:
            _logger.warning(
                f"Stage 4/4: {GADM_filename} of country {two_digits_2_name_country(country_code)} does not exist, downloading to {GADM_inputfile_gpkg}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_gpkg), exist_ok=True)

        with requests.get(GADM_url, stream=True) as r:
            with open(GADM_inputfile_gpkg, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return GADM_inputfile_gpkg, GADM_filename


def get_GADM_layer(country_list, layer_id, geo_crs, update=False, outlogging=False):
    """
    Function to retrive a specific layer id of a geopackage for a selection of countries

    Parameters
    ----------
    country_list : str
        List of the countries
    layer_id : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested

    """
    # initialization of the geoDataFrame
    geodf_list = []

    for country_code in country_list:
        # download file gpkg
        file_gpkg, name_file = download_GADM(country_code, update, outlogging)

        # # get layers of a geopackage
        list_layers = fiona.listlayers(file_gpkg)

        # # get layer name
        if (layer_id < 0) | (layer_id >= len(list_layers)):
            # when layer id is negative or larger than the number of layers, select the last layer
            layer_id = len(list_layers) - 1

        # read gpkg file
        geodf_temp = gpd.read_file(file_gpkg, layer="ADM_ADM_" + str(layer_id)).to_crs(
            geo_crs
        )

        # convert country name representation of the main country (GID_0 column)
        geodf_temp["GID_0"] = [
            three_2_two_digits_country(twoD_c) for twoD_c in geodf_temp["GID_0"]
        ]

        # TODO Check results for India (currently fixed manually)
        # { "type": "Feature", "properties": { "GADM_ID": "Z04.13_1", "country": null, "pop": 0.0, "gdp": 25389020.0 }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 78.651345210000045, 32.092282259000172 ], [ 78.7314071670001, 32.00764847 ], [ 78.771026610000035, 31.999811171000033 ], [ 78.725717640000141, 31.886429265000061 ], [ 78.729984178, 31.911489485000175 ], [ 78.69931415800005, 31.965185464000058 ], [ 78.631230238000171, 32.025268473000153 ], [ 78.601231257000052, 32.077213504000099 ], [ 78.522762238000155, 32.120071503000133 ], [ 78.450207298000123, 32.189766424000027 ], [ 78.439732199, 32.20926645500009 ], [ 78.449314317000074, 32.234630435000156 ], [ 78.464118958000086, 32.245960239 ], [ 78.539421082000047, 32.24998092800007 ], [ 78.592910768000195, 32.218669891000104 ], [ 78.593696591000139, 32.182418821000056 ], [ 78.574676512000167, 32.157730101000141 ], [ 78.597335817000157, 32.122032166999986 ], [ 78.639610289000132, 32.112121581000167 ], [ 78.651345210000045, 32.092282259000172 ] ] ] } },
        # { "type": "Feature", "properties": { "GADM_ID": "Z09.13_1", "country": null, "pop": 0.0, "gdp": 4672460.0 }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 78.777839659, 31.401430130999984 ], [ 78.827217100000155, 31.38128089900016 ], [ 78.84028625100018, 31.354440690999979 ], [ 78.829010009000115, 31.330200197000181 ], [ 78.858245849000014, 31.312980649000053 ], [ 78.90685200300004, 31.246172194000053 ], [ 78.848414582000089, 31.308830324000155 ], [ 78.819583443000113, 31.294107315000076 ], [ 78.782832481000128, 31.298885326000118 ], [ 78.759280562000185, 31.316608305000045 ], [ 78.746196542000121, 31.361664376000192 ], [ 78.782336582000198, 31.443024376000096 ], [ 78.740642462000096, 31.488691365000193 ], [ 78.795478820000142, 31.452762600000085 ], [ 78.792480469999987, 31.417858121000165 ], [ 78.777839659, 31.401430130999984 ] ] ] } },
        # { "type": "Feature", "properties": { "GADM_ID": "Z01.14_1", "country": null, "pop": 0.0, "gdp": 61683032064.0 }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 75.071609497000111, 32.482963563 ], [ 75.030166626000096, 32.491760255000088 ], [ 74.979347229000041, 32.447769166000114 ], [ 74.933631897000112, 32.467571259000124 ], [ 74.899658203000172, 32.466129304000106 ], [ 74.835327147999976, 32.496212006000121 ], [ 74.758026123000093, 32.479938507000043 ], [ 74.698417663000157, 32.484939575999988 ], [ 74.681381226000155, 32.503810883000085 ], [ 74.680511474000127, 32.538188934000175 ], [ 74.65239715600012, 32.560161591000167 ], [ 74.642211914000086, 32.61104965300018 ], [ 74.676208497000175, 32.658000945000026 ], [ 74.650138855000193, 32.718410491000043 ], [ 74.694580078000172, 32.809520723000048 ], [ 74.697021485000107, 32.839714050000055 ], [ 74.64199829100005, 32.82289886500007 ], [ 74.641540527000018, 32.791889191000109 ], [ 74.614234925000119, 32.757732391000161 ], [ 74.531532287000118, 32.740730287000076 ], [ 74.466056824000134, 32.781360626000037 ], [ 74.430061340000066, 32.844654084000069 ], [ 74.428489684000056, 32.883602142000143 ], [ 74.413894653000114, 32.903949738000108 ], [ 74.318031311000141, 32.919700623000097 ], [ 74.320716859000129, 32.941429139000093 ], [ 74.353179932000046, 32.98313140800002 ], [ 74.349830628000177, 33.016250610000156 ], [ 74.254966735000039, 33.048774719000107 ], [ 74.20524597300016, 33.043300630000033 ], [ 74.176246644000173, 33.107151032000161 ], [ 74.111701965000066, 33.16444015400009 ], [ 74.014274597000167, 33.1982498160001 ], [ 74.020301819000053, 33.258781432999967 ], [ 74.102981568000189, 33.270542146000082 ], [ 74.170745849000014, 33.346649170000035 ], [ 74.18508148300009, 33.383510590000128 ], [ 74.187072754000098, 33.455520629000091 ], [ 74.178138733000139, 33.48168182500018 ], [ 74.095710755000141, 33.57039642400008 ], [ 74.028442383000197, 33.57468032800017 ], [ 73.956253053000069, 33.720291138000107 ], [ 74.065727235000168, 33.820098878000124 ], [ 74.143341065000072, 33.830139161000147 ], [ 74.219673156000056, 33.867176056000119 ], [ 74.261123656000109, 33.92479705900007 ], [ 74.262741089000087, 33.974288941000168 ], [ 74.248687744000165, 34.013820648000149 ], [ 74.214675902000124, 34.038555145000032 ], [ 74.161422729000094, 34.038871765000124 ], [ 74.123695374000022, 34.055179596000016 ], [ 74.00607299700016, 34.022441864000086 ], [ 73.898170471000014, 34.03282165600001 ], [ 73.883216858000083, 34.049949646000016 ], [ 73.904083253000067, 34.122249604000103 ], [ 73.976577758000133, 34.211803436000082 ], [ 73.976226806000057, 34.264808654000149 ], [ 73.89940643400007, 34.359813691000056 ], [ 73.841751100000067, 34.338321686000029 ], [ 73.807395935000102, 34.356719971000075 ], [ 73.773391724000078, 34.349750519000168 ], [ 73.764320374000192, 34.371070861000135 ], [ 73.77698516800001, 34.416301728 ], [ 73.898742675000108, 34.49576950200003 ], [ 73.886688233000143, 34.514999390000071 ], [ 73.895339966000051, 34.546669007 ], [ 73.948066712000127, 34.573089600000174 ], [ 73.927330017000031, 34.645889283000088 ], [ 73.960220336000134, 34.696681976000093 ], [ 74.012893677000136, 34.709854125999982 ], [ 74.049133300000051, 34.684959412000069 ], [ 74.141792297000165, 34.691200257000162 ], [ 74.20295715300017, 34.73657989600008 ], [ 74.279747009000118, 34.768268586000147 ], [ 74.307296753000173, 34.800304412000116 ], [ 74.341262817000143, 34.78969192500017 ], [ 74.375541688000169, 34.803443909000066 ], [ 74.580177306999985, 34.770099639000193 ], [ 74.671302795000031, 34.70069885200013 ], [ 74.703002930000196, 34.693080903000123 ], [ 74.83560180600017, 34.677551270000151 ], [ 74.910507203000066, 34.683208465000064 ], [ 75.018745423000041, 34.641498566000166 ], [ 75.142295837000063, 34.662410736000083 ], [ 75.249710083000082, 34.647312165000073 ], [ 75.266059875, 34.63954925600018 ], [ 75.261856080000086, 34.610881806000179 ], [ 75.377822876000039, 34.552421569000103 ], [ 75.615562439000144, 34.538539886000024 ], [ 75.749206544000174, 34.515830994000112 ], [ 75.843276977000073, 34.574630737000064 ], [ 75.918647766000106, 34.604839325000171 ], [ 75.946250915, 34.631179810000106 ], [ 75.991851807000103, 34.629970550000053 ], [ 76.036460877000138, 34.670269013000052 ], [ 76.074440002000074, 34.676311494000061 ], [ 76.121932983000079, 34.646652222000057 ], [ 76.157913208000195, 34.64340973000003 ], [ 76.260047912000061, 34.684371949000024 ], [ 76.316139220000139, 34.731632234000188 ], [ 76.385932923000041, 34.735950471000137 ], [ 76.474342346000128, 34.79437255900018 ], [ 76.561813355000197, 34.758338928000114 ], [ 76.681816101000152, 34.759590149000189 ], [ 76.743240356000172, 34.840099334000115 ], [ 76.738861084000177, 34.901927949000026 ], [ 76.749366760999976, 34.926601409999989 ], [ 76.811050414000192, 34.935455323000099 ], [ 76.870040893000066, 34.972431183000083 ], [ 76.970359802000189, 34.935581208000144 ], [ 76.997489930000029, 34.939151763000041 ], [ 77.010139466000169, 34.956287384000063 ], [ 77.018096924000076, 35.038463592000028 ], [ 77.068367004000095, 35.024120330000187 ], [ 77.15682220400015, 35.054698945000155 ], [ 77.650566101000152, 35.367080688000101 ], [ 77.823928834000071, 35.501331329000038 ], [ 77.846778873, 35.500782009000091 ], [ 77.879417419000163, 35.436248779000096 ], [ 77.9029693600001, 35.427589416000103 ], [ 77.961723328000062, 35.480060577000131 ], [ 78.07935333100005, 35.492198942000186 ], [ 78.10398864900003, 35.484458922000158 ], [ 78.107696534000127, 35.46456909200009 ], [ 78.021202090000145, 35.416881560000036 ], [ 77.990287780000017, 35.379032140000049 ], [ 78.00202200200016, 35.29462100000012 ], [ 78.043186190000085, 35.201829910000185 ], [ 78.124176030000115, 35.135421751000081 ], [ 78.117095947000053, 35.111640931000125 ], [ 78.0828704810001, 35.082489011000177 ], [ 78.184646607000104, 34.881519318000073 ], [ 78.207321167000032, 34.766109467000092 ], [ 78.248031617000095, 34.695098877000135 ], [ 78.300338746000136, 34.661930085000051 ], [ 78.398361207000164, 34.648632050000117 ], [ 78.41037633000019, 34.61252513300019 ], [ 78.42511749300013, 34.60346222000004 ], [ 78.56933594000003, 34.611949921000189 ], [ 78.660057067000082, 34.529720307000105 ], [ 78.73007965100004, 34.496650700000146 ], [ 78.82495117100018, 34.414131164000025 ], [ 78.922546387000125, 34.375389098000085 ], [ 78.990050999000061, 34.317829000000131 ], [ 78.993400573000088, 34.28874206800009 ], [ 78.975006101000133, 34.242858888000114 ], [ 78.94211578200003, 34.204792022000106 ], [ 78.871147001000054, 34.157230001000073 ], [ 78.729233001000125, 34.095698999000035 ], [ 78.706382999000084, 34.072181964000151 ], [ 78.712120000000198, 34.05048 ], [ 78.76064300000013, 34.036098 ], [ 78.770361926000078, 34.019725797000149 ], [ 78.763626100000181, 34.001510621000079 ], [ 78.71579000000014, 33.964531001000125 ], [ 78.711563111000032, 33.922500611000089 ], [ 78.726218500000186, 33.87006694400003 ], [ 78.783988909000072, 33.77043175100016 ], [ 78.706787002000056, 33.669478990000187 ], [ 78.764022235000198, 33.561838556000055 ], [ 78.823257001000172, 33.483490000000131 ], [ 79.070793, 33.29092 ], [ 79.074966000000188, 33.228611 ], [ 79.18547094000013, 33.191223030000117 ], [ 79.19366508000013, 33.142970970000079 ], [ 79.148865059000173, 33.07480596 ], [ 79.15466303900007, 33.023777039000151 ], [ 79.248695041000076, 32.99347296000002 ], [ 79.280220060000147, 32.944331969000075 ], [ 79.278778079000119, 32.882888970000124 ], [ 79.249779000000103, 32.832137970000076 ], [ 79.28891802000004, 32.793166981000127 ], [ 79.291083959000105, 32.759055991000082 ], [ 79.264580039000066, 32.67419400100016 ], [ 79.293662999000162, 32.667416010000011 ], [ 79.264052999000171, 32.589890008999987 ], [ 79.260645060000172, 32.532858990000022 ], [ 79.126525878000166, 32.464309692000143 ], [ 79.104469299000129, 32.371711732000165 ], [ 79.049217225, 32.383369446000188 ], [ 78.993759154000031, 32.364299775000177 ], [ 78.977203368000119, 32.342769624000027 ], [ 78.914672851000034, 32.367092133000142 ], [ 78.874572754000098, 32.414619446000131 ], [ 78.817901611000082, 32.438739776000148 ], [ 78.776672364000149, 32.483928681000123 ], [ 78.763710022000055, 32.567810059000067 ], [ 78.780326843000125, 32.59832000800003 ], [ 78.763916016, 32.678489685000102 ], [ 78.74592590300017, 32.696079255000086 ], [ 78.654426575000173, 32.649230957000043 ], [ 78.625671386000079, 32.607051850000062 ], [ 78.544929505000084, 32.612308502000133 ], [ 78.485610962000067, 32.573219300000119 ], [ 78.452720642000088, 32.57950210700011 ], [ 78.412239074000126, 32.556060791000164 ], [ 78.405227661000026, 32.525230409000073 ], [ 78.416526795000038, 32.515098572000113 ], [ 78.382431030000077, 32.529651641000157 ], [ 78.310409545000027, 32.47592926100009 ], [ 78.27851867600009, 32.507930756000178 ], [ 78.299171448000038, 32.529270173000043 ], [ 78.296623229000147, 32.577510834000066 ], [ 78.331756591000101, 32.58486938599998 ], [ 78.388816834000124, 32.623229981000065 ], [ 78.37092590300017, 32.638160706000065 ], [ 78.363151549000122, 32.671989441000164 ], [ 78.377822877000085, 32.719310761000088 ], [ 78.371353149000186, 32.761829376000094 ], [ 78.34474182200006, 32.763431547999971 ], [ 78.28765869200015, 32.737831116000109 ], [ 78.295722962000184, 32.712810516000047 ], [ 78.236587524000072, 32.694671631000176 ], [ 78.212936401000093, 32.669189454000048 ], [ 78.183952332000104, 32.663871765000124 ], [ 78.152542115000131, 32.679019928 ], [ 78.121498108000026, 32.642349244000059 ], [ 78.093376160000105, 32.662151337000068 ], [ 78.065330506, 32.624279022999986 ], [ 78.033897400000171, 32.635608673000036 ], [ 78.035797120000154, 32.591949463000049 ], [ 77.982559205000143, 32.586528779000105 ], [ 77.952659606999987, 32.640399932999969 ], [ 77.903251649000083, 32.693099975000166 ], [ 77.918472290000068, 32.730789186000095 ], [ 77.915130615000066, 32.766059876000043 ], [ 77.880897521, 32.775409699000136 ], [ 77.847602844000107, 32.828411103000178 ], [ 77.767723083000192, 32.858470917000091 ], [ 77.76597595100003, 32.880290986000091 ], [ 77.790168762000064, 32.905429840000181 ], [ 77.712486267000088, 32.971801757000094 ], [ 77.703689575000169, 32.966129304000162 ], [ 77.67086029100011, 32.986270905000026 ], [ 77.653312684000127, 32.959400176000088 ], [ 77.612960816000054, 32.941410065000071 ], [ 77.583923340000183, 32.943801880000024 ], [ 77.519561767000141, 32.887378694000176 ], [ 77.484550476000095, 32.886539459000062 ], [ 77.455223084000181, 32.86174011200012 ], [ 77.43617248400011, 32.865890504000163 ], [ 77.425720214000194, 32.883670807000044 ], [ 77.387466432000167, 32.885520934000112 ], [ 77.359680175000165, 32.859580994 ], [ 77.356475830000193, 32.825271606000058 ], [ 77.327896118000183, 32.819881439000028 ], [ 77.305160523000097, 32.854049682000095 ], [ 77.281280518000187, 32.855651855000076 ], [ 77.258316040000068, 32.876361847000169 ], [ 77.23887634200014, 32.865161895000028 ], [ 77.22482299800015, 32.880020141000159 ], [ 77.228256226000155, 32.895568848000153 ], [ 77.188957215000073, 32.910060882000153 ], [ 77.173858642000027, 32.931129456000065 ], [ 77.146186828, 32.914020539000148 ], [ 77.13639068700013, 32.980461121000189 ], [ 77.103157044000056, 32.987201691000109 ], [ 77.074897765000173, 32.973350525000171 ], [ 77.03531646800019, 32.999679566000054 ], [ 76.998466491000045, 32.989570617000027 ], [ 76.917510985000092, 33.033679962000178 ], [ 76.877197266000053, 33.114238739000143 ], [ 76.844436646000133, 33.112300873000038 ], [ 76.803459168000188, 33.155689238999969 ], [ 76.803306579000093, 33.193851472000176 ], [ 76.818588257000101, 33.20557022100013 ], [ 76.811080933000142, 33.215461732000108 ], [ 76.828346252000188, 33.24639892700003 ], [ 76.802757263000103, 33.236660003000168 ], [ 76.778266908000148, 33.255710601000089 ], [ 76.755340575999981, 33.249488832000168 ], [ 76.759880066000107, 33.224670411000091 ], [ 76.729866027000128, 33.179431916000055 ], [ 76.627006530000131, 33.162479400000052 ], [ 76.599586487000067, 33.173759461000145 ], [ 76.581901551000158, 33.206871032000038 ], [ 76.547477723000156, 33.209529877000193 ], [ 76.469062806000125, 33.179908752000074 ], [ 76.391677857000047, 33.186691284000119 ], [ 76.275001525000164, 33.103988649000087 ], [ 76.237289428000111, 33.03258895800019 ], [ 76.09368896400008, 33.004829407000045 ], [ 76.082206726000038, 32.96879959100005 ], [ 76.028076172000169, 32.937438965000126 ], [ 76.017189026000153, 32.917640687000187 ], [ 75.945770263000043, 32.885139466000169 ], [ 75.867118835000099, 32.930210113000101 ], [ 75.807952882000052, 32.929790496000123 ], [ 75.789833068000064, 32.88140106100019 ], [ 75.821060181000178, 32.842761993000067 ], [ 75.873786927000083, 32.813690186000144 ], [ 75.912483215000123, 32.757122040000183 ], [ 75.896575927000015, 32.691059112000119 ], [ 75.923622132000105, 32.64369964500014 ], [ 75.804336548000151, 32.487060547000056 ], [ 75.745925903000057, 32.476348878000124 ], [ 75.711532594000175, 32.418159484000114 ], [ 75.684959411000193, 32.400619507000044 ], [ 75.545150757000158, 32.360820770000032 ], [ 75.524269104000041, 32.340049744000055 ], [ 75.501235962000123, 32.276050567000141 ], [ 75.490226746000133, 32.310211183000035 ], [ 75.467758178000167, 32.317562103000057 ], [ 75.47223663300008, 32.340049744000055 ], [ 75.41499328600014, 32.324272157000109 ], [ 75.33177947900009, 32.344699859000116 ], [ 75.322326661000091, 32.337326051 ], [ 75.231529236000142, 32.399188995000031 ], [ 75.195472717000087, 32.40531158400006 ], [ 75.188529969000115, 32.42448043800016 ], [ 75.147636412999987, 32.41337967000004 ], [ 75.127929687000176, 32.421169282000108 ], [ 75.109786987000064, 32.454040528000121 ], [ 75.071609497000111, 32.482963563 ] ] ] } },
        # create a subindex column that is useful
        # in the GADM processing of sub-national zones
        geodf_temp["GADM_ID"] = geodf_temp[f"GID_{layer_id}"]

        # append geodataframes
        geodf_list.append(geodf_temp)

    geodf_GADM = gpd.GeoDataFrame(pd.concat(geodf_list, ignore_index=True))
    geodf_GADM.set_crs(geo_crs)

    return geodf_GADM


def _simplify_polys(polys, minarea=0.0001, tolerance=0.008, filterremote=False):
    "Function to simplify the shape polygons"
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon(
                [
                    p
                    for p in takewhile(lambda p: p.area > minarea, polys)
                    if not filterremote or (mainpoly.distance(p) < mainlength)
                ]
            )
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)


def countries(countries, geo_crs, update=False, out_logging=False):
    "Create country shapes"

    if out_logging:
        _logger.info("Stage 1 of 4: Create country shapes")

    # download data if needed and get the layer id 0, corresponding to the countries
    df_countries = get_GADM_layer(countries, 0, geo_crs, update, out_logging)

    # select and rename columns
    df_countries = df_countries[["GID_0", "geometry"]].copy()
    df_countries.rename(columns={"GID_0": "name"}, inplace=True)

    # set index and simplify polygons
    ret_df = df_countries.set_index("name")["geometry"].map(_simplify_polys)

    return ret_df


def country_cover(country_shapes, eez_shapes=None, out_logging=False, distance=0.1):

    if out_logging:
        _logger.info("Stage 3 of 4: Merge country shapes to create continent shape")

    shapes = country_shapes.apply(lambda x: x.buffer(distance))
    shapes_list = list(shapes)
    if eez_shapes is not None:
        shapes_list += list(eez_shapes)

    africa_shape = unary_union(shapes_list)

    return africa_shape


def save_to_geojson(df, fn):
    if os.path.exists(fn):
        os.unlink(fn)  # remove file if it exists
    if not isinstance(df, gpd.GeoDataFrame):
        df = gpd.GeoDataFrame(dict(geometry=df))

    # save file if the GeoDataFrame is non-empty
    if df.shape[0] > 0:
        df = df.reset_index()
        schema = {**gpd.io.file.infer_schema(df), "geometry": "Unknown"}
        df.to_file(fn, driver="GeoJSON", schema=schema)
    else:
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass


def load_EEZ(countries_codes, geo_crs, EEZ_gpkg="./data/eez/eez_v11.gpkg"):
    """
    Function to load the database of the Exclusive Economic Zones.
    The dataset shall be downloaded independently by the user (see guide) or
    together with pypsa-earth package.
    """
    if not os.path.exists(EEZ_gpkg):
        raise Exception(
            f"File EEZ {EEZ_gpkg} not found, please download it from https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip and copy it in {os.path.dirname(EEZ_gpkg)}"
        )

    geodf_EEZ = gpd.read_file(EEZ_gpkg).to_crs(geo_crs)
    geodf_EEZ.dropna(axis=0, how="any", subset=["ISO_TER1"], inplace=True)
    # [["ISO_TER1", "TERRITORY1", "ISO_SOV1", "ISO_SOV2", "ISO_SOV3", "geometry"]]
    geodf_EEZ = geodf_EEZ[["ISO_TER1", "geometry"]]
    selected_countries_codes_3D = [
        two_2_three_digits_country(x) for x in countries_codes
    ]
    geodf_EEZ = geodf_EEZ[
        [any([x in selected_countries_codes_3D]) for x in geodf_EEZ["ISO_TER1"]]
    ]
    geodf_EEZ["ISO_TER1"] = geodf_EEZ["ISO_TER1"].map(
        lambda x: three_2_two_digits_country(x)
    )
    geodf_EEZ.reset_index(drop=True, inplace=True)

    geodf_EEZ.rename(columns={"ISO_TER1": "name"}, inplace=True)

    return geodf_EEZ


def eez(countries, geo_crs, country_shapes, EEZ_gpkg, out_logging=False, distance=0.01):
    """
    Creates offshore shapes by
    - buffer smooth countryshape (=offset country shape)
    - and differ that with the offshore shape
    Leads to for instance a 100m non-build coastline

    """

    if out_logging:
        _logger.info("Stage 2 of 4: Create offshore shapes")

    # load data
    df_eez = load_EEZ(countries, geo_crs, EEZ_gpkg)

    ret_df = df_eez[["name", "geometry"]]
    # create unique shape if country is described by multiple shapes
    for c_code in countries:
        selection = ret_df.name == c_code
        n_offshore_shapes = selection.sum()

        if n_offshore_shapes > 1:
            # when multiple shapes per country, then merge polygons
            geom = ret_df[selection].geometry.unary_union
            ret_df.drop(ret_df[selection].index, inplace=True)
            ret_df = ret_df.append(
                {"name": c_code, "geometry": geom}, ignore_index=True
            )

    ret_df = ret_df.set_index("name")["geometry"].map(
        lambda x: _simplify_polys(x, minarea=0.001, tolerance=0.0001)
    )

    ret_df = ret_df.apply(lambda x: make_valid(x))
    country_shapes = country_shapes.apply(lambda x: make_valid(x))

    country_shapes_with_buffer = country_shapes.buffer(distance)
    ret_df_new = ret_df.difference(country_shapes_with_buffer)

    # repeat to simplify after the buffer correction
    ret_df_new = ret_df_new.map(
        lambda x: x
        if x is None
        else _simplify_polys(x, minarea=0.001, tolerance=0.0001)
    )
    ret_df_new = ret_df_new.apply(lambda x: x if x is None else make_valid(x))

    # Drops empty geometry
    ret_df = ret_df_new.dropna()

    return ret_df


def download_WorldPop(
    country_code,
    worldpop_method,
    year=2020,
    update=False,
    out_logging=False,
    size_min=300,
):
    """
    Download Worldpop using either the standard method or the API method.
        Parameters
        ----------
        worldpop_method: str
             worldpop_method = "api" will use the API method to access the WorldPop 100mx100m dataset.  worldpop_method = "standard" will use the standard method to access the WorldPop 1KMx1KM dataset.
        country_code : str
            Two letter country codes of the downloaded files.
            Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
        year : int
            Year of the data to download
        update : bool
            Update = true, forces re-download of files
        size_min : int
            Minimum size of each file to download
    """
    if worldpop_method == "api":
        return download_WorldPop_API(country_code, year, update, out_logging, size_min)

    elif worldpop_method == "standard":
        return download_WorldPop_standard(
            country_code, year, update, out_logging, size_min
        )


def download_WorldPop_standard(
    country_code,
    year=2020,
    update=False,
    out_logging=False,
    size_min=300,
):
    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    WorldPop_filename : str
        Name of the file
    """
    if out_logging:
        _logger.info("Stage 3/4: Download WorldPop datasets")

    if country_code == "XK":
        WorldPop_filename = f"srb_ppp_{year}_UNadj_constrained.tif"
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/SRB/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/SRB/{WorldPop_filename}",
        ]
    else:
        WorldPop_filename = f"{two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
        # Urls used to possibly download the file
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        ]

    WorldPop_inputfile = os.path.join(
        os.getcwd(), "data", "WorldPop", WorldPop_filename
    )  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            _logger.warning(
                f"Stage 4/4: {WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break
        if not loaded:
            _logger.error(f"Stage 4/4: Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def download_WorldPop_API(
    country_code, year=2020, update=False, out_logging=False, size_min=300
):
    """
    Download tiff file for each country code using the api method from worldpop API with 100mx100m resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    WorldPop_filename : str
        Name of the file
    """
    if out_logging:
        _logger.info("Stage 3/4: Download WorldPop datasets")

    WorldPop_filename = f"{two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
    # Request to get the file
    WorldPop_inputfile = os.path.join(
        os.getcwd(), "data", "WorldPop", WorldPop_filename
    )  # Input filepath tif
    os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)
    year_api = int(str(year)[2:])
    loaded = False
    WorldPop_api_urls = [
        f"https://www.worldpop.org/rest/data/pop/wpgp?iso3={two_2_three_digits_country(country_code)}",
    ]
    for WorldPop_api_url in WorldPop_api_urls:
        with requests.get(WorldPop_api_url, stream=True) as r:
            WorldPop_tif_url = r.json()["data"][year_api]["files"][0]

        with requests.get(WorldPop_tif_url, stream=True) as r:
            with open(WorldPop_inputfile, "wb") as f:
                if float(r.headers["Content-length"]) > size_min:
                    shutil.copyfileobj(r.raw, f)
                    loaded = True
                    break
    if not loaded:
        _logger.error(f"Stage 4/4: Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def convert_GDP(name_file_nc, year=2015, out_logging=False):
    """
    Function to convert the nc database of the GDP to tif, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or toghether with pypsa-earth package.
    """

    if out_logging:
        _logger.info("Stage 4/4: Access to GDP raster data")

    # tif namefile
    name_file_tif = name_file_nc[:-2] + "tif"

    # path of the nc file
    GDP_nc = os.path.join(os.getcwd(), "data", "GDP", name_file_nc)  # Input filepath nc

    # path of the tif file
    GDP_tif = os.path.join(
        os.getcwd(), "data", "GDP", name_file_tif
    )  # Input filepath nc

    # Check if file exists, otherwise throw exception
    if not os.path.exists(GDP_nc):
        raise Exception(
            f"File {name_file_nc} not found, please download it from https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0 and copy it in {os.path.dirname(GDP_nc)}"
        )

    # open nc dataset
    GDP_dataset = xr.open_dataset(GDP_nc)

    # get the requested year of data or its closest one
    list_years = GDP_dataset["time"]
    if year not in list_years:
        if out_logging:
            _logger.warning(
                f"Stage 3/4 GDP data of year {year} not found, selected the most recent data ({int(list_years[-1])})"
            )
        year = float(list_years[-1])

    # subset of the database and conversion to dataframe
    GDP_dataset = GDP_dataset.sel(time=year).drop("time")
    GDP_dataset.rio.to_raster(GDP_tif)

    return GDP_tif, name_file_tif


def load_GDP(
    countries_codes,
    year=2015,
    update=False,
    out_logging=False,
    name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
):
    """
    Function to load the database of the GDP, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or toghether with pypsa-earth package.
    """

    if out_logging:
        _logger.info("Stage 4/4: Access to GDP raster data")

    # path of the nc file
    name_file_tif = name_file_nc[:-2] + "tif"
    GDP_tif = os.path.join(
        os.getcwd(), "data", "GDP", name_file_tif
    )  # Input filepath tif

    if update | (not os.path.exists(GDP_tif)):
        if out_logging:
            _logger.warning(
                f"Stage 4/4: File {name_file_tif} not found, the file will be produced by processing {name_file_nc}"
            )
        convert_GDP(name_file_nc, year, out_logging)

    return GDP_tif, name_file_tif


def generalized_mask(src, geom, **kwargs):
    "Generalize mask function to account for Polygon and MultiPolygon"
    if geom.geom_type == "Polygon":
        return mask(src, [geom], **kwargs)
    elif geom.geom_type == "MultiPolygon":
        return mask(src, geom.geoms, **kwargs)
    else:
        return mask(src, geom, **kwargs)


def _sum_raster_over_mask(shape, img):
    """
    Function to sum the raster value within a shape
    """
    # select the desired area of the raster corresponding to each polygon
    # Approximation: the population is measured including the pixels
    #   where the border of the shape lays. This leads to slightly overestimate
    #   the output, but the error is limited and it enables halving the
    #   computational time
    out_image, out_transform = generalized_mask(
        img, shape, all_touched=True, invert=False, nodata=0.0
    )
    # calculate total output in the selected geometry
    out_image[np.isnan(out_image)] = 0
    out_sum = out_image.sum()
    # out_sum = out_image.sum()/2 + out_image_int.sum()/2

    return out_sum


def add_gdp_data(
    df_gadm,
    year=2020,
    update=False,
    out_logging=False,
    name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
    nprocesses=2,
    disable_progressbar=False,
):
    """
    Function to add gdp data to arbitrary number of shapes in a country

    Inputs:
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Outputs:
    --------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["gdp"]
    """
    if out_logging:
        _logger.info("Stage 4/4: Add gdp data to GADM GeoDataFrame")

    # initialize new gdp column
    df_gadm["gdp"] = 0.0

    GDP_tif, name_tif = load_GDP(year, update, out_logging, name_file_nc)

    with rasterio.open(GDP_tif) as src:
        # resample data to target shape
        tqdm_kwargs = dict(
            ascii=False,
            unit=" geometries",
            total=df_gadm.shape[0],
            desc="Compute GDP ",
        )
        for i in tqdm(df_gadm.index, **tqdm_kwargs):
            df_gadm.loc[i, "gdp"] = _sum_raster_over_mask(df_gadm.geometry.loc[i], src)
    return df_gadm

    # for index, row in tqdm(df_gadm.iterrows(), index=df_gadm.shape[0]):
    #     # select the desired area of the raster corresponding to each polygon
    #     # Approximation: the gdp is measured excluding the pixels
    #     #   where the border of the shape lays. This may affect the computation
    #     #   but it is conservative and avoids considering multiple times the same
    #     #   pixels
    #     out_image, out_transform = generalized_mask(src,
    #                                                 row["geometry"],
    #                                                 all_touched=True,
    #                                                 invert=False,
    #                                                 nodata=0.0)
    #     # out_image_int, out_transform = mask(src,
    #     #                                row["geometry"],
    #     #                                all_touched=False,
    #     #                                invert=False,
    #     #                                nodata=0.0)

    #     # calculate total gdp in the selected geometry
    #     gdp_by_geom = np.nansum(out_image)
    #     # gdp_by_geom = out_image.sum()/2 + out_image_int.sum()/2

    #     if out_logging == True:
    #         _logger.info("Stage 4/4 GDP: shape: " + str(index) +
    #                      " out of " + str(df_gadm.shape[0]))

    #     # update the gdp data in the dataset
    #     df_gadm.loc[index, "gdp"] = gdp_by_geom


def _init_process_pop(df_gadm_, year_, worldpop_method_):
    global df_gadm, year, worldpop_method
    df_gadm, year, worldpop_method = df_gadm_, year_, worldpop_method_


def _process_func_pop(c_code):

    # get subset by country code
    country_rows = df_gadm.loc[df_gadm["country"] == c_code].copy()

    # get worldpop image
    WorldPop_inputfile, WorldPop_filename = download_WorldPop(
        c_code, worldpop_method, year, False, False
    )

    with rasterio.open(WorldPop_inputfile) as src:

        for i in country_rows.index:
            country_rows.loc[i, "pop"] = _sum_raster_over_mask(
                country_rows.geometry.loc[i], src
            )

    return country_rows


def add_population_data(
    df_gadm,
    country_codes,
    worldpop_method,
    year=2020,
    update=False,
    out_logging=False,
    nprocesses=2,
    disable_progressbar=False,
):
    """
    Function to add population data to arbitrary number of shapes in a country

    Inputs:
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Outputs:
    --------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["pop"]
    """

    if out_logging:
        _logger.info("Stage 4/4 POP: Add population data to GADM GeoDataFrame")

    # initialize new population column
    df_gadm["pop"] = 0.0

    tqdm_kwargs = dict(
        ascii=False,
        unit=" countries",
        # total=len(country_codes),
        desc="Compute population ",
    )
    if (nprocesses is None) or (nprocesses == 1):

        with tqdm(total=df_gadm.shape[0], **tqdm_kwargs) as pbar:

            for c_code in country_codes:
                # get subset by country code
                country_rows = df_gadm.loc[df_gadm["country"] == c_code]

                # get worldpop image
                WorldPop_inputfile, WorldPop_filename = download_WorldPop(
                    worldpop_method, c_code, year, update, out_logging
                )

                with rasterio.open(WorldPop_inputfile) as src:

                    for i, row in country_rows.iterrows():
                        df_gadm.loc[i, "pop"] = _sum_raster_over_mask(row.geometry, src)
                        pbar.update(1)

    else:

        kwargs = {
            "initializer": _init_process_pop,
            "initargs": (df_gadm, year, worldpop_method),
            "processes": nprocesses,
        }
        with mp.get_context("spawn").Pool(**kwargs) as pool:
            if disable_progressbar:
                _ = list(pool.map(_process_func_pop, country_codes))
                for elem in _:
                    df_gadm.loc[elem.index, "pop"] = elem["pop"]
            else:
                _ = list(
                    tqdm(
                        pool.imap(_process_func_pop, country_codes),
                        total=len(country_codes),
                        **tqdm_kwargs,
                    )
                )
                for elem in _:
                    df_gadm.loc[elem.index, "pop"] = elem["pop"]


def gadm(
    worldpop_method,
    gdp_method,
    countries,
    geo_crs,
    layer_id=2,
    update=False,
    out_logging=False,
    year=2020,
    nprocesses=None,
):

    if out_logging:
        _logger.info("Stage 4/4: Creation GADM GeoDataFrame")

    # download data if needed and get the desired layer_id
    df_gadm = get_GADM_layer(countries, layer_id, geo_crs, update)

    # select and rename columns
    df_gadm.rename(columns={"GID_0": "country"}, inplace=True)

    # drop useless columns
    df_gadm.drop(
        df_gadm.columns.difference(["country", "GADM_ID", "geometry"]),
        axis=1,
        inplace=True,
    )

    if worldpop_method != False:
        # add the population data to the dataset
        add_population_data(
            df_gadm,
            countries,
            worldpop_method,
            year,
            update,
            out_logging,
            nprocesses=nprocesses,
        )

    if gdp_method != False:
        # add the gdp data to the dataset
        add_gdp_data(
            df_gadm,
            year,
            update,
            out_logging,
            name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
        )

    # set index and simplify polygons
    df_gadm.set_index("GADM_ID", inplace=True)
    df_gadm["geometry"] = df_gadm["geometry"].map(_simplify_polys)

    return df_gadm


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_shapes")
        sets_path_to_root("pypsa-earth")
    configure_logging(snakemake)

    out = snakemake.output

    countries_list = snakemake.config["countries"]
    layer_id = snakemake.config["build_shape_options"]["gadm_layer_id"]
    update = snakemake.config["build_shape_options"]["update_file"]
    out_logging = snakemake.config["build_shape_options"]["out_logging"]
    year = snakemake.config["build_shape_options"]["year"]
    nprocesses = snakemake.config["build_shape_options"]["nprocesses"]
    EEZ_gpkg = snakemake.input["eez"]
    worldpop_method = snakemake.config["build_shape_options"]["worldpop_method"]
    gdp_method = snakemake.config["build_shape_options"]["gdp_method"]
    geo_crs = snakemake.config["crs"]["geo_crs"]
    distance_crs = snakemake.config["crs"]["distance_crs"]

    country_shapes = countries(countries_list, geo_crs, update, out_logging)

    country_shapes.reset_index().to_file(snakemake.output.country_shapes)

    offshore_shapes = eez(
        countries_list, geo_crs, country_shapes, EEZ_gpkg, out_logging
    )

    offshore_shapes.reset_index().to_file(snakemake.output.offshore_shapes)

    africa_shape = gpd.GeoDataFrame(
        geometry=[country_cover(country_shapes, offshore_shapes.geometry)]
    )
    africa_shape.reset_index().to_file(snakemake.output.africa_shape)

    gadm_shapes = gadm(
        worldpop_method,
        gdp_method,
        countries_list,
        geo_crs,
        layer_id,
        update,
        out_logging,
        year,
        nprocesses=nprocesses,
    )
    save_to_geojson(gadm_shapes, out.gadm_shapes)
