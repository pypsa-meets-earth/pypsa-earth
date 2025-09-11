# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pandas as pd

REF_OBJ_DF = pd.DataFrame(
    {
        "folder": [
            "NG",
            "custom",
            "landlock",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "monte-carlo",
            "tutorial",
        ],
        "file": [
            "elec_s_5_ec_lcopt_Co2L-4H_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m0_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m1_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m2_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m3_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m4_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m5_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m6_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m7_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_m8_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_python.log",
        ],
        "objective": [
            2665231623.0,
            2548087889.0,
            1098081954.0,
            1825115005.0,
            10076445520.0,
            873800591.8,
            1848771010.0,
            45024568060.0,
            2756737870.0,
            6436852411.0,
            1447523254.0,
            406621644.0,
            2823517462.0,
        ],
    }
)

SCALE = 1_000_000


def extract_obj_for_path(dir, fl, ref_df=REF_OBJ_DF):
    mask = (ref_df["folder"] == dir) & (ref_df["file"] == fl)

    obj_value = ref_df["objective"][mask]

    if obj_value.empty:
        obj_value = "NA"
    else:
        obj_value = round(obj_value.values[0] / SCALE, 2)

    return obj_value


# print("test 1")
# print(
#    extract_obj_for_path(
#        "custom",
#        "elec_s_6_ec_lcopt_Co2L-4H_python.log",
#        ref_df=REF_OBJ_DF
#    )
# )
#
# print("test 2")
# print(
#    extract_obj_for_path(
#        "custom2",
#        "elec_s_6_ec_lcopt_Co2L-4H_python.log",
#        ref_df=REF_OBJ_DF
#    )
# )
