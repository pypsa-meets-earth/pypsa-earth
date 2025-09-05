# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pandas as pd

REF_OBJ_DF = pd.DataFrame(
    {
        "folder": ["custom", "tutorial"],
        "file": [
            "elec_s_6_ec_lcopt_Co2L-4H_python.log",
            "elec_s_6_ec_lcopt_Co2L-4H_python.log",
        ],
        "objective": [2548130080.0, 2823517462.0],
    }
)


def extract_obj_for_path(dir, fl, ref_df=REF_OBJ_DF):
    mask = (ref_df["folder"] == dir) & (ref_df["file"] == fl)

    obj_value = ref_df["objective"][mask]

    if obj_value.empty:
        obj_value = "NA"
    else:
        obj_value = obj_value.values[0]

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
