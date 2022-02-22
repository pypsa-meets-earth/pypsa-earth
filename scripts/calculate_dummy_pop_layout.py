import os

import numpy as np
import pandas as pd
import pypsa
from helpers import mock_snakemake

def calculate_dummy_pop_layout(n, inhabitants):
    """
    Function to create dummy pop layout (clustered)
    """
    # Get pop_layout of PyPSA-Eur-Sec (to get column names)
    pop_layout_eur = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    # Create new pop_layout for morocco
    pop_layout_dummy = pd.DataFrame(columns=pop_layout_eur.columns, index = n.buses.index)

    # Calculate the rural and urban population for morocco, assuming an equal distribution between nodes and urban/rural
    pop_per_node = inhabitants / len(pop_layout_dummy.index) / 2

    # Add population to pop_layout_dummy
    pop_layout_dummy['urban'] = pop_per_node
    pop_layout_dummy['rural'] = pop_per_node
    pop_layout_dummy['total'] = pop_layout_dummy['urban'] + pop_layout_dummy['rural']

    # Add country and fraction of nodes population to countries total population
    pop_layout_dummy['ct'] = 'MA'
    pop_layout_dummy['fraction'] = 1 / len(pop_layout_dummy[pop_layout_dummy['ct'] == 'MA'].index)

    return pop_layout_dummy


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake("calculate_dummy_pop_layout",
                                   simpl="",
                                   clusters="4")

    n = pypsa.Network(snakemake.input.network)

    # Get pop_layout 
    inhabitants = 37000 # Value for Morocco in thousands
    pop_layout = calculate_dummy_pop_layout(n, inhabitants)

    # Save pop_layout
    pop_layout.to_csv(snakemake.output[0])
