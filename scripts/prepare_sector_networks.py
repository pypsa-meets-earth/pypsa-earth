import pypsa
import os

if __name__ == "__main__":

    # TODO add mock_snakemake func

    os.path.dirname(os.getcwd())
    # TODO fetch from config
    n = pypsa.Network(
        "{}/networks/elec_s_4.nc".format(os.path.dirname(os.getcwd())))

    n.buses

    # TODO logging

    # TODO fetch options from the config file

    # TODO define spatial (for biomass and co2)

    # TODO changes in case of myopic oversight

    # TODO add co2 tracking function 

    # TODO add generation

    # TODO add storage  HERE THE H2 CARRIER IS ADDED IN PYPSA-EUR-SEC

    # TODO add options as in PyPSA-EUR-SEC