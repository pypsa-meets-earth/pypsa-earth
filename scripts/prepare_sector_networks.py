import pypsa



add_generation(n, costs):
    print(adding electricity generation)

if __name__ == "__main__":

    n = pypsa.Network("networks/elec_s_4.nc")

    n.buses()