
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

SDIR = config['summary_dir'] + '/' + config['run']
RDIR = config['results_dir'] + config['run']
CDIR = config['costs_dir']


rule prepare_sector_network:
    input:
        network='networks/elec_s{simpl}_{clusters}.nc',
        costs=CDIR + "costs_2030.csv",
        h2_cavern="data/hydrogen_salt_cavern_potentials.csv"
        

    output: RDIR + '/prenetworks/elec_s{simpl}_{clusters}.nc'

    script: "scripts/prepare_sector_network.py"