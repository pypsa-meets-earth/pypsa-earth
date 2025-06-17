READ_ME GloFAS-ERA5

Script modificati sono: build_renewable_profiles, add_electricity
build_renewable_profiles ha flag (if/else) metodo GloFAS o alternativa versione corrente PyPSA
add_electricity ha flag che porta ad utilizzare attach_hydro_GloFAS o attach_hydro (uguale a prima) perché inserire if/else all'interno della funzione (quindi non ripeterla) per ora è complicato

Necessarie modifiche anche nel config (Config_SAPPGloFAS) che deve essere salvato come config.yaml e Snakefile

Modificato custom_powerplants.csv con già centrali GloHydroRes + centrali non hydro di PyPSA