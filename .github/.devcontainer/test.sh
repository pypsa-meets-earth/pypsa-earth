#!/bin/bash

cp config.tutorial.yaml config.yaml

snakemake -j 1 solve_all_networks
