#!/bin/bash

source activate env3.5

snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -c {cluster.cores}  -t {cluster.time} --mem={cluster.mem} -e {cluster.error} -o {cluster.output}"
