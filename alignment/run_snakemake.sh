#!/bin/bash

#set +eu
#source /n/macklis_lab/users/priyav/environments/conda/envs/rnaseq_env_R3_6_3/bin/activate
snakemake --unlock --cores 1
#snakemake --cleanup-metadata /n/holyscratch01/macklis_lab/priyav/lab_notebook/amalgam_rnaseq/merged_fastq/CPN*
snakemake --rerun-incomplete --jobs 999 --cluster "sbatch -p shared -c {resources.cores}  -t {resources.time_min} --mem={resources.memory_mb} -o logs/slurm/%j.out -e logs/slurm/%j.err"
