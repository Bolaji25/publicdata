#! /bin/bash

## Allocate resources
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=all

## job name
#SBATCH --job-name="kranzData"

source ${CONDA_ACTIVATE}  MC-HiC-env


Rscript ./kranz_data.R

