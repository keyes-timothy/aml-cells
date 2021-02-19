#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem 50GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tkeyes@stanford.edu

ml R/4.0.2

Rscript feature_extraction_sherlock.R
