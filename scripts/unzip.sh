#!/bin/bash
#SBATCH --job-name=unzip
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=50gb
#SBATCH --time=7:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.error

bzip2 -d *.bz2
