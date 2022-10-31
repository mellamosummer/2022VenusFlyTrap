#!/bin/bash
#SBATCH --partition batch
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=3gb

for file in data/[namelist.txt] ##change list depending on what you want to BLAST
  do
  base=$(basename ${file} .lst)
  grep -f ${file} results/Dm_proteins_BLAST_combined/combined_blast_results_oneHSP.out > results/DEGlistsWithBLASTdescriptions/swissprotBLAST/${base}_blast_description.out
  done
