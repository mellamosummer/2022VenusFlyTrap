#!/bin/bash
#SBATCH --partition batch
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=3gb

 grep -f data/traps_prey_names_sed.txt results/Dm_proteins_BLAST_combined/combined_blast_results_oneHSP.out > Mfuzz_traps_prey_BLAST_swissprot.txt
 grep -f data/traps_prey_names_sed.txt results/Dm_proteins_BLAST_combined/combined_blast_results_tairdb.out > Mfuzz_traps_prey_BLAST_tair.txt
 grep -f data/traps_no_prey_names_sed.txt results/Dm_proteins_BLAST_combined/combined_blast_results_oneHSP.out > Mfuzz_traps_no_prey_BLAST_swissprot.txt
 grep -f data/traps_no_prey_names_sed.txt results/Dm_proteins_BLAST_combined/combined_blast_results_tairdb.out > Mfuzz_traps_no_prey_BLAST_tair.txt

