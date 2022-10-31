#!/bin/bash
#SBATCH --partition batch
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=3gb

awk -F "," '{ print $1}' geneclusters2.csv > traps_no_prey_names.txt
grep -f traps_no_prey_names.txt geneclusters8.csv > overlapgeneclusters.txt
awk -F "," '{ print $1}' overlapgeneclusters.txt > overlapgenenames.txt
grep -f overlapgenenames.txt geneclusters2.csv > overlapgenestrapsnopreyclusters.txt
sort -t"," -k1 overlapgeneclusters.txt > sortedoverlaptrapswithpreygenes.txt
sort -t"," -k1 overlapgenestrapsnopreyclusters.txt > sortedoverlaptrapsnopreygenes.txt
paste sortedoverlaptrapsnopreygenes.txt sortedoverlaptrapswithpreygenes.txt > overlap_comparison.txt
grep -P 1'\t' overlap_comparison.txt > overlapcluster1.txt
grep -P 2'\t' overlap_comparison.txt > overlapcluster2.txt
