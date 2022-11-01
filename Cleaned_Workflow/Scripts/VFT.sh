#!/bin/bash
#SBATCH --job-name=VFT                                    # Job name
#SBATCH --partition=batch                                 # Partition (queue) name
#SBATCH --ntasks=1			                                  # Single task job
#SBATCH --cpus-per-task=8	                                # Number of cores per taskT
#SBATCH --mem=50gb                                       # Total memory for job
#SBATCH --time=36:00:00  		                              # Time limit hrs:min:sec
#SBATCH --output="/home/srb67793/2022VenusFlyTrap/log.%j" # Location of standard output and error log files
#SBATCH --mail-user=srb67793@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                          # Mail events (BEGIN, END, FAIL, ALL)

##################################

#SUMMER BLANCO
#PHD STUDENT, PLANT BIOLOGY
#LEEBENS-MACK & CHANG LABS
#UNIVERSITY OF GEORGIA

#Venus Fly Trap RNAseq
#SAPELO2

##################################

#THIS SCRIPT:
#1) TRIMS VFT ILLUMINA SHORT READS WITH TRIMMMOMATIC -- DONE
#2) QC'S VFT ILLUMINA SHORT READS WITH FASTQC & MULTIQC -- DONE
#3) MAKES INDEX WITH KALLISTO INDEX -- DONE
#4) QUANTIFIES TRANSCRIPT ABUNDANCE WITH KALLISTO QUANT - DONE
#5) DIFFERENTIAL GENEEXPRESSION ANALYSIS WITH SLEUTH

##################################
# SET UP
##################################

#set output directory variable
OUTDIR="/scratch/srb67793/2022VenusFlyTrap"

#if output directory doesn't exist, create it

# if [ ! -d $OUTDIR ]
# then
#     mkdir -p $OUTDIR
# fi

##################################
#  DO IN AN INTERACTIVE SESSION
#################################
## 3on xfer node
## qlogin
## scp -r /project/jlmlab/VFT.tar.gz $OUTDIR

## tar -xf $OUTDIR/VFT.tar.gz
## mkdir $OUTDIR/rawreads
## mv $OUTDIR/VFT/*.bz2 $OUTDIR/rawreads

##################################
#  UNZIP READS
##################################

# bzip2 -d $OUTDIR/rawreads/*bz2

##################################
#  LOAD MODULES
##################################

# module load FastQC/0.11.9-Java-11
# module load MultiQC/1.8-foss-2019b-Python-3.7.4
# module load Trimmomatic/0.39-Java-1.8.0_144
# module load kallisto/0.46.1-foss-2019b

##################################
# 1) TRIMS VFT READS (more than 12 hrs)
##################################

#QC pre-trim with FASTQC & MultiQC
# mkdir $OUTDIR/FastQC
# mkdir $OUTDIR/FastQC/pretrim
# fastqc -o $OUTDIR/FastQC/pretrim/ $OUTDIR/rawreads/*
# multiqc $OUTDIR/FastQC/pretrim/*.zip -o $OUTDIR/FastQC/pretrim/


# mkdir $OUTDIR/trimmedreads
# for infile in $OUTDIR/rawreads/*1.fastq; do
#   base=$(basename ${infile} 1.fastq);
#   java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 \
#   ${infile} \
#   $OUTDIR/rawreads/${base}2.fastq \
#   $OUTDIR/trimmedreads/${base}1_paired.fastq \
#   $OUTDIR/trimmedreads/${base}1_unpaired.fastq \
#   $OUTDIR/trimmedreads/${base}2_paired.fastq \
#   $OUTDIR/trimmedreads/${base}2_unpaired.fastq \
#   ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 \
#   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# done

# mkdir $OUTDIR/trimmedreads/paired
# mkdir $OUTDIR/trimmedreads/unpaired
# mv $OUTDIR/trimmedreads/*_paired.fastq $OUTDIR/trimmedreads/paired
# mv $OUTDIR/trimmedreads/*_unpaired.fastq $OUTDIR/trimmedreads/unpaired

####################################################################
# 2) QC'S G MACULATUM ILLUMINA SHORT READS
####################################################################

# QC post-trim with FASTQC & MultiQC
# mkdir $OUTDIR/FastQC/trimmed
# fastqc -o $OUTDIR/FastQC/trimmed $OUTDIR/trimmedreads/paired/*
# multiqc $OUTDIR/FastQC/trimmed/*.zip -o $OUTDIR/FastQC/trimmed

####################################################################
# 3) MAKE KALLISTO INDEX
####################################################################
#
# mkdir $OUTDIR/kallisto
# kallisto index -i $OUTDIR/kallisto/VFT.idx $OUTDIR/VFT/Dm_transcripts.fa

####################################################################
# 4) KALLISTO QUANT
####################################################################

# ls $OUTDIR/trimmedreads/paired | awk -F _ '{print $1}' | uniq > $OUTDIR/trimmedreads/paired/librarynames.txt
# mkdir $OUTDIR/kallisto/quant
# cat $OUTDIR/trimmedreads/paired/librarynames.txt | while read i; do
#   kallisto quant -i $OUTDIR/kallisto/VFT.idx -o $OUTDIR/kallisto/quant/${i} -b 100 \
#   $OUTDIR/trimmedreads/paired/${i}_RNAseq*
# done
#
# ####################################################################
# # 4) SLEUTH
# ####################################################################

# mkdir $OUTDIR/sleuth

source activate R
R --no-save < /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Scripts/sleuth_timeseries.R
R --no-save < /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Scripts/sleuth.R
