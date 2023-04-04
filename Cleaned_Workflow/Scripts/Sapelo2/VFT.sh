#!/bin/bash
#SBATCH --job-name=VFT                                    # Job name
#SBATCH --partition=batch                                 # Partition (queue) name
#SBATCH --ntasks=1			                                  # Single task job                           #
#SBATCH --mem=10gb                                      # Total memory for job
#SBATCH --time=24:00:00  		                              # Time limit hrs:min:sec
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

# source activate R
# # R --no-save < /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Scripts/sleuth_timeseries.R
# # R --no-save < /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Scripts/sleuth.R
# source deactivate R


# ####################################################################
# # 4) MASIGPRO -- ran on local comp.
# ####################################################################

# # mkdir $OUTDIR/masigpro
# R --no-save < /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Scripts/MaSigPro.R


####################################################################
# 5) BLAST LISTS
####################################################################


ml BLAST+/2.12.0-gompi-2020b
ml HpcGridRunner/1.0.2
# ml seqkit/0.16.1
#
# mkdir $OUTDIR/BLAST
# mkdir $OUTDIR/BLAST/db
# mkdir $OUTDIR/BLAST/DmProteinsBLAST
# mkdir $OUTDIR/BLAST/DmProteinsBLAST/

#maketairdb
# curl -s https://www.arabidopsis.org/download_files/Proteins/Araport11_protein_lists/Araport11_pep_20220914.gz | gunzip -c > $OUTDIR/BLAST/db/Araport11_pep_20220914
#
# makeblastdb -in $OUTDIR/BLAST/db/Araport11_pep_20220914 -parse_seqids -blastdb_version 5 -title "db_Araport11_pep_20220914" -dbtype prot

#translate transcripts to proteins
# seqkit translate $OUTDIR/VFT/Dm_transcripts.fa > $OUTDIR/BLAST/Dm_proteins.fa

#blast Dmproteins file
# hpc_FASTA_GridRunner.pl --grid_conf=$HPCGR_CONF_DIR/sapelo_1c_3g.conf \
# --cmd_template "blastp \
# -query __QUERY_FILE__ \
# -db /db/ncbiblast/20221102/swissprot \
# -max_target_seqs 1 \
# -max_hsps 1 \
# -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
# -evalue 1e-20 -num_threads 1" \
# --seqs_per_bin 100 \
# --query_fasta $OUTDIR/BLAST/Dm_proteins.fa \
# --out_dir $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprot
#
# cat $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprot/*/*.OUT > $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprotBLASTconcat.txt

# hpc_FASTA_GridRunner.pl --grid_conf=$HPCGR_CONF_DIR/sapelo_1c_3g.conf \
# --cmd_template "blastp \
# -query __QUERY_FILE__ \
# -db $OUTDIR/BLAST/db/Araport11_pep_20220914 \
# -max_target_seqs 1 \
# -max_hsps 1 \
# -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
# -evalue 1e-20 -num_threads 1" \
# --seqs_per_bin 100 \
# --query_fasta $OUTDIR/BLAST/Dm_proteins.fa \
# --out_dir $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdb

# cat $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdb/*/*.OUT > $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdbBLASTconcat.txt

#cleans CSV output into individual files for each cluster with 1 gene name per line for grep BLAST

# mkdir $OUTDIR/BLAST/MaSigProClusters
# mkdir $OUTDIR/BLAST/MaSigProClusters/NoPreyTimeSeries
# mkdir $OUTDIR/BLAST/MaSigProClusters/PreyTimeSeries
# mkdir $OUTDIR/BLAST/MaSigProClusters/BothTimeSeries
#
# #no prey
#
# for k in 1 2; do
#    sed '1d' /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Results/maSigPro/NoPreyTimeSeries/MaSigProNoPreyCluster${k}.csv | awk -F "," '{ print $2}' | tr -d '"' > $OUTDIR/BLAST/MaSigProClusters/NoPreyTimeSeries/MaSigProNoPreyCluster${k}.txt
# done
#
# #prey
#
# for k in 1 2 3 4 5; do
#    sed '1d' /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Results/maSigPro/PreyTimeSeries/MaSigProPreyCluster${k}.csv |    > $OUTDIR/BLAST/MaSigProClusters/PreyTimeSeries/PreyCluster${k}.txt
# done
#
# #shared time point
# for k in 1 2 3 4 5 6; do
#    sed '1d' /home/srb67793/2022VenusFlyTrap/Cleaned_Workflow/Results/maSigPro/SharedTimeSeries/MaSigProSharedCluster${k}.csv | awk -F "," '{ print $2}' | tr -d '"' > $OUTDIR/BLAST/MaSigProClusters/BothTimeSeries/SharedCluster${k}.txt
# done
#
# mkdir $OUTDIR/BLAST/Results
#
# #grep blast swissprot
# for file in $OUTDIR/BLAST/MaSigProClusters/PreyTimeSeries/PreyCluster*; do
#   base=$(basename ${file} .txt)
#   grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprotBLASTconcat.txt > $OUTDIR/BLAST/Results/${base}_swissprotBLAST.txt
# done
#
# for file in $OUTDIR/BLAST/MaSigProClusters/BothTimeSeries/SharedCluster*; do
#   base=$(basename ${file} .txt)
#   grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprotBLASTconcat.txt > $OUTDIR/BLAST/Results/${base}_swissprotBLAST.txt
# done
#
# for file in $OUTDIR/BLAST/MaSigProClusters/NoPreyTimeSeries/MaSigProNoPreyCluster*; do
#   base=$(basename ${file} .txt)
#   grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsSwissprotBLASTconcat.txt > $OUTDIR/BLAST/Results/${base}_swissprotBLAST.txt
# done

#grep blast tairdb
# for file in $OUTDIR/BLAST/MaSigProClusters/PreyTimeSeries/PreyCluster*; do
#   base=$(basename ${file} .txt)
#   grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdbBLASTconcat.txt > $OUTDIR/BLAST/Results/${base}_tairBLAST.txt
# done
#
# for file in $OUTDIR/BLAST/MaSigProClusters/BothTimeSeries/SharedCluster*; do
#   base=$(basename ${file} .txt)
#   grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdbBLASTconcat.txt > $OUTDIR/BLAST/Results/${base}_tairBLAST.txt
# done
#

#BLAST sleuth LISTS

mkdir $OUTDIR/BLAST/Sleuth
sed '1d' /scratch/srb67793/2022VenusFlyTrap/sleuth/petioleVstraps_sleuth_q_0.05.csv | awk -F "," '{ print $1}' | tr -d '"' > $OUTDIR/BLAST/Sleuth/SleuthPetiolesVsTraps.txt

sed '1d' /scratch/srb67793/2022VenusFlyTrap/sleuth/SleuthPreyNoPrey60minResults_q_0.05.csv | awk -F "," '{ print $1}' | tr -d '"' > $OUTDIR/BLAST/Sleuth/SleuthPreyNoPrey60min.txt

sed '1d' /scratch/srb67793/2022VenusFlyTrap/sleuth/SleuthPreyNoPrey1440min_q_0.05.csv | awk -F "," '{ print $1}' | tr -d '"' > $OUTDIR/BLAST/Sleuth/SleuthPreyNoPrey1440min.txt

for file in $OUTDIR/BLAST/Sleuth/Sleuth*; do
  base=$(basename ${file} .txt)
  grep -f ${file} $OUTDIR/BLAST/DmProteinsBLAST/DmProteinsTairdbBLASTconcat.txt > $OUTDIR/BLAST/Sleuth/${base}_tairBLAST.txt
done

for i in $OUTDIR/BLAST/Sleuth/*BLAST*; do
  base=$(basename ${i} _tairBLAST.txt)
  awk '{print $2}' ${i} | cut -c 1-9 > ${base}_ATnames.txt
done


#on local computer#
#
# for i in /Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/BLAST/TairBLAST/*.txt; do
#   base=$(basename ${i} _tairBLAST.txt)
#   awk '{print $2}' ${i} | cut -c 1-9 > ${base}_ATnames.txt
# done

#make background set of tair proteins
# awk '{print $2}' DmProteinsTairdbBLASTconcat.txt | cut -c 1-9 > AllATnames.txt
