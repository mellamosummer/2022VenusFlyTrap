- 1. copy over reads to scratch directory on xfer node
  - cp -r /project/jlmlab/VFT.tar.gz /scratch/srb67793/VFT_Analysis2/

- 2. untar and unzip files
  - tar -xf
  - qlogin
  - bzip2 -d /scratch/srb67793/VFT_Analysis2/VFT/Raw_reads/*.bz2

- 3. QC reads before trimming with FASTQC script
- 4. summarize FASTQC results with MultiQC
- 5. Trim reads using trimmomatic script
- 6. QC reads after trimming with FASTQC script
- 7. summarize FASTQC results with MultiQC
- 8. Create index for Kallist with Kallisto index script
- 9. Quantify transcript abundance with Kallist quant script
- 10. Download kallisto abundance files to local computer & run Sleuth in R to generate pairwise comparison lists OR MaSigPro/Mfuzz script to generate cluster lists**
- 11. Use seqtk to generate Dm_proteins file (DNA sequences -> AA sequences), then use HPCgridrunner script to BLAST
