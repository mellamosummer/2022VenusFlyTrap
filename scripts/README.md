- copy over reads to scratch directory on xfer node
  - cp -r /project/jlmlab/VFT.tar.gz /scratch/srb67793/VFT_Analysis2/

- untar and unzip files
  - tar -xf
  - qlogin
  - bzip2 -d /scratch/srb67793/VFT_Analysis2/VFT/Raw_reads/*.bz2

- QC reads before trimming with FASTQC script
- summarize FASTQC results with MultiQC
- Trim reads using trimmomatic script
- QC reads after trimming with FASTQC script
- summarize FASTQC results with MultiQC
- Create index for Kallist with Kallisto index script
- Quantify transcript abundance with Kallist quant script
- Download kallisto abundance files to local computer & run Sleuth in R to generate pairwise comparison lists OR MaSigPro/Mfuzz script to generate cluster lists**
- Use seqtk to generate Dm_proteins file (DNA sequences -> AA sequences), then use HPCgridrunner script to BLAST entire Dm proteins file
- use grep commands to pull BLAST hits from Dm proteins BLAST results for each Mfuzz cluster
