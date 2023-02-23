## Methods ##

Trimmomatic v0.39 (Bolger et al., 2014) was used for quality trimming and adapter clipping. Sequencing reads were visualized for quality using FASTQC v0.11.9 (Andrews, 2010) & MultiQC v1.8 (Ewels et al., 2016) before and after trimming.

Trimmed reads were pseudo-aligned to the _D. muscipula_ reference transcriptome index then quantified using Kallisto v0.46.1 (Bray et al., 2016). A principal component analysis (PCA) was generated using the R software package Sleuth (Pimentel et al., 2017) to visualize the variance between sample groups and replicates. Kallisto transcript abundances were then used for pair-wise differential gene expression analyses between petioles and all traps at all time points, and between traps given prey and mechanically triggered at 5 min, 1 hr, & 24 hr.

To further explore genes involved in the trap mechanism, time series analysis was performed using the R software package maSigPro (Conesa et al., 2006) for prey (5 min,1 hr, 3 hr, 12 hr, 24 hr, 48 hr, 72 hr) and no prey (0 min, 5 min, 1 hr, 24 hr) treatments. A multidimensional scaling (MDS) plot was generated to identify outliers. A single library (JMHK) for a leaf trap collected at 24 hrs for the prey treatment was removed from the analysis. 

maSigPro uses a two-step regression strategy Ð first, identifying differentially expressed genes, then clustering DEG's into significant time-structured profiles. A k-means clustering approach was used to group differentially expressed genes by similar profiles. Cluster number was determined by visually inspecting profiles for a range of k values.

Differentially expressed transcripts were then translated using SeqKit (Shen et al., 2016) and protein sequences were BLASTed against a custom database of annotated _Arabidopsis thaliana_ proteins (Araport11_pep_20220914) from The Arabidopsis Information Resource (Berardini et al., 2015) and the Swiss-Prot database (Bairoch & Boeckmann, 1994). Functional enrichment analysis was conducted using the web-based application, GOrilla (Eden et al., 2009), with _A. thaliana_ protein homologs of the _D. muscipula_ transcriptome as a background set and the _A. thaliana_ protein homologs for each time point by treatment as the target set.

Scripts and results for these analyses can be found at: www.github.com/mellamosummer/2022VenusFlyTrap/tree/main/Cleaned_Workflow












## Citations ##

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Bairoch, A., & Boeckmann, B. (1994). The SWISS-PROT protein sequence data bank: Current status. Nucleic Acids Research, 22(17), 3578Ð3580.

Berardini, T. Z., Reiser, L., Li, D., Mezheritsky, Y., Muller, R., Strait, E., & Huala, E. (2015). The arabidopsis information resource: Making and mining the Ògold standardÓ annotated reference plant genome: Tair: Making and Mining the ÒGold StandardÓ Plant Genome. Genesis, 53(8), 474Ð485. https://doi.org/10.1002/dvg.22877

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114Ð2120. https://doi.org/10.1093/bioinformatics/btu170

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), 525Ð527. https://doi.org/10.1038/nbt.3519

Conesa, A., Nueda, M. J., Ferrer, A., & Tal—n, M. (2006). maSigPro: A method to identify significantly differential expression profiles in time-course microarray experiments. Bioinformatics, 22(9), 1096Ð1102. https://doi.org/10.1093/bioinformatics/btl056

Eden, E., Navon, R., Steinfeld, I., Lipson, D., & Yakhini, Z. (2009). GOrilla: A tool for discovery and visualization of enriched GO terms in ranked gene lists. BMC Bioinformatics, 10(1), 48. https://doi.org/10.1186/1471-2105-10-48

Ewels, P., Magnusson, M., Lundin, S., & KŠller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047Ð3048. https://doi.org/10.1093/bioinformatics/btw354

Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature Methods, 14(7), 687Ð690. https://doi.org/10.1038/nmeth.4324

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962
