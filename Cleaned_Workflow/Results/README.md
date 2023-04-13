# 1. MultiQC

Quality control checks on raw sequence data, aggreagrated for all libraries.

# 2. Kallisto Quant

Abundances of transcripts quantified for each library.

# 3. Sleuth

Differential gene expression analysis.
Includes prey time series & pair-wise time point DE analysis.

# 4. maSigPro

Time series DE analysis.
Inlcudes mechanical trigger time series, prey time series, & share time point between treatments. 

# 5. BLAST

Blast results of all D muscipula proteins using Tair & SwissProt databases.
Also includes annotated gene lists from maSigPro & Sleuth analyses.

# 6. OrthoFinder

Orthogroup analysis including the following species:
Aquilegia coerulea 
Amaranthus hypochondriacus
Beta vulgaris
Portulaca amilis
Arabidopsis thaliana
Glycine max
Vitis vinifera
Solanum lycopersicum
Lactuca sativa
Aldrovanda vesiculosa
Dionaea muscipula
Drosera spatulata

# 7. Gene Co-expression analysis

Gene co-expression analysis that identifies modules of coexpressed genes. Analysis conducted with all trap samples (prey & no prey), as well a prey time series alone. BLAST annotations & orthogroups appended to result files. Enzymes from literature review are extracted from the modules and visualized as line plots. 

