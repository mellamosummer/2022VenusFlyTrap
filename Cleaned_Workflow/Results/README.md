# 1. MultiQC

Quality control checks on raw sequence data, aggreagrated for all libraries.

# 2. Kallisto Quant

Abundances of transcripts quantified for each library.

# 3. Sleuth

Differential gene expression analysis in R. 
Analyses included for:
- Petioles vs Traps
- Prey time series
- Prey vs No Prey 1 hr
- Prey vs No Prey 24 hr

# 4. maSigPro

Time series DE analysis in R.
Analyses included for:
- prey time series
- no prey time series
- only shared time point between treatments (5 min, 1 hr, 24 hr)

# 5. BLAST

Blast results of all Dionaea muscipula proteins using Tair & SwissProt databases.
Includes annotated gene lists from maSigPro & Sleuth analyses.

# 6. OrthoFinder

Orthogroup analysis including the following species:
- Aquilegia coerulea 
- Amaranthus hypochondriacus
- Beta vulgaris
- Portulaca amilis
- Arabidopsis thaliana
- Glycine max
- Vitis vinifera
- Solanum lycopersicum
- Lactuca sativa
- Aldrovanda vesiculosa
- Dionaea muscipula
- Drosera spatulata

# 7. Gene Co-expression analysis

Gene co-expression analysis that identifies modules of coexpressed genes in R.
Analysis included for:
- all trap samples (prey & no prey)
- prey time series alone

BLAST annotations & orthogroups appended to result files. 
Enzymes from literature review are extracted from the modules and visualized as line plots. 

