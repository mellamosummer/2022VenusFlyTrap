There are 48 folders here corresponding to each sample & replicate. Each folder contains 3 files: 1) abundance.h5, 2) abundance.tsv, and 3) abundance.tsv file. These are used for the Sleuth & MaSigPro DEG analyses in R.

1) "abundance.h5 is a HDF5 binary file containing run info, abundance esimates, bootstrap estimates, and transcript length information length. This file can be read in by sleuth."
2) "abundance.tsv is a plaintext file of the abundance estimates. It does not contains bootstrap estimates. Please use the --plaintext mode to output plaintext abundance estimates. Alternatively, kallisto h5dump can be used to output an HDF5 file to plaintext. The first line contains a header for each column, including estimated counts, TPM, effective length."
3) "abundance.tsv is a json file containing information about the run"

See kallisto manual here for more information: https://pachterlab.github.io/kallisto/manual
