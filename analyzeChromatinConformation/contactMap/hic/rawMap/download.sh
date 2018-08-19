#!/bin/bash

# Download the raw HiC data from the following database
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525

# Change the link if a different cell type is selected 

wget -O GM12878.tar.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz
tar -xavf GM12878.tar.gz
