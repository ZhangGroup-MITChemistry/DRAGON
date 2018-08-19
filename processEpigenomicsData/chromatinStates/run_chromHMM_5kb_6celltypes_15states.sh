#!/bin/bash

# binarize the data
java -jar /path_to_ChromHMM/ChromHMM.jar BinarizeBed -b 5000 /path_to_ChromHMM/CHROMSIZES/hg19.txt /path_to_bedfiles ./cellmarktable.txt ./binarizedData_5kb_6celltype_15states

# build states
java -jar /path_to_ChromHMM/ChromHMM.jar LearnModel -b 5000 -r 1000 ./binarizedData_ucsc_5kb_6celltype_15states ./OUTPUTSAMPLE_5kb_6celltype_15states 15 hg19
