#!/usr/bin/env bash



## Create .csv file using the R implementation of read_ensembl_gtf()
singularity exec ../../singularity_container/article_analysis_rnapysoforms_2024-10-23.sif ../bin/test_r_function.R > /dev/null 2>&1


## Test if the R and Python implementations produce the same results
singularity exec ../../singularity_container/article_analysis_rnapysoforms_2024-10-23.sif ../bin/check_r_and_python_implementation_convergence.py


## Delete R output file
rm -f ./r_output.csv
