#!/usr/bin/env bash


## Pull singularity container and move it to proper folder
singularity pull --arch amd64 library://ebbertlab/rnapysoforms/article_analysis_rnapysoforms:2024-10-23
mv article_analysis_rnapysoforms_2024-10-23.sif ../../singularity_container/ 
