#!/bin/bash
#SBATCH --time=7-00:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=timing_script_submission_%j    # Job name
#SBATCH --ntasks 1 ## One task
#SBATCH --cpus-per-task 1
#SBATCH --nodes 1 ## Execute everything in only one node
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal         # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e timing_script_submission_%j.err             # Error file for this job.
#SBATCH -o timing_script_submission_%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr     # Project allocation account name (REQUIRED)

## Submit Python timing script
singularity exec ../../singularity_container/article_analysis_rnapysoforms_2024-10-23.sif ../bin/time_python_script.py

## Submit R timing script
singularity exec ../../singularity_container/article_analysis_rnapysoforms_2024-10-23.sif ../bin/time_r_script.R

