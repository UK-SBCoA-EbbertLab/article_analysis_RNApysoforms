#!/usr/bin/env python3

import time
import polars as pl
import RNApysoforms as RNApy

"""
Script Purpose:
This script processes genomic annotation data from an ENSEMBL GTF file to generate plots of transcript structures for each gene. It repeats the plotting process 100 times to measure the execution time of the entire operation, saving the timing results to a CSV file.

Workflow:
1. Read the GTF file and create a DataFrame containing genomic annotations.
2. Extract a list of unique gene IDs from the DataFrame.
3. Create a list of DataFrames, each containing data for a single gene.
4. For 100 iterations:
   a. Start timing the iteration.
   b. For each gene DataFrame:
       i. Generate intron annotations.
       ii. Create plotting traces for the gene's transcripts.
       iii. Generate and save a plot for the gene.
   c. Record the elapsed time for the iteration.
5. Save the timing results to a CSV file.
"""

# Step 1: Read the ENSEMBL GTF file and create a DataFrame with genomic annotations
df = RNApy.read_ensembl_gtf("../../data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

# Step 2: Extract a list of unique gene IDs from the DataFrame, maintaining the original order
gene_ids_list = df["gene_id"].unique(maintain_order=True).to_list()

# Step 3: Create a list of DataFrames, each filtered by a unique gene_id
df_list = []
for gene_id in gene_ids_list:
    # Filter the main DataFrame to include only records for the current gene_id
    df_gene = df.filter(pl.col("gene_id") == gene_id)
    # Add the gene-specific DataFrame to the list
    df_list.append(df_gene)

# Initialize a list to store timing results for each iteration
timings = []

# Step 4: Repeat the plotting process 100 times, timing each iteration
for i in range(1, 101):

    # a. Record the start time of the iteration
    start_time = time.time()

    # Initialize a counter for naming output plot files
    j = 1

    # b. Process each gene DataFrame in the list
    for gene_df in df_list:

        # i. Generate intron annotations from exon data for the gene
        gene_df = RNApy.to_intron(gene_df)

        # ii. Create plotting traces for the gene's transcripts
        traces = RNApy.make_traces(
            annotation=gene_df,                # DataFrame with gene annotations
            x_start="start",                   # Column name for the start position
            x_end="end",                       # Column name for the end position
            y='transcript_id',                 # Column name for the y-axis grouping
            annotation_hue="transcript_biotype",  # Column name for color coding
        )

        # iii. Generate the plot using the traces, adding a subplot title
        fig = RNApy.make_plot(
            traces=traces, 
            subplot_titles=["Transcript Structure"]
        )

        # Save the plot to an HTML file with a unique filename
        fig.write_html(("../../figures/RNApysoforms/plot_" + str(j) + ".html"))

        # Increment the counter for the next plot
        j += 1

    # c. Record the end time of the iteration and calculate the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time  # Note: 'end_time_2' is undefined; intended to calculate elapsed time

    # Add the iteration number and elapsed time to the timings list
    timings.append({'iteration': i, 'elapsed_time_python': elapsed_time})

# Step 5: Create a DataFrame from the timing results
timing_df = pl.DataFrame(timings)

# Save the timing results to a CSV file for analysis
timing_df.write_csv("../../data/python_timing_results.csv")

