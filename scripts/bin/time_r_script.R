#!/usr/bin/env Rscript

# Load necessary libraries
library(readr)         # For reading and writing data
library(dplyr)         # For data manipulation
library(stringr)       # For string operations
library(purrr)         # For iterating over lists
library(ggtranscript)  # For creating transcript plots
library(microbenchmark) # For timing execution
library(ggplot2)       # For additional plotting capabilities


# Step 1: Define a function to read and process the ENSEMBL GTF file
read_ensembl_gtf <- function(path) {

    # Check if the file exists at the specified path
    if (!file.exists(path)) {
        stop(sprintf("File '%s' does not exist. Please provide a valid file path.", path))
    }

    # Check if the path points to a file and not a directory
    if (file.info(path)$isdir) {
        stop(sprintf("'%s' is not a file. Please provide a valid file path.", path))
    }

    # Ensure the file has a '.gtf' extension (case-insensitive)
    if (tolower(tools::file_ext(path)) != "gtf") {
        stop("File must have a '.gtf' extension.")
    }

    # Specify the expected column names for the GTF file
    column_names <- c(
        "seqnames",    # Chromosome or sequence name
        "source",      # Annotation source
        "type",        # Feature type (e.g., exon, CDS)
        "start",       # Start position of the feature
        "end",         # End position of the feature
        "score",       # Score value (usually '.')
        "strand",      # Strand information ('+' or '-')
        "phase",       # Reading frame phase
        "attributes"   # Additional attributes in key-value pairs
    )

    # Define the data types for each column using readr's col_types
    col_types <- cols(
        seqnames = col_character(),
        source = col_character(),
        type = col_character(),
        start = col_integer(),
        end = col_integer(),
        score = col_character(),
        strand = col_character(),
        phase = col_character(),
        attributes = col_character()
    )

    # Use read_tsv to read the GTF file into a data frame
    df <- read_tsv(
        file = path,
        col_names = column_names,
        comment = "#",   # Skip comment lines starting with '#'
        col_types = col_types
    )

    # Retain only rows where the 'type' is either 'exon' or 'CDS'
    df_filtered <- df %>%
        filter(type %in% c("exon", "CDS"))

    # Extract key-value pairs from the 'attributes' column using regular expressions
    df_extracted <- df_filtered %>%
        mutate(
            # Extract 'gene_id' from attributes
            gene_id = str_match(attributes, 'gene_id "([^"]+)"')[,2],

            # Extract 'gene_name'; may be missing and result in NA
            gene_name = str_match(attributes, 'gene_name "([^"]+)"')[,2],

            # Extract 'transcript_id'
            transcript_id = str_match(attributes, 'transcript_id "([^"]+)"')[,2],

            # Extract 'transcript_name'; may be missing
            transcript_name = str_match(attributes, 'transcript_name "([^"]+)"')[,2],

            # Extract 'transcript_biotype'
            transcript_biotype = str_match(attributes, 'transcript_biotype "([^"]+)"')[,2],

            # Extract 'exon_number'; may be missing or non-integer
            exon_number = str_match(attributes, 'exon_number "([^"]+)"')[,2]
        )

    # Replace missing 'gene_name' and 'transcript_name' with their IDs
    df_filled <- df_extracted %>%
        mutate(
            gene_name = coalesce(gene_name, gene_id),
            transcript_name = coalesce(transcript_name, transcript_id)
        )

    # Select and arrange necessary columns for the final output
    df_result <- df_filled %>%
        select(
            gene_id,
            gene_name,
            transcript_id,
            transcript_name,
            transcript_biotype,
            seqnames,
            strand,
            type,
            start,
            end,
            exon_number
        )

    # Check for any missing values in the data frame
    if (anyNA(df_result)) {
        stop("This GTF file is not consistent with the 2024 ENSEMBL GTF format. See the vignette for handling other GTF formats.")
    }

    # Convert 'exon_number' to integer type for downstream analysis
    df_result <- df_result %>%
        mutate(exon_number = as.integer(exon_number))

    return(df_result)
}



# Step 2: Read the GTF file and create a data frame with genomic annotations
df <- read_ensembl_gtf("../../data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

# Step 3: Extract a list of unique gene IDs while maintaining the original order
gene_ids_list <- df %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    pull(gene_id)

# Step 4: Create a list of data frames, each filtered by a unique gene_id
df_list <- gene_ids_list %>%
    map(~ df %>% filter(gene_id == .x))

# Step 5: Define a function to loop through gene data frames and generate plots, render and save
loop_through_gene_dfs <- function(df_list) {

    # Initialize a counter for naming output plot files
    j <- 1

    # Loop through each gene data frame in the list
    for (gene_df in df_list) {

        # Open a PDF device to save the plot
        pdf(paste("../../figures/ggtranscript/plot_", j, ".pdf", sep=""))

        # Generate the plot for the gene's transcripts and print it to the PDF
        print(
            gene_df %>%
                ggplot(aes(
                    xstart = start,          # Start position for exons
                    xend = end,              # End position for exons
                    y = transcript_id        # Grouping by transcript ID on the y-axis
                )) +
                geom_range(
                    aes(fill = transcript_biotype)  # Color coding by transcript biotype
                ) +
                geom_intron(
                    data = to_intron(gene_df, "transcript_id"),  # Generate intron data
                    aes(strand = strand)                         # Strand information for introns
                )
        )

        # Close the PDF device
        dev.off()

        # Increment the counter for the next plot
        j <- j + 1

    }
}

# Step 7: Iterate 100 times, timing each iteration, and collect timing results
timing_df <- data.frame(iteration = integer(), elapsed_time_r = numeric())

# Outer loop using system.time to measure each iteration's timing
for (i in 1:100) {
    
    # Time each iteration of the plotting process with saving
    time_taken <- system.time({loop_through_gene_dfs(df_list)})

    # Append the iteration number and elapsed time to the timing data frame
    timing_df <- rbind(timing_df, data.frame(iteration = i, elapsed_time = time_taken["elapsed"]))
}

# Step 8: Save the timing results to a CSV file for analysis
write.csv(timing_df, "../../data/r_timing_results.csv", row.names = FALSE)
