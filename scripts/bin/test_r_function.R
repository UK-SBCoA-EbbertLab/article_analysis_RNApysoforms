#!/usr/bin/env Rscript

# Load necessary libraries from the Tidyverse suite
library(readr)    # For reading and writing data
library(dplyr)    # For data manipulation
library(stringr)  # For string operations

#' Read an ENSEMBL GTF File and Extract Relevant Genomic Features
#'
#' This function reads a Gene Transfer Format (GTF) file from the ENSEMBL database and processes it
#' to extract specific genomic features, namely 'exon' and 'CDS' (Coding DNA Sequence). It parses
#' the 'attributes' column to retrieve key information such as gene and transcript identifiers,
#' names, biotypes, and exon numbers. The function performs validation checks on the input file
#' to ensure correctness and consistency with the expected GTF format.
#'
#' @param path A string specifying the file path to the ENSEMBL GTF file. The file must have a '.gtf' extension.
#'
#' @return A tibble (data frame) containing the following columns:
#' \describe{
#'   \item{gene_id}{Identifier for the gene.}
#'   \item{gene_name}{Name of the gene. If missing, filled with `gene_id`.}
#'   \item{transcript_id}{Identifier for the transcript.}
#'   \item{transcript_name}{Name of the transcript. If missing, filled with `transcript_id`.}
#'   \item{transcript_biotype}{Biotype classification of the transcript.}
#'   \item{seqnames}{Chromosome or sequence name.}
#'   \item{strand}{Strand information ('+' or '-').}
#'   \item{type}{Feature type ('exon' or 'CDS').}
#'   \item{start}{Start position of the feature.}
#'   \item{end}{End position of the feature.}
#'   \item{exon_number}{Exon number within the transcript.}
#' }
#'
#' @details
#' - **File Validation:**
#'   - Checks if the specified file exists.
#'   - Ensures the path points to a file and not a directory.
#'   - Verifies that the file has a '.gtf' extension.
#'
#' - **Data Reading:**
#'   - Reads the GTF file using `read_tsv` from the `readr` package.
#'   - Assigns appropriate column names and data types.
#'   - Skips lines starting with '#' as they are comments.
#'
#' - **Data Processing:**
#'   - Filters the data to include only 'exon' and 'CDS' feature types.
#'   - Extracts relevant attributes from the 'attributes' column using regular expressions.
#'   - Fills missing `gene_name` with `gene_id` and `transcript_name` with `transcript_id`.
#'   - Selects and orders the necessary columns for the final output.
#'
#' - **Validation:**
#'   - Checks for any missing values in the resulting data frame to ensure consistency.
#'   - Converts the `exon_number` column to integer type.
#'
#' @examples
#' \dontrun{
#' # Read a GTF file and write the processed data to a CSV file
#' df <- read_ensembl_gtf("../data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")
#' write_csv(df, "r_output.csv")
#' }
#'
#' @export
read_ensembl_gtf <- function(path) {
  
  # ------------------------------
  # 1. File Validation
  # ------------------------------
  
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
  
  # ------------------------------
  # 2. Define Column Names and Types
  # ------------------------------
  
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
  
  # ------------------------------
  # 3. Read the GTF File
  # ------------------------------
  
  # Use read_tsv to read the GTF file
  # - `col_names` assigns the predefined column names
  # - `comment` skips lines starting with '#'
  # - `col_types` enforces the specified data types
  df <- read_tsv(
    file = path,
    col_names = column_names,
    comment = "#",
    col_types = col_types
  )
  
  # ------------------------------
  # 4. Filter for Relevant Features
  # ------------------------------
  
  # Retain only rows where the 'type' is either 'exon' or 'CDS'
  df_filtered <- df %>%
    filter(type %in% c("exon", "CDS"))
  
  # ------------------------------
  # 5. Extract Attributes from 'attributes' Column
  # ------------------------------
  
  # Use regular expressions to extract key-value pairs from the 'attributes' column
  df_extracted <- df_filtered %>%
    mutate(
      # Extract 'gene_id' using regex and capture the value inside quotes
      gene_id = str_match(attributes, 'gene_id "([^"]+)"')[,2],
      
      # Extract 'gene_name'; may be missing, hence can result in NA
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
  
  # ------------------------------
  # 6. Handle Missing Values
  # ------------------------------
  
  # Replace missing 'gene_name' with 'gene_id' and 'transcript_name' with 'transcript_id'
  df_filled <- df_extracted %>%
    mutate(
      gene_name = coalesce(gene_name, gene_id),
      transcript_name = coalesce(transcript_name, transcript_id)
    )
  
  # ------------------------------
  # 7. Select and Reorder Relevant Columns
  # ------------------------------
  
  # Choose only the necessary columns for the final output and arrange them in order
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
  
  # ------------------------------
  # 8. Validate the Final Data Frame
  # ------------------------------
  
  # Check for any missing (NA) values in the resulting data frame
  if (anyNA(df_result)) {
    stop("This GTF file is not consistent with the 2024 ENSEMBL GTF format. See the vignette for handling other GTF formats.")
  }
  
  # ------------------------------
  # 9. Cast 'exon_number' to Integer
  # ------------------------------
  
  # Convert the 'exon_number' column to integer type
  # This ensures that exon numbers are stored as integers for downstream analysis
  df_result <- df_result %>%
    mutate(exon_number = as.integer(exon_number))
  
  # ------------------------------
  # 10. Return the Processed Data Frame
  # ------------------------------
  
  return(df_result)
}

# ------------------------------
# Main Execution
# ------------------------------

# Read the ENSEMBL GTF file using the defined function
df <- read_ensembl_gtf("../../data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

# Write the processed data frame to a CSV file named 'r_output.csv' in the current directory
write_csv(df, "r_output.csv")

