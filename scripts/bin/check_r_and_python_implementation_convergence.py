#!/usr/bin/env python3

# Import necessary modules:
# RNApysoforms is assumed to be a package that contains the function to read GTF files.
# Polars is a fast DataFrame library similar to pandas, which will be used to handle CSV files.
import RNApysoforms as RNApy
import polars as pl

# Read the ENSEMBL GTF file using the RNApysoforms package and store it as a Polars DataFrame.
# The GTF file contains genomic feature data, and it is loaded from the relative path "../data/"
df_py = RNApy.read_ensembl_gtf("../../data/Homo_sapiens_chr21_and_Y.GRCh38.110.gtf")

# Read the CSV file produced by an R script into a Polars DataFrame.
# The CSV file is located in the current directory and named 'r_output.csv'.
df_R = pl.read_csv("./r_output.csv")

# Compare the two DataFrames using the `equals` method, which checks for equality of content.
# This method returns True if both DataFrames have identical data, else it returns False.
are_equal = df_py.equals(df_R)

# If the DataFrames are equal, print a success message indicating they match.
# Otherwise, print a message indicating a mismatch between the two DataFrames.
if are_equal:
    print("\nTHE DATAFRAMES PRODUCED BY THE R AND THE PYTHON FUNCTIONS ARE THE SAME")
else:
    print("\nTHE DATAFRAMES PRODUCED BY THE R AND THE PYTHON FUNCTION ARE DIFFERENT")

# Print the first 5 rows of the Python-generated DataFrame for inspection.
print("\nHere are the first 5 lines of the Python generated dataframe\n")
print(df_py.head())

# Print the first 5 rows of the R-generated DataFrame for inspection.
print("\nHere are the first 5 lines of the R generated dataframe\n")
print(df_R.head())

