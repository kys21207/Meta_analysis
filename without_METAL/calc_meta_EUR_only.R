#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_meta_analysis.R <input_file.tsv.gz>")
}

input_file <- args[1]

# Extract basename for output file naming
library(tools)
base_name <- file_path_sans_ext(file_path_sans_ext(basename(input_file))) # Remove .tsv.gz
output_file <- paste0(base_name, "_EUR.tsv")

cat("Processing", input_file, "-> Output:", output_file, "\n")

# Load necessary packages
library(data.table)

# Start timing
start_time <- Sys.time()

# Read the data efficiently
cat("Reading data file...\n")
gwas_data <- fread(input_file, header = TRUE)

# Check if required columns exist
required_cols <- c("fg_beta", "fg_sebeta", "MVP_EUR_beta", "MVP_EUR_sebeta", 
                  "ukbb_beta", "ukbb_sebeta")
                  
missing_cols <- required_cols[!required_cols %in% names(gwas_data)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse=", "))
}

# Perform meta-analysis using vectorized operations
cat("Performing meta-analysis...\n")

# Create columns to track non-missing values
gwas_data[, `:=`(
  fg_valid = !is.na(fg_beta) & !is.na(fg_sebeta),
  MVP_EUR_valid = !is.na(MVP_EUR_beta) & !is.na(MVP_EUR_sebeta),
  ukbb_valid = !is.na(ukbb_beta) & !is.na(ukbb_sebeta)
)]

# Calculate the number of valid studies for each row
gwas_data[, meta_EUR_N := fg_valid + MVP_EUR_valid + ukbb_valid]

# Calculate weights using vectorized operations
gwas_data[, `:=`(
  fg_weight = ifelse(fg_valid, 1 / (fg_sebeta^2), 0),
  MVP_EUR_weight = ifelse(MVP_EUR_valid, 1 / (MVP_EUR_sebeta^2), 0),
  ukbb_weight = ifelse(ukbb_valid, 1 / (ukbb_sebeta^2), 0)
)]

# Calculate total weight
gwas_data[, total_weight := fg_weight + MVP_EUR_weight + ukbb_weight]

# Calculate meta-analysis beta (weighted average) - FIXED
gwas_data[, meta_beta := {
  sum_weighted_beta <- 0
  if (fg_valid) sum_weighted_beta <- sum_weighted_beta + fg_weight * fg_beta
  if (MVP_EUR_valid) sum_weighted_beta <- sum_weighted_beta + MVP_EUR_weight * MVP_EUR_beta
  if (ukbb_valid) sum_weighted_beta <- sum_weighted_beta + ukbb_weight * ukbb_beta
  sum_weighted_beta / total_weight
}, by = seq_len(nrow(gwas_data))]

# Calculate meta-analysis standard error
gwas_data[, meta_se := sqrt(1 / total_weight)]

# Calculate meta p-value
gwas_data[, `:=`(
  meta_z = meta_beta / meta_se,
  meta_p = 2 * pnorm(-abs(meta_beta / meta_se))
)]

# Cochran's Q test - vectorized
gwas_data[, Q := {
  q_sum <- 0
  if(fg_valid) q_sum <- q_sum + fg_weight * (fg_beta - meta_beta)^2
  if(MVP_EUR_valid) q_sum <- q_sum + MVP_EUR_weight * (MVP_EUR_beta - meta_beta)^2
  if(ukbb_valid) q_sum <- q_sum + ukbb_weight * (ukbb_beta - meta_beta)^2
  q_sum
}, by = seq_len(nrow(gwas_data))]

# Calculate Q p-value, accounting for df
gwas_data[, Q_pval := pchisq(Q, df = meta_EUR_N - 1, lower.tail = FALSE)]

# Set meta-analysis results to NA where there are fewer than 2 studies
gwas_data[meta_EUR_N < 2, `:=`(
  meta_beta = NA_real_,
  meta_se = NA_real_,
  meta_p = NA_real_,
  Q = NA_real_,
  Q_pval = NA_real_
)]

# Filter rows with meta_EUR_N > 2 and non-missing meta_beta and meta_se
cat("Filtering results...\n")
filtered_gwas_data <- gwas_data[meta_EUR_N > 2 & !is.na(meta_beta) & !is.na(meta_se)]

# Select only the requested columns for the final output
# Check if these columns exist in the data
final_columns <- c("#CHR", "POS", "REF", "ALT", "SNP", "rsid", 
                   "meta_EUR_N", "meta_beta", "meta_se", "meta_p", "Q", "Q_pval")
                   
present_columns <- final_columns[final_columns %in% names(filtered_gwas_data)]
missing_cols <- final_columns[!final_columns %in% names(filtered_gwas_data)]

if (length(missing_cols) > 0) {
  cat("Warning: Some requested columns are missing:", paste(missing_cols, collapse=", "), "\n")
  cat("Using available columns:", paste(present_columns, collapse=", "), "\n")
}

# Create final results table with only selected columns
final_results <- filtered_gwas_data[, ..present_columns]

# Remove intermediate calculation columns
gwas_data[, c("fg_valid", "MVP_EUR_valid", "ukbb_valid", "fg_weight", 
              "MVP_EUR_weight", "ukbb_weight", "total_weight") := NULL]

# Save final results to a new file
cat("Writing results to", output_file, "\n")
fwrite(final_results, output_file, sep = "\t")

# Report timing and statistics
end_time <- Sys.time()
cat("Processing completed in", difftime(end_time, start_time, units = "secs"), "seconds\n")
cat("Original number of rows:", nrow(gwas_data), "\n")
cat("Filtered number of rows:", nrow(final_results), "\n")
cat("Distribution of studies used per analysis:\n")
print(table(final_results$meta_EUR_N))
