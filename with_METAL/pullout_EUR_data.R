#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_meta_analysis.R <input_file.tsv.gz>")
}

input_file <- args[1]
output_dir <- args[2]

#input_file="/mnt/project/publically_available_summary_statistics/meta/L12_PSORIASIS_meta_out.tsv.gz"

# Extract basename for output file naming
library(tools)
library(dplyr)
base_name <- file_path_sans_ext(file_path_sans_ext(basename(input_file))) # Remove .tsv.gz
pheno_code <- gsub("_meta_out", "", base_name)
output_file1 <- paste0(pheno_code, "_ukb_EUR.tsv")
output_file2 <- paste0(pheno_code, "_mvp_EUR.tsv")
output_file3 <- paste0(pheno_code, "_fg_EUR.tsv")


# Load necessary packages
library(data.table)

# Start timing
start_time <- Sys.time()

# Read the data efficiently
cat("Reading data file...\n")
gwas_data <- fread(input_file, header = TRUE)

# Load endpoints file as data.table for better performance
endpoints <- fread("/mnt/project/publically_available_summary_statistics/meta/endpoints.tsv",
                   header = TRUE, sep = "\t", fill = TRUE, quote = "", na.strings = c("NA", ""))

# Get phenotype size info (using data.table syntax)
pheno_size <- endpoints[phenocode == pheno_code]
if (nrow(pheno_size) < 1) {
  stop("Missing required phenotype sample size for phenocode: ", pheno_code)
}

# Check if required columns exist
required_cols <- c("fg_beta", "fg_sebeta", "MVP_EUR_beta", "MVP_EUR_sebeta", 
                  "ukbb_beta", "ukbb_sebeta")
                  
missing_cols <- required_cols[!required_cols %in% names(gwas_data)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse=", "))
}

# Check if AF columns exist (optional)
af_cols <- c("fg_af_alt", "MVP_EUR_af_alt", "ukbb_af_alt")
missing_af_cols <- af_cols[!af_cols %in% names(gwas_data)]
if (length(missing_af_cols) > 0) {
  cat("Warning: Some allele frequency columns are missing:", paste(missing_af_cols, collapse=", "), "\n")
  cat("Weighted allele frequency calculation may be incomplete or skipped\n")
}

# Perform meta-analysis using vectorized operations
cat("Performing meta-analysis...\n")

# Create columns to track non-missing values (vectorized)
gwas_data[, `:=`(
  fg_valid = !is.na(fg_beta) & !is.na(fg_sebeta),
  MVP_EUR_valid = !is.na(MVP_EUR_beta) & !is.na(MVP_EUR_sebeta),
  ukbb_valid = !is.na(ukbb_beta) & !is.na(ukbb_sebeta)
)]

# Add phenotype size info once instead of per-row
# Extract as numeric values to avoid data.table reference issues
fg_n_cases <- as.numeric(pheno_size$fg_n_cases[1])
fg_n_controls <- as.numeric(pheno_size$fg_n_controls[1])
MVP_n_cases <- as.numeric(pheno_size$MVP_EUR_n_cases[1])
MVP_n_controls <- as.numeric(pheno_size$MVP_EUR_n_controls[1])
ukbb_n_cases <- as.numeric(pheno_size$ukbb_n_cases[1])
ukbb_n_controls <- as.numeric(pheno_size$ukbb_n_controls[1])

gwas_data[, `:=`(
  fg_n = fg_n_cases + fg_n_controls,
  MVP_n = MVP_n_cases + MVP_n_controls,
  ukbb_n = ukbb_n_cases + ukbb_n_controls
)]

ukb <- gwas_data %>% filter(ukbb_valid) %>% select(`#CHR`,POS,REF,ALT,SNP, ukbb_af_alt, ukbb_beta, ukbb_sebeta,ukbb_pval, ukbb_n) %>% rename(CHR=`#CHR`) 
mvp <- gwas_data %>% filter(MVP_EUR_valid) %>% select(`#CHR`,POS,REF,ALT,SNP, MVP_EUR_af_alt, MVP_EUR_beta, MVP_EUR_sebeta,MVP_EUR_pval, MVP_n) %>% rename(CHR=`#CHR`) 
fg <- gwas_data %>% filter(fg_valid) %>% select(`#CHR`,POS,REF,ALT,SNP, fg_af_alt, fg_beta, fg_sebeta, fg_pval, fg_n) %>% rename(CHR=`#CHR`) 

fwrite(ukb, file = paste0(output_dir,"/",output_file1, ".gz"), sep = "\t", compress = "gzip")
fwrite(mvp, file = paste0(output_dir,"/",output_file2, ".gz"), sep = "\t", compress = "gzip")
fwrite(fg, file = paste0(output_dir,"/",output_file3, ".gz"), sep = "\t", compress = "gzip")

# Report timing and statistics
end_time <- Sys.time()
cat("Processing completed in", difftime(end_time, start_time, units = "secs"), "seconds\n")
