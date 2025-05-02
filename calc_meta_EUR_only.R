# Load necessary packages
library(data.table)

# Start timing
start_time <- Sys.time()

# Read the data
gwas_data <- fread("/mnt/project/publically_available_summary_statistics/meta/G6_MS_meta_out.tsv.gz", header = TRUE)

# Perform meta-analysis using vectorized operations
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

# Select only the requested columns for the final output
# Assuming the first column is named "#CHR" with the hash symbol
final_columns <- c("#CHR", "POS", "REF", "ALT", "SNP", "rsid", 
                   "meta_EUR_N", "meta_beta", "meta_se", "meta_p", "Q", "Q_pval")

# Create final results table with only selected columns
final_results <- filtered_gwas_data[, ..final_columns]

# Save final results to a new file
fwrite(final_results, "gwas_meta_analysis_final_results.tsv", sep = "\t")

# Report timing and statistics
end_time <- Sys.time()
cat("Processing completed in", difftime(end_time, start_time, units = "secs"), "seconds\n")
cat("Original number of rows:", nrow(gwas_data), "\n")
cat("Filtered number of rows:", nrow(final_results), "\n")
cat("Distribution of studies used per analysis:\n")
print(table(final_results$meta_EUR_N))
