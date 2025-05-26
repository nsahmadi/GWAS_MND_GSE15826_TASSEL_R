# Suppose your MLM data is in a data frame called MLM_tassel

# 1. Extract the p-values (overall SNP test)
p_values <- MLM_tassel$p

# 2. Apply Bonferroni correction
# Number of tests = number of non-NA p-values
num_tests <- sum(!is.na(p_values))
bonferroni_p <- p.adjust(p_values, method = "bonferroni", n = num_tests)

# 3. Add corrected p-values back to your data frame
MLM_tassel$bonferroni_p <- bonferroni_p

# 4. Filter SNPs with significant corrected p-values (e.g., alpha=0.05)
significant_snps <- MLM_tassel[!is.na(MLM_tassel$bonferroni_p) & MLM_tassel$bonferroni_p < 0.05, ]

# 5. View significant SNPs
print(significant_snps)

# Using the same p-values vector
fdr_p <- p.adjust(MLM_tassel$p, method = "fdr", n = num_tests)

# Add to your data frame
MLM_tassel$fdr_p <- fdr_p

# Filter significant SNPs by FDR threshold (e.g., 0.05)
significant_snps_fdr <- MLM_tassel[!is.na(MLM_tassel$fdr_p) & MLM_tassel$fdr_p < 0.05, ]

print(significant_snps_fdr)