# Set working directory
setwd("C:/Users/Negar/Desktop/GWAS/GSE15826-motor neuron disease-164 samples")

# Read files
platform <- read.delim("platform.txt", skip = 10, header = TRUE, stringsAsFactors = FALSE)
series <- read.delim("series.txt", skip = 62, header = TRUE, stringsAsFactors = FALSE)

# Subset to first 10,000 SNPs directly
series <- series[1:10000, ]

# Match IDs
series$ID_REF <- gsub("AFFX-SNP-", "SNP_A-", series$ID_REF)
matched_platform <- platform[match(series$ID_REF, platform$ID), ]

# Metadata columns with position as integer, replace "---" or NA by 0
result <- data.frame(
  rs = matched_platform$SNP_ID,
  alleles = paste0(matched_platform$Allele.A, "/", matched_platform$Allele.B),
  chrom = matched_platform$Chromosome,
  pos = ifelse(is.na(matched_platform$Physical.Position) | matched_platform$Physical.Position == "---",
               0,
               as.integer(matched_platform$Physical.Position)),
  strand = matched_platform$STRAND,
  assembly = NA,
  center = NA,
  protLSID = NA,
  assayLSID = NA,
  panelLSID = NA,
  Qcode = NA,
  stringsAsFactors = FALSE
)

# Prepare genotypes matrix
genotypes <- as.matrix(series[, -1])
allele_a <- matched_platform$Allele.A
allele_b <- matched_platform$Allele.B

# Genotype patterns
AA_matrix <- matrix(paste0(allele_a, allele_a), nrow = nrow(genotypes), ncol = ncol(genotypes))
BB_matrix <- matrix(paste0(allele_b, allele_b), nrow = nrow(genotypes), ncol = ncol(genotypes))
AB_matrix <- matrix(paste0(allele_a, allele_b), nrow = nrow(genotypes), ncol = ncol(genotypes))

# Replace genotypes
geno_final <- genotypes
geno_final[genotypes == "AA"] <- AA_matrix[genotypes == "AA"]
geno_final[genotypes == "BB"] <- BB_matrix[genotypes == "BB"]
geno_final[genotypes %in% c("AB", "BA")] <- AB_matrix[genotypes %in% c("AB", "BA")]
geno_final[genotypes == "--" | genotypes == ""] <- NA

# Combine metadata + genotypes
final_df <- cbind(result, as.data.frame(geno_final, stringsAsFactors = FALSE))
colnames(final_df)[1] <- "rs#"

# Remove rows with missing or invalid chrom or pos (pos = 0 means invalid)
final_df <- final_df[!is.na(final_df$chrom) & !is.na(final_df$pos) & final_df$chrom != "---" & final_df$pos != 0, ]

# Sort by chromosome numerically and position
suppressWarnings(final_df$chrom_numeric <- as.numeric(final_df$chrom))
final_df <- final_df[order(final_df$chrom_numeric, final_df$pos), ]
final_df$chrom_numeric <- NULL

# Write to file
write.table(
  final_df,
  file = "tassel_genotype_input_10000snps.hmp.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Check dimensions
dim(final_df)