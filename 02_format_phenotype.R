# Read the full file as plain text
lines <- readLines("series.txt")

# Extract line 27 (GEO accession IDs)
geo_line <- lines[27]
geo_values <- unlist(strsplit(geo_line, "\t"))[-1]  # Remove the initial label
geo_values <- gsub('"', '', geo_values)  # Remove quotes

# Extract line 40 (Phenotype info)
pheno_line <- lines[40]
pheno_values <- unlist(strsplit(pheno_line, "\t"))[-1]
pheno_values <- gsub('"', '', pheno_values)  # Remove quotes
pheno_values <- gsub("phenotype: ", "", pheno_values)

# Convert phenotype to 0 (Ctrl) or 1 (MND)
trait_numeric <- ifelse(pheno_values == "Ctrl", 0, 1)

# Create TASSEL-compatible phenotype data frame
phenodata <- data.frame(
  taxa = geo_values,
  `<Trait>` = trait_numeric,
  stringsAsFactors = FALSE
)

# Save to tab-delimited file
write.table(phenodata, file = "phenotype_for_tassel.txt", sep = "\t", row.names = FALSE, quote = FALSE)
View(phenodata)