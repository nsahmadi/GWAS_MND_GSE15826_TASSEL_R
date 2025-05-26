# Load required libraries
library(qqman)
library(dplyr)

# Step 1: Import TASSEL result file
# Replace "MLM_tassel.txt" with the correct path to your TASSEL output file
MLM_tassel <- read.table("MLM_tassel.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 2: Clean and prepare data
MLM_clean <- MLM_tassel %>%
  filter(Marker != "None") %>%                        # Remove header-like rows
  mutate(
    Chr = as.numeric(Chr),                            # Convert chromosome to numeric
    Pos = as.numeric(Pos),                            # Convert position to numeric
    p = as.numeric(p)                                 # Convert p-values to numeric
  ) %>%
  filter(!is.na(Chr) & !is.na(Pos) & !is.na(p))       # Remove rows with NA in key columns

# Step 3: Calculate Bonferroni correction threshold
bonferroni_threshold <- -log10(0.05 / nrow(MLM_clean))

# Step 4: Plot Manhattan and QQ plots side by side
par(mfrow = c(1, 2))  # Set layout: 1 row, 2 columns

# Manhattan plot
manhattan(
  MLM_clean,
  chr = "Chr",
  bp = "Pos",
  snp = "Marker",
  p = "p",
  col = c("steelblue", "darkorange"),
  genomewideline = bonferroni_threshold,
  suggestiveline = FALSE,
  main = "Manhattan Plot"
)

# QQ plot
qq(MLM_clean$p, main = "QQ Plot")