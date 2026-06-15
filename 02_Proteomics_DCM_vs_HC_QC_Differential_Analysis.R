############################################################
# Project: Serum Proteomics Analysis in DCM vs HC
# Script: 02_Proteomics_DCM_vs_HC_QC_Differential_Analysis.R
#
# Workflow:
#   1. Data import
#   2. Sample QC
#   3. QC correlation analysis
#   4. Protein filtering
#   5. Missing value imputation
#   6. CV-based filtering
#   7. Median normalization
#   8. log2 transformation
#   9. PCA & PLS-DA
#  10. Differential expression analysis
#  11. Volcano plot
############################################################

# ==================================================
# Step 0. Load required packages
# ==================================================

library(impute)
library(preprocessCore)
library(readxl)
library(dplyr)
library(tibble)
library(stringr)
library(writexl)
library(corrplot)
library(ggplot2)
library(factoextra)

# ==================================================
# Step 1. Load data
# ==================================================

data_raw <- read_excel(
  "data/serum_DCM_Report_processed.xlsx"
)

data_raw <- as.data.frame(data_raw)

cat("Raw data dimensions: proteins =",
    nrow(data_raw),
    ", samples =",
    ncol(data_raw) - 2, "\n")

protein_names <- data_raw[[1]]
descriptions  <- data_raw[[2]]
num_data <- data_raw[, 3:ncol(data_raw)]

num_data[num_data == "NaN"] <- NA
num_data <- apply(num_data, 2, as.numeric)
num_data <- as.data.frame(num_data)

rownames(num_data) <- protein_names

# ==================================================
# Step 2. Define groups
# ==================================================

dcm_cols <- grep("^DCM", colnames(num_data))
hc_cols  <- grep("^HC", colnames(num_data))
qc_cols  <- grep("^QC", colnames(num_data))

# ==================================================
# Step 3. Sample QC (protein detection rate)
# ==================================================

protein_counts <- colSums(!is.na(num_data))
threshold <- mean(protein_counts) - 3 * sd(protein_counts)

keep_samples <- names(protein_counts[protein_counts >= threshold])

num_data <- num_data[, keep_samples]

cat("After sample QC: proteins =",
    nrow(num_data),
    ", samples =",
    ncol(num_data), "\n")

# Update group indices
dcm_cols <- grep("^DCM", colnames(num_data))
hc_cols  <- grep("^HC", colnames(num_data))
qc_cols  <- grep("^QC", colnames(num_data))

# ==================================================
# Step 4. QC correlation analysis
# ==================================================

qc_data <- num_data[, qc_cols]

corr_qc <- cor(qc_data, use = "pairwise.complete.obs")

png("results/QC_correlation_plot.png",
    width = 1000, height = 1000)

corrplot(
  corr_qc,
  method = "pie",
  type = "upper",
  tl.col = "black",
  tl.cex = 1.2
)

dev.off()

# ==================================================
# Step 5. Protein-level missing value filtering
# ==================================================

keep_proteins <- (
  apply(num_data[, dcm_cols], 1, function(x) mean(is.na(x)) <= 0.5) |
    apply(num_data[, hc_cols], 1, function(x) mean(is.na(x)) <= 0.5)
)

num_data <- num_data[keep_proteins, ]

protein_names <- protein_names[keep_proteins]
descriptions  <- descriptions[keep_proteins]

cat("After protein filtering:",
    nrow(num_data), "\n")

# ==================================================
# Step 6. Missing value imputation
# ==================================================

impute_block <- function(df_block) {
  
  na_ratio <- rowMeans(is.na(df_block))
  
  to_knn <- which(na_ratio < 0.5)
  to_min <- which(na_ratio >= 0.5)
  
  if (length(to_knn) > 0) {
    df_block[to_knn, ] <- impute.knn(
      as.matrix(df_block[to_knn, ])
    )$data
  }
  
  if (length(to_min) > 0) {
    min_val <- min(df_block[to_min, ], na.rm = TRUE)
    df_block[to_min, ][is.na(df_block[to_min, ])] <- min_val / 2
  }
  
  return(df_block)
}

num_data[, dcm_cols] <- impute_block(num_data[, dcm_cols])
num_data[, hc_cols]  <- impute_block(num_data[, hc_cols])
num_data[, qc_cols]  <- impute_block(num_data[, qc_cols])

# ==================================================
# Step 7. QC-based CV filtering
# ==================================================

qc_data <- num_data[, qc_cols]

qc_cv <- apply(qc_data, 1, function(x) sd(x) / mean(x))

keep_proteins_qc <- qc_cv < 0.5

num_data <- num_data[keep_proteins_qc, ]

protein_names <- protein_names[keep_proteins_qc]
descriptions  <- descriptions[keep_proteins_qc]

cat("After CV filtering:", nrow(num_data), "\n")

# ==================================================
# Step 8. Median normalization
# ==================================================

medians <- apply(num_data, 2, median, na.rm = TRUE)
medians[medians == 0] <- 1

num_data <- t(t(num_data) / medians)

# ==================================================
# Step 9. log2 transformation
# ==================================================

num_data <- log2(num_data + 1e-6)

# ==================================================
# Step 10. PCA (DCM vs HC)
# ==================================================

pca_data <- num_data[, c(dcm_cols, hc_cols)]

var_proteins <- apply(pca_data, 1, var, na.rm = TRUE)
pca_data <- pca_data[var_proteins > 0, ]

pca_res <- prcomp(t(pca_data), scale. = TRUE)

group <- c(
  rep("DCM", length(dcm_cols)),
  rep("HC", length(hc_cols))
)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = factor(group)
)

explained <- round(
  100 * summary(pca_res)$importance[2, 1:2],
  1
)

p_pca <- ggplot(
  pca_df,
  aes(PC1, PC2, color = Group)
) +
  stat_ellipse(level = 0.95) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA: DCM vs HC",
    x = paste0("PC1 (", explained[1], "%)"),
    y = paste0("PC2 (", explained[2], "%)")
  )

ggsave("results/DCM_vs_HC_PCA.png", p_pca)

# ==================================================
# Step 11. PLS-DA
# ==================================================

library(mixOmics)

X <- t(pca_data)

Y <- factor(group)

plsda_res <- plsda(X, Y, ncomp = 2)

scores <- plsda_res$variates$X

plsda_df <- data.frame(
  Dim1 = scores[, 1],
  Dim2 = scores[, 2],
  Group = Y
)

p_plsda <- ggplot(
  plsda_df,
  aes(Dim1, Dim2, color = Group)
) +
  stat_ellipse(level = 0.95) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PLS-DA: DCM vs HC")

ggsave("results/DCM_vs_HC_PLSDA.png", p_plsda)

# ==================================================
# Step 12. Differential expression analysis
# ==================================================

pvals <- apply(num_data, 1, function(x) {
  t.test(x[dcm_cols], x[hc_cols])$p.value
})

log2fc <- rowMeans(num_data[, dcm_cols]) -
  rowMeans(num_data[, hc_cols])

fc <- 2^log2fc

fdr <- p.adjust(pvals, method = "BH")

result_df <- data.frame(
  Protein_ID = protein_names,
  Description = descriptions,
  log2FC = log2fc,
  FC = fc,
  p_value = pvals,
  FDR = fdr
)

sig_proteins <- result_df %>%
  filter(FDR < 0.05 & (FC > 1.5 | FC < 0.67))

# ==================================================
# Step 13. Volcano plot
# ==================================================

result_df$negLogFDR <- -log10(result_df$FDR)

result_df$Group <- "NotSig"
result_df$Group[result_df$FDR < 0.05 &
                  result_df$log2FC > log2(1.5)] <- "Up"
result_df$Group[result_df$FDR < 0.05 &
                  result_df$log2FC < log2(0.67)] <- "Down"

p_volcano <- ggplot(
  result_df,
  aes(log2FC, negLogFDR, color = Group)
) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Volcano Plot: DCM vs HC")

ggsave("results/DCM_vs_HC_Volcano.png", p_volcano)

cat("Analysis completed.\n")