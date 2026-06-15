############################################################
# Project: Integrated CMR and Plasma Proteomics in iDCM
# Script: 01_CMR_Phenotyping_PCA.R
#
# Purpose:
#   - Assess variable distributions
#   - Perform log-transformation when necessary
#   - Standardize CMR variables
#   - Conduct PCA with Varimax rotation
#   - Generate correlation heatmap
#   - Generate scree plot
#   - Export PCA loadings and PC scores
############################################################

# ==================================================
# Step 0. Install and load required packages
# ==================================================

packages <- c(
  "readxl",
  "e1071",
  "ggcorrplot",
  "psych",
  "dplyr",
  "reshape2",
  "ggplot2"
)

lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

# ==================================================
# Step 1. Import CMR dataset
# ==================================================

file_path <- "data/CMR_dataset.xlsx"

data <- read_excel(file_path)

# ==================================================
# Step 2. Select 27 CMR variables for PCA
# ==================================================

vars <- c(
  "LVEDVI", "LVESVI", "LVEF", "LVMI",
  "RVEDVI", "RVESVI", "RVEF", "RVMI",
  "LVGLS", "LVGCS", "LVGRS",
  "RVGLS", "RVGCS", "RVGRS",
  "LAVImax", "LAEFtotal", "sLA", "eLA", "aLA",
  "RAVImax", "RAEFtotal", "sRA", "eRA", "aRA",
  "preT1", "ECV", "LGEmass"
)

pca_data <- data[, vars]

# ==================================================
# Step 3. Skewness assessment and log1p transformation
# ==================================================

negative_cols <- sapply(
  pca_data,
  function(x) any(x < 0, na.rm = TRUE)
)

cat("Variables containing negative values:\n")
print(names(negative_cols[negative_cols]))

positive_data <- pca_data[, !negative_cols]

skew_vals <- sapply(
  positive_data,
  function(x) {
    e1071::skewness(
      x,
      na.rm = TRUE,
      type = 2
    )
  }
)

skewed_vars <- names(skew_vals[skew_vals > 1])

cat("Variables selected for log1p transformation (skewness > 1):\n")
print(skewed_vars)

if (length(skewed_vars) > 0) {
  pca_data[, skewed_vars] <- log1p(pca_data[, skewed_vars])
}

write.csv(
  data.frame(
    Variable = names(skew_vals),
    Skewness = skew_vals
  ),
  file = "results/PCA_skewness_check.csv",
  row.names = FALSE
)

# ==================================================
# Step 4. Z-score normalization
# ==================================================

pca_scaled <- as.data.frame(scale(pca_data))

# ==================================================
# Step 5. Correlation matrix visualization
# ==================================================

corr_matrix <- cor(
  pca_scaled,
  use = "pairwise.complete.obs"
)

p_corr <- ggcorrplot(
  corr_matrix,
  method = "square",
  type = "full",
  lab = TRUE,
  lab_size = 3,
  colors = c("blue", "white", "red"),
  outline.col = "gray"
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      angle = 0
    )
  )

print(p_corr)

ggsave(
  filename = "results/Correlation_Heatmap_Full.png",
  plot = p_corr,
  width = 10,
  height = 10,
  dpi = 300
)

# ==================================================
# Step 6. Principal Component Analysis
# (4 Components with Varimax Rotation)
# ==================================================

pca_result <- psych::principal(
  pca_scaled,
  nfactors = 4,
  rotate = "varimax",
  scores = TRUE
)

# ==================================================
# Step 7. Variance Explained and Scree Plot
# ==================================================

eig_values <- eigen(cor(pca_scaled))$values

var_explained <- eig_values / sum(eig_values) * 100
cum_var <- cumsum(var_explained)

df <- data.frame(
  PC = 1:length(eig_values),
  Variance = var_explained,
  Cumulative = cum_var
)

p_scree <- ggplot(df, aes(x = PC)) +
  geom_bar(
    aes(y = Variance),
    stat = "identity",
    fill = "steelblue"
  ) +
  geom_line(
    aes(y = Cumulative),
    color = "red",
    linewidth = 1
  ) +
  geom_point(
    aes(y = Cumulative),
    color = "red",
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(eig_values)
  ) +
  ylab("Variance Explained (%)") +
  xlab("Principal Components") +
  ggtitle(
    "Scree Plot: Variance Explained and Cumulative Variance"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )

print(p_scree)

ggsave(
  filename = "results/Scree_Plot_CMR_PCA.png",
  plot = p_scree,
  width = 8,
  height = 6,
  dpi = 300
)

# ==================================================
# Step 8. Export PCA Loadings
# ==================================================

loadings_df <- as.data.frame(
  unclass(pca_result$loadings)
)

write.csv(
  loadings_df,
  file = "results/PCA_Loadings_4PC.csv",
  row.names = TRUE
)

# ==================================================
# Step 9. Export Principal Component Scores
# ==================================================

pc_scores <- as.data.frame(
  pca_result$scores
)

final_data <- cbind(
  data,
  pc_scores
)

write.csv(
  final_data,
  file = "results/CMR_PCA_Results_With_PC_Scores.csv",
  row.names = FALSE
)

# ==================================================
# Step 10. PCA Loading Heatmap
# ==================================================

loadings_df$Variable <- rownames(loadings_df)

loadings_df_sub <- loadings_df[
  ,
  c("RC1", "RC2", "RC3", "RC4", "Variable")
]

loadings_long <- reshape2::melt(
  loadings_df_sub,
  id.vars = "Variable",
  variable.name = "PC",
  value.name = "Loading"
)

loadings_long$AbsLoading <- abs(
  loadings_long$Loading
)

loadings_long$PC <- gsub(
  "RC",
  "PC",
  loadings_long$PC
)

p_loading <- ggplot(
  loadings_long,
  aes(
    x = PC,
    y = Variable,
    fill = AbsLoading
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(Loading, 2)),
    size = 3,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#00BFFF",
    mid = "white",
    high = "#FF4500",
    midpoint = 0.4,
    limits = c(
      0,
      max(loadings_long$AbsLoading)
    )
  ) +
  labs(
    title = "PCA Varimax Loadings Heatmap (PC1-PC4)",
    x = "Principal Components",
    y = "Variables",
    fill = "Absolute Loading"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )

print(p_loading)

ggsave(
  filename = "results/PCA_Loadings_Heatmap_4PC.png",
  plot = p_loading,
  width = 10,
  height = 8,
  dpi = 300
)

# ==================================================
# End of Script
# ==================================================