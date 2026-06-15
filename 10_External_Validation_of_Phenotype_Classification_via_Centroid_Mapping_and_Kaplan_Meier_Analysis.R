############################################################
# 10 External Validation of Phenotype Classification via Centroid Mapping and Kaplan–Meier Analysis
############################################################

# ==================================================
# Step 0. Load required packages
# ==================================================

packages <- c(
  "openxlsx",
  "dplyr",
  "survival",
  "survminer"
)

lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

# ==================================================
# Step 1. Read datasets
# ==================================================

df_discovery <- openxlsx::read.xlsx("data/discovery_centroids.xlsx")
df_validation <- openxlsx::read.xlsx("data/validation_pc_scores.xlsx")

# ==================================================
# Step 2. Compute centroids (PC1–PC4 mean by phenotype)
# ==================================================

centroids <- df_discovery %>%
  group_by(PhenotypeGroup) %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    PC3 = mean(PC3, na.rm = TRUE),
    PC4 = mean(PC4, na.rm = TRUE)
  )

print("Discovery cohort centroids:")
print(centroids)

# ==================================================
# Step 3. Euclidean distance function
# ==================================================

euclid <- function(x, center) sqrt(sum((x - center)^2))

center1 <- as.numeric(centroids[1, 2:5])
center2 <- as.numeric(centroids[2, 2:5])

# ==================================================
# Step 4. Assign phenotype in validation cohort
# ==================================================

df_validation$dist_to_1 <- apply(
  df_validation[, c("PC1_score", "PC2_score", "PC3_score", "PC4_score")],
  1,
  function(x) euclid(as.numeric(x), center1)
)

df_validation$dist_to_2 <- apply(
  df_validation[, c("PC1_score", "PC2_score", "PC3_score", "PC4_score")],
  1,
  function(x) euclid(as.numeric(x), center2)
)

df_validation$PredictedSubtype <- ifelse(
  df_validation$dist_to_1 < df_validation$dist_to_2,
  centroids$PhenotypeGroup[1],
  centroids$PhenotypeGroup[2]
)

# ==================================================
# Step 5. Export validation classification results
# ==================================================

df_out <- df_validation %>%
  dplyr::select(
    ID,
    PC1_score,
    PC2_score,
    PC3_score,
    PC4_score,
    dist_to_1,
    dist_to_2,
    PredictedSubtype
  )

openxlsx::write.xlsx(
  df_out,
  file = "External_Validation_Phenotype_Results.xlsx",
  rowNames = FALSE
)

cat("Validation classification completed.\n")

# ==================================================
# Step 6. Kaplan–Meier survival analysis
# ==================================================

df_val_km <- df_out

df_val_km$Time <- as.numeric(df_val_km$Time)
df_val_km$Event <- as.numeric(df_val_km$Event)

df_val_km$PredictedSubtype <- factor(
  df_val_km$PredictedSubtype,
  levels = centroids$PhenotypeGroup
)

surv_obj <- Surv(time = df_val_km$Time, event = df_val_km$Event)

fit_km <- survfit(surv_obj ~ PredictedSubtype, data = df_val_km)

km_plot <- ggsurvplot(
  fit_km,
  data = df_val_km,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Survival probability",
  legend.title = "Phenotype",
  legend.labs = levels(df_val_km$PredictedSubtype),
  palette = c("#E64B35", "#4DBBD5"),
  surv.median.line = "hv",
  ggtheme = theme_bw(base_size = 14)
)

print(km_plot)

ggsave(
  filename = "External_Validation_KM_Curve.png",
  plot = km_plot$plot,
  width = 6,
  height = 5,
  dpi = 300
)

cat("KM curve saved.\n")

# ==================================================
# Step 7. Follow-up time summary
# ==================================================

median_followup <- median(df_val_km$Time, na.rm = TRUE)
q_followup <- quantile(df_val_km$Time, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

cat("Follow-up summary (days):\n")
cat("Median:", median_followup, "\n")
cat("Q1:", q_followup[1], "\n")
cat("Q3:", q_followup[3], "\n")

# ==================================================
# End of script
# ==================================================