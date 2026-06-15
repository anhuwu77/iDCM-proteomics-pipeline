############################################################
# LASSO regression for PCA-derived phenotypes
# Protein signature construction
############################################################

library(readxl)
library(glmnet)
library(tidyverse)

# =========================
# Path
# =========================

save_path <- "results/lasso/"
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load data
# =========================

data <- read_excel("data/clinical_proteomics_lasso_input.xlsx")

id_col <- data[[1]]
pca_scores <- data[, 2:5]
protein_data <- data[, 6:(5 + 1309)]

x <- as.matrix(protein_data)

# =========================
# Output containers
# =========================

lasso_results <- list()

protein_scores <- data.frame(ID = id_col)

# =========================
# Main loop
# =========================

for (i in 1:4) {
  
  cat("\nProcessing PCA", i, "\n")
  
  y <- as.numeric(pca_scores[[i]])
  
  # remove NA pairs
  valid_idx <- complete.cases(y)
  y <- y[valid_idx]
  x_sub <- x[valid_idx, ]
  
  set.seed(123)
  
  # =========================
  # Cross-validated LASSO
  # =========================
  
  cv_model <- cv.glmnet(
    x_sub,
    y,
    alpha = 1,
    standardize = TRUE
  )
  
  lambda_min <- cv_model$lambda.min
  
  final_model <- glmnet(
    x_sub,
    y,
    alpha = 1,
    lambda = lambda_min
  )
  
  # =========================
  # Extract coefficients
  # =========================
  
  coef_mat <- as.matrix(coef(final_model))
  
  coef_df <- data.frame(
    Protein = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat)
  ) %>%
    filter(Protein != "(Intercept)",
           Coefficient != 0)
  
  lasso_results[[paste0("PCA", i)]] <- coef_df
  
  write.csv(
    coef_df,
    paste0(save_path, "LASSO_PCA", i, "_coefficients.csv"),
    row.names = FALSE
  )
  
  # =========================
  # Protein score calculation
  # =========================
  
  predicted_score <- predict(
    final_model,
    newx = x_sub,
    s = lambda_min
  )
  
  protein_scores[[paste0("ProteinScore_PCA", i)]] <- NA
  protein_scores[valid_idx, paste0("ProteinScore_PCA", i)] <- as.numeric(predicted_score)
  
  # =========================
  # Plot: CV curve
  # =========================
  
  pdf(paste0(save_path, "LASSO_PCA", i, "_CV_plot.pdf"))
  
  plot(cv_model)
  abline(v = log(lambda_min), col = "red", lty = 2)
  
  n_selected <- nrow(coef_df)
  
  title(main = paste0(
    "PCA", i,
    " LASSO (Selected proteins = ",
    n_selected, ")"
  ))
  
  dev.off()
  
  # =========================
  # Plot: coefficient path
  # =========================
  
  pdf(paste0(save_path, "LASSO_PCA", i, "_Path_plot.pdf"))
  
  plot(cv_model$glmnet.fit, xvar = "lambda")
  abline(v = log(lambda_min), col = "red", lty = 2)
  
  title(main = paste0("LASSO Path - PCA", i))
  
  dev.off()
}

# =========================
# Save outputs
# =========================

write.csv(
  protein_scores,
  paste0(save_path, "Patient_Protein_Score_LASSO.csv"),
  row.names = FALSE
)

cat("\nLASSO analysis completed.\n")