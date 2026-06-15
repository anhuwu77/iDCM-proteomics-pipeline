############################################################
# 09 LASSO-Derived Plasma Protein Signature for Clinical Phenotype Classification
############################################################

# ==================================================
# Step 0. Load required packages
# ==================================================

packages <- c(
  "readxl",
  "tidyverse",
  "glmnet",
  "pROC",
  "caret",
  "writexl",
  "tibble"
)

lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
})

# ==================================================
# Step 1. Import dataset
# ==================================================

file_path <- "data/phenotype_protein_dataset.xlsx"
data <- read_excel(file_path)

# protein matrix (adjust column range if needed)
protein_matrix <- as.matrix(data[, 6:(6 + 1309 - 1)])

# phenotype label (K-means derived clinical phenotype)
data$group_binary <- ifelse(
  data$PhenotypeGroup2 == "Advanced Remodeling subtype",
  1, 0
)

# ==================================================
# Step 2. LASSO regression for feature selection
# ==================================================

set.seed(123)

cvfit <- cv.glmnet(
  x = protein_matrix,
  y = data$group_binary,
  alpha = 1,
  family = "binomial",
  standardize = TRUE,
  nfolds = 10
)

best_lambda <- cvfit$lambda.min

lasso_model <- glmnet(
  x = protein_matrix,
  y = data$group_binary,
  alpha = 1,
  family = "binomial",
  lambda = best_lambda
)

# ==================================================
# Step 3. Extract LASSO-selected proteins
# ==================================================

coef_mat <- as.matrix(coef(lasso_model))

coef_df <- data.frame(
  Protein = rownames(coef_mat),
  Coefficient = as.numeric(coef_mat)
)

lasso_selected <- coef_df %>%
  filter(Coefficient != 0 & Protein != "(Intercept)") %>%
  arrange(desc(abs(Coefficient)))

write_xlsx(
  lasso_selected,
  "LASSO_Selected_Proteins_Phenotype.xlsx"
)

# ==================================================
# Step 4. Build phenotype classification model
# ==================================================

selected_proteins <- lasso_selected$Protein

model_data <- as.data.frame(
  protein_matrix[, selected_proteins, drop = FALSE]
)

model_data$group_binary <- data$group_binary

glm_model <- glm(group_binary ~ ., data = model_data, family = binomial)

pred_prob <- predict(glm_model, type = "response")

# ==================================================
# Step 5. ROC analysis
# ==================================================

roc_obj <- roc(data$group_binary, pred_prob)

auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

plot(
  roc_obj,
  print.auc = FALSE,
  col = "#E41A1C",
  lwd = 2,
  main = "ROC Curve: LASSO-Derived Phenotype Signature"
)

auc_label <- sprintf(
  "AUC = %.3f (95%% CI: %.3f - %.3f)",
  auc_value, auc_ci[1], auc_ci[3]
)

legend(
  "bottomright",
  legend = auc_label,
  bty = "n",
  text.col = "#E41A1C",
  cex = 0.9
)

# ==================================================
# Step 6. Save prediction results
# ==================================================

data$PredProb <- pred_prob
data$PredClass <- ifelse(pred_prob > 0.5, 1, 0)

write_xlsx(
  data,
  "Phenotype_LASSO_Prediction_Results.xlsx"
)
# ==================================================
# Step 7. Prognostic value of LASSO-derived signature
# ==================================================

data$SignatureScore <- pred_prob

median_score <- median(data$SignatureScore, na.rm = TRUE)

data$RiskGroup <- ifelse(
  data$SignatureScore >= median_score,
  "High",
  "Low"
)

data$RiskGroup <- factor(data$RiskGroup, levels = c("Low", "High"))

surv_obj <- Surv(data$Time, data$Event)

fit <- survfit(surv_obj ~ RiskGroup, data = data)

cox_model <- coxph(surv_obj ~ RiskGroup, data = data)
cox_sum <- summary(cox_model)

hr <- cox_sum$coefficients[1, "exp(coef)"]
ci_low <- cox_sum$conf.int[1, "lower .95"]
ci_up <- cox_sum$conf.int[1, "upper .95"]

hr_text <- sprintf("HR = %.2f (95%% CI: %.2f - %.2f)", hr, ci_low, ci_up)

logrank <- survdiff(surv_obj ~ RiskGroup, data = data)
p_val <- 1 - pchisq(logrank$chisq, df = 1)

p_text <- sprintf("Log-rank p = %.3f", p_val)

km_plot <- ggsurvplot(
  fit,
  data = data,
  risk.table = TRUE,
  pval = FALSE,
  conf.int = TRUE,
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Risk Group",
  legend.labs = c("Low", "High"),
  palette = c("#377EB8", "#E41A1C")
)

km_plot$plot <- km_plot$plot +
  annotate("text", x = 10, y = 0.3, label = hr_text, size = 5) +
  annotate("text", x = 10, y = 0.2, label = p_text, size = 5)

print(km_plot)

# ==================================================
# Step 8. End
# ==================================================

cat("Analysis completed: LASSO phenotype signature + prognostic validation finished.\n")