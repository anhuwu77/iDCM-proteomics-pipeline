############################################################
# Project: ELISA-Based 5-Protein Prognostic Model in iDCM
# Script: 11_ELISA_5Protein_Cox_KM.R
#
# Purpose:
#   - Build multivariable Cox model using 5 ELISA proteins
#   - Construct protein-based risk score
#   - Perform survival stratification (median split)
#   - Generate Kaplan–Meier curves
#   - Calculate HR and log-rank p value
#   - Export survival plots
############################################################

# ============================================================
# Step 0. Load required packages
# ============================================================
library(readxl)
library(dplyr)
library(survival)
library(survminer)

# ============================================================
# Step 1. Import dataset
# ============================================================
file_path <- "data/ELISA_protein_survival_data.xlsx"
data <- read_excel(file_path)

# ============================================================
# Step 2. Define selected proteins (5-protein signature)
# ============================================================
selected_proteins <- c(
  "ALDH6A1",
  "ANP",
  "MYL3",
  "PKIA",
  "PPIA"
)

# ============================================================
# Step 3. Construct Cox dataset
# ============================================================
cox_data <- data %>%
  dplyr::select(all_of(selected_proteins), Time, Event) %>%
  na.omit()

# ============================================================
# Step 4. Multivariable Cox regression
# ============================================================
cox_model <- coxph(Surv(Time, Event) ~ ., data = cox_data)
summary(cox_model)

cox_coef <- coef(cox_model)

# ============================================================
# Step 5. Construct protein risk score
# ============================================================
cox_data$ProteinScore <- as.matrix(cox_data[, selected_proteins]) %*% cox_coef

# map back to original dataset
data$ProteinScore <- NA
data$ProteinScore[as.numeric(rownames(cox_data))] <- cox_data$ProteinScore

# ============================================================
# Step 6. Median split risk grouping
# ============================================================
median_score <- median(cox_data$ProteinScore, na.rm = TRUE)

cox_data$RiskGroup <- ifelse(
  cox_data$ProteinScore >= median_score,
  "High",
  "Low"
)

cox_data$RiskGroup <- factor(cox_data$RiskGroup, levels = c("Low", "High"))

# ============================================================
# Step 7. Time conversion (days → months)
# ============================================================
cox_data$Time_month <- cox_data$Time / 30.44

# ============================================================
# Step 8. Kaplan–Meier + Cox HR + log-rank test
# ============================================================
fit <- survfit(Surv(Time_month, Event) ~ RiskGroup, data = cox_data)

km_cox <- coxph(Surv(Time_month, Event) ~ RiskGroup, data = cox_data)
summary_km <- summary(km_cox)

hr <- summary_km$coefficients[1, "exp(coef)"]
lower_ci <- summary_km$conf.int[1, "lower .95"]
upper_ci <- summary_km$conf.int[1, "upper .95"]

hr_text <- sprintf("HR = %.2f (95%% CI: %.2f - %.2f)", hr, lower_ci, upper_ci)

logrank <- survdiff(Surv(Time_month, Event) ~ RiskGroup, data = cox_data)
p_val <- 1 - pchisq(logrank$chisq, df = 1)
p_text <- sprintf("Log-rank p = %.3f", p_val)

# ============================================================
# Step 9. KM plot
# ============================================================
km_plot <- ggsurvplot(
  fit,
  data = cox_data,
  pval = FALSE,
  conf.int = TRUE,
  risk.table = TRUE,
  break.time.by = 6,
  xlim = c(0, 36),
  palette = c("#377EB8", "#E41A1C"),
  legend.title = "Protein Risk Group",
  legend.labs = c("Low", "High"),
  title = "Kaplan–Meier Curve Based on 5-Protein Signature",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  risk.table.height = 0.25
)

km_plot$plot <- km_plot$plot +
  annotate("text", x = 10, y = 0.35, label = hr_text, size = 5) +
  annotate("text", x = 10, y = 0.25, label = p_text, size = 5)

print(km_plot)

# ============================================================
# Step 10. Export figures
# ============================================================
output_dir <- "results/"

ggsave(
  filename = paste0(output_dir, "KM_5Protein_Cox.png"),
  plot = km_plot$plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = paste0(output_dir, "KM_5Protein_Cox.tiff"),
  plot = km_plot$plot,
  width = 6,
  height = 5,
  dpi = 300
)

cat("Analysis completed: 5-protein survival model finished.\n")