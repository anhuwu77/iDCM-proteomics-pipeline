############################################################
# Project: iDCM CMR-Proteomics Study
# Script: 03_PCA_Survival_KM_Cox_Analysis.R
#
# Purpose:
#   - Cox regression for PC1–PC4
#   - Kaplan–Meier survival analysis using optimal cutpoints
#   - Visualization of survival curves with confidence intervals
############################################################

# ==================================================
# Load required packages
# ==================================================

library(readxl)
library(survival)
library(survminer)

# ==================================================
# Load dataset
# ==================================================

data <- read_excel(
  "data/clinical_data_final.xlsx"
)

# ==================================================
# Preprocessing
# ==================================================

# Convert days to months
data$Time <- data$Time / 30.44

# Remove missing values in key variables
data <- data[!is.na(data$Time) &
               !is.na(data$Event), ]

# Build survival object
surv_obj <- Surv(time = data$Time, event = data$Event)

# ==================================================
# Cox regression (continuous PC variables)
# ==================================================

cox_model <- coxph(
  surv_obj ~ PC1 + PC2 + PC3 + PC4,
  data = data
)

cat("\n=== Cox regression results ===\n")
print(summary(cox_model))

# ==================================================
# Function: KM analysis for each PC
# ==================================================

run_km <- function(data, pc_name) {
  
  # Optimal cutpoint
  cut <- surv_cutpoint(
    data,
    time = "Time",
    event = "Event",
    variables = pc_name
  )
  
  group <- surv_categorize(cut)[[pc_name]]
  
  # Ensure factor order
  group <- factor(group, levels = c("Low", "High"))
  
  data[[paste0(pc_name, "_group")]] <- group
  
  # KM model
  fit <- survfit(
    Surv(Time, Event) ~ group,
    data = data
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = data,
    conf.int = TRUE,
    pval = TRUE,
    risk.table = TRUE,
    title = paste0("Kaplan-Meier Curve: ", pc_name),
    legend.title = "",
    legend.labs = c(paste0(pc_name, " = Low"),
                    paste0(pc_name, " = High")),
    xlab = "Time (months)",
    ylab = "Event-free survival probability",
    break.time.by = 6,
    xlim = c(0, 36)
  )
  
  return(list(fit = fit, plot = p))
}

# ==================================================
# PC1–PC4 KM analysis
# ==================================================

res_PC1 <- run_km(data, "PC1")
res_PC2 <- run_km(data, "PC2")
res_PC3 <- run_km(data, "PC3")
res_PC4 <- run_km(data, "PC4")

# ==================================================
# Save plots
# ==================================================

ggsave("results/KM_PC1.png", res_PC1$plot$plot)
ggsave("results/KM_PC2.png", res_PC2$plot$plot)
ggsave("results/KM_PC3.png", res_PC3$plot$plot)
ggsave("results/KM_PC4.png", res_PC4$plot$plot)

# ==================================================
# End of script
# ==================================================