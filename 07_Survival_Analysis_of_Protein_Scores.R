############################################
# Survival Analysis of Protein Scores
# Kaplan–Meier and Cox regression
# Adjusted for Age, Sex, BMI, and NT-proBNP
############################################

# =============================
# Load packages
# =============================
library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(writexl)

# =============================
# Data input
# =============================
data <- read_excel("DATA_FILE.xlsx")
names(data) <- make.names(names(data))

data$Sex <- factor(data$Sex, levels = c(0,1), labels = c("Male", "Female"))

# Protein score variables
pc_vars <- c("ProteinScore_PCA1", "ProteinScore_PCA3", "ProteinScore_PCA4")

# =============================
# Output container
# =============================
results_list <- list()

# =============================
# Main analysis loop
# =============================
for (pc in pc_vars) {
  
  if (!(pc %in% names(data))) next
  
  cat("\nProcessing:", pc, "\n")
  
  # =============================
  # Optimal cutpoint grouping
  # =============================
  cut_res <- surv_cutpoint(
    data,
    time = "Time",
    event = "Event",
    variables = pc
  )
  
  group_opt <- surv_categorize(cut_res)[[pc]]
  group_opt <- factor(group_opt, levels = c("low", "high"), labels = c("Low", "High"))
  
  data$group_opt <- group_opt
  
  # Safety check
  if (length(unique(group_opt)) < 2) next
  
  # =============================
  # Kaplan–Meier analysis
  # =============================
  fit_km <- survfit(Surv(Time, Event) ~ group_opt, data = data)
  
  km_plot <- ggsurvplot(
    fit_km,
    data = data,
    conf.int = TRUE,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time",
    ylab = "Event-free survival probability",
    legend.title = "",
    legend.labs = c("Low", "High")
  )
  
  ggsave(
    filename = paste0("KM_", pc, ".png"),
    plot = km_plot$plot,
    width = 8,
    height = 6
  )
  
  # =============================
  # Cox regression
  # =============================
  
  # Unadjusted model
  cox_unadj <- coxph(Surv(Time, Event) ~ group_opt, data = data)
  s_unadj <- summary(cox_unadj)
  
  hr_unadj <- exp(coef(cox_unadj))
  ci_unadj <- exp(confint(cox_unadj))
  p_unadj <- s_unadj$coefficients[, "Pr(>|z|)"]
  
  # Adjusted model
  cox_adj <- coxph(Surv(Time, Event) ~ group_opt + Age + Sex + BMI + NT_proBNP, data = data)
  s_adj <- summary(cox_adj)
  
  hr_adj <- exp(coef(cox_adj))["group_optHigh"]
  ci_adj <- exp(confint(cox_adj))["group_optHigh", ]
  p_adj <- s_adj$coefficients["group_optHigh", "Pr(>|z|)"]
  
  # =============================
  # Store results
  # =============================
  results_list[[pc]] <- data.frame(
    Variable = pc,
    Model = c("Unadjusted", "Adjusted"),
    HR = c(as.numeric(hr_unadj), as.numeric(hr_adj)),
    CI_Lower = c(ci_unadj[1,1], ci_adj[1]),
    CI_Upper = c(ci_unadj[1,2], ci_adj[2]),
    P_value = c(p_unadj, p_adj)
  )
}

# =============================
# Export results
# =============================
final_results <- bind_rows(results_list)

write_xlsx(final_results, "ProteinScore_KM_Cox_Results.xlsx")

cat("\nAnalysis completed.\n")