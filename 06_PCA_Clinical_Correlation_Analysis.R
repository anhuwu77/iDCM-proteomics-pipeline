############################################################
# Correlation analysis between PCA-derived traits
# and clinical variables (Spearman + Kruskal-Wallis)
############################################################

library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ==============================
# Load data
# ==============================

file_path <- "data/PCA_clinical_combined.xlsx"
df <- read_excel(file_path)

# ==============================
# Define variables
# ==============================

pc_vars <- c("PC1", "PC3", "PC4",
             "ProteinScore_PCA1",
             "ProteinScore_PCA3",
             "ProteinScore_PCA4")

clinical_vars <- c("Sex", "Age", "BMI", "SBP", "DBP", "HR",
                   "NYHA", "eGFR", "NP",
                   "diabetes", "hypertension",
                   "Atrial fibrillation/flutter",
                   "Ventricular arrhythmias", "LBBB")

df <- df[, c(pc_vars, clinical_vars)]

cat_vars <- c("Sex", "NYHA", "diabetes", "hypertension",
              "Atrial fibrillation/flutter",
              "Ventricular arrhythmias", "LBBB")

cont_vars <- setdiff(clinical_vars, cat_vars)

# ==============================
# Data formatting
# ==============================

df$Sex <- factor(df$Sex, levels = c(0, 1), labels = c("Male", "Female"))
df$diabetes <- factor(df$diabetes, levels = c(0, 1), labels = c("No", "Yes"))
df$hypertension <- factor(df$hypertension, levels = c(0, 1), labels = c("No", "Yes"))
df$`Atrial fibrillation/flutter` <- factor(df$`Atrial fibrillation/flutter`, levels = c(0, 1), labels = c("No", "Yes"))
df$`Ventricular arrhythmias` <- factor(df$`Ventricular arrhythmias`, levels = c(0, 1), labels = c("No", "Yes"))
df$LBBB <- factor(df$LBBB, levels = c(0, 1), labels = c("No", "Yes"))

df$NYHA <- factor(df$NYHA, levels = c(1,2,3,4),
                  labels = c("I","II","III","IV"),
                  ordered = TRUE)

df[cont_vars] <- lapply(df[cont_vars], as.numeric)

# ==============================
# 1. Spearman correlation (continuous)
# ==============================

cor_mat <- matrix(NA, length(pc_vars), length(cont_vars),
                  dimnames = list(pc_vars, cont_vars))

p_mat <- cor_mat

for (i in seq_along(pc_vars)) {
  for (j in seq_along(cont_vars)) {
    
    x <- df[[pc_vars[i]]]
    y <- df[[cont_vars[j]]]
    
    if (all(is.na(x)) | all(is.na(y))) next
    
    test <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    cor_mat[i, j] <- test$estimate
    p_mat[i, j] <- test$p.value
  }
}

stars <- ifelse(p_mat < 0.001, "***",
                ifelse(p_mat < 0.01, "**",
                       ifelse(p_mat < 0.05, "*", "")))

label_mat <- paste0(sprintf("%.2f", cor_mat), stars)

pdf("Spearman_PC_vs_Clinical.pdf", width = 10, height = 6)

Heatmap(cor_mat,
        name = "Spearman r",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.text(label_mat[i, j], x, y, gp = gpar(fontsize = 10))
        })

dev.off()

# ==============================
# 2. PC vs categorical variables
#    Kruskal-Wallis effect size (eta²)
# ==============================

eta_mat <- matrix(NA, length(pc_vars), length(cat_vars),
                  dimnames = list(pc_vars, cat_vars))

p_cat <- eta_mat

for (i in seq_along(pc_vars)) {
  for (j in seq_along(cat_vars)) {
    
    x <- df[[pc_vars[i]]]
    g <- df[[cat_vars[j]]]
    
    if (length(unique(g)) < 2) next
    
    kw <- kruskal.test(x ~ g)
    
    H <- as.numeric(kw$statistic)
    n <- length(x)
    
    eta2 <- H / (n - 1)
    
    eta_mat[i, j] <- eta2
    p_cat[i, j] <- kw$p.value
  }
}

cat_stars <- ifelse(p_cat < 0.001, "***",
                    ifelse(p_cat < 0.01, "**",
                           ifelse(p_cat < 0.05, "*", "")))

eta_label <- paste0(sprintf("%.2f", eta_mat), cat_stars)

pdf("Kruskal_PC_vs_Categorical.pdf", width = 10, height = 6)

Heatmap(eta_mat,
        name = "Eta²",
        col = colorRamp2(c(0, 0.1, 0.3, 1),
                         c("white", "pink", "red", "darkred")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.text(eta_label[i, j], x, y, gp = gpar(fontsize = 10))
        })

dev.off()

cat("\nCorrelation analysis completed.\n")