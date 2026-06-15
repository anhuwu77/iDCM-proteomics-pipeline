############################################################
# Machine Learning Feature Selection per PCA axis
# RF + XGBoost + Union visualization
############################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(randomForest)
library(xgboost)
library(openxlsx)

# =========================
# Output path
# =========================

save_path <- "results/mechanism_proteins/"
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# =========================
# Load data
# =========================

data <- read_excel("data/clinical_proteomics_pca_input.xlsx")

pca_scores <- data[, 2:5]
x_all <- data[, 6:1314] %>% mutate_all(as.numeric)

wb <- createWorkbook()

# =========================
# Loop PCA axes
# =========================

for (i in 1:4) {
  
  axis_name <- paste0("PCA", i)
  cat("\nProcessing:", axis_name, "\n")
  
  # phenotype (PCA score)
  y <- as.numeric(pca_scores[[i]])
  x <- as.matrix(x_all)
  
  # remove NA
  valid_idx <- complete.cases(y)
  y <- y[valid_idx]
  x <- x[valid_idx, ]
  
  # =========================
  # XGBoost
  # =========================
  
  dtrain <- xgb.DMatrix(data = x, label = y)
  
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse"
  )
  
  set.seed(123)
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 200,
    verbose = 0
  )
  
  xgb_importance <- xgb.importance(model = xgb_model)
  
  xgb_top <- head(xgb_importance$Feature, 100)
  
  # =========================
  # Random Forest
  # =========================
  
  set.seed(123)
  rf_model <- randomForest(
    x = x,
    y = y,
    importance = TRUE,
    ntree = 1000
  )
  
  rf_importance <- importance(rf_model, type = 1)
  rf_ranked <- sort(rf_importance[, 1], decreasing = TRUE)
  
  rf_top <- names(rf_ranked)[1:100]
  
  # =========================
  # UNION proteins
  # =========================
  
  union_proteins <- union(rf_top, xgb_top)
  
  cat("Union size:", length(union_proteins), "\n")
  
  addWorksheet(wb, axis_name)
  writeData(wb, axis_name,
            data.frame(Protein = union_proteins))
  
  # =========================
  # Prepare UNION plot data
  # =========================
  
  xgb_df <- xgb_importance %>%
    filter(Feature %in% union_proteins) %>%
    select(Feature, Gain) %>%
    mutate(Source = ifelse(Feature %in% rf_top, "Both", "XGBoost"))
  
  rf_df <- data.frame(
    Feature = names(rf_ranked),
    Importance = as.numeric(rf_ranked)
  ) %>%
    filter(Feature %in% union_proteins) %>%
    mutate(Source = ifelse(Feature %in% xgb_top, "Both", "RF")) %>%
    rename(Gain = Importance)
  
  combined_df <- bind_rows(xgb_df, rf_df) %>%
    group_by(Feature, Source) %>%
    summarise(Importance = mean(Gain), .groups = "drop")
  
  combined_df <- combined_df %>%
    arrange(desc(Importance)) %>%
    mutate(Feature = factor(Feature,
                            levels = rev(unique(Feature))))
  
  # =========================
  # UNION plot
  # =========================
  
  union_plot <- ggplot(
    combined_df,
    aes(x = Feature, y = Importance, fill = Source)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(
      "RF" = "#56B4E9",
      "XGBoost" = "#E69F00",
      "Both" = "#009E73"
    )) +
    theme_bw() +
    labs(
      title = paste0(axis_name, " - RF & XGBoost Union Proteins"),
      x = "Protein",
      y = "Importance",
      fill = "Source"
    )
  
  ggsave(
    filename = paste0(save_path, axis_name, "_UNION_plot.png"),
    plot = union_plot,
    width = 8,
    height = 12,
    dpi = 300
  )
}

# =========================
# Save Excel
# =========================

save_file <- paste0(save_path,
                    "PCA_RF_XGB_union.xlsx")

saveWorkbook(wb, save_file, overwrite = TRUE)

cat("\nDONE. Results saved.\n")