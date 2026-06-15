############################################
# Phenotype clustering based on Cox-derived risk score
# K-means clustering + KM analysis + cluster validation
############################################

# =============================
# Load packages
# =============================
library(readxl)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(cluster)
library(factoextra)
library(writexl)

# =============================
# Read data
# =============================
data <- read_excel("DATA_FILE.xlsx")
names(data) <- make.names(names(data))

# =============================
# Cox-derived risk score
# =============================
surv_obj <- Surv(time = data$Time, event = data$Event)

cox_model <- coxph(surv_obj ~ PC1 + PC2 + PC3 + PC4, data = data)
data$risk_score <- predict(cox_model, type = "lp")

# =============================
# K-means clustering (K = 2)
# =============================
set.seed(123)

kmeans_model <- kmeans(scale(data$risk_score), centers = 2, nstart = 50)

data$cluster <- factor(kmeans_model$cluster,
                       labels = c("Cluster1", "Cluster2"))

# =============================
# Cluster summary
# =============================
cluster_summary <- data.frame(
  Metric = c("Cluster1_size", "Cluster2_size",
             "Within_SS", "Between_SS", "Total_SS", "Explained_variance"),
  Value = c(kmeans_model$size[1],
            kmeans_model$size[2],
            kmeans_model$tot.withinss,
            kmeans_model$betweenss,
            kmeans_model$totss,
            kmeans_model$betweenss / kmeans_model$totss)
)

print(cluster_summary)

# =============================
# Rename phenotype groups
# =============================
data$Phenotype <- recode(data$cluster,
                         "Cluster1" = "Advanced remodeling subtype",
                         "Cluster2" = "Compensated function subtype")

# =============================
# Survival analysis (KM)
# =============================
data$Time_months <- data$Time / 30.44

fit <- survfit(Surv(Time_months, Event) ~ Phenotype, data = data)

km_plot <- ggsurvplot(
  fit,
  data = data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Event-free survival probability",
  legend.title = "",
  legend.labs = levels(factor(data$Phenotype)),
  break.time.by = 6,
  xlim = c(0, 36)
)

ggsave("KM_Phenotype_Clusters.png", km_plot$plot, width = 6, height = 5)

# =============================
# PCA profile of clusters
# =============================
pc_summary <- data %>%
  group_by(Phenotype) %>%
  summarise(across(starts_with("PC"), mean), .groups = "drop")

pc_long <- pc_summary %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "PC",
               values_to = "MeanScore")

bar_plot <- ggplot(pc_long, aes(x = PC, y = MeanScore, fill = PC)) +
  geom_col(width = 0.6) +
  facet_wrap(~Phenotype, ncol = 1) +
  geom_hline(yintercept = 0) +
  labs(title = "PCA profile across phenotype clusters",
       x = "PCA axis",
       y = "Mean score") +
  theme_minimal()

ggsave("PC_Profile_by_Phenotype.png", bar_plot, width = 6, height = 6)

# =============================
# Elbow method
# =============================
wss <- sapply(1:6, function(k) {
  kmeans(scale(data$risk_score), centers = k, nstart = 20)$tot.withinss
})

elbow_df <- data.frame(K = 1:6, WSS = wss)

elbow_plot <- ggplot(elbow_df, aes(x = K, y = WSS)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow method for optimal K",
       x = "K",
       y = "Within-cluster SS") +
  theme_minimal()

ggsave("Elbow_plot.png", elbow_plot, width = 6, height = 4)

# =============================
# Silhouette analysis
# =============================
dist_mat <- dist(scale(data$risk_score))
sil <- silhouette(kmeans_model$cluster, dist_mat)

sil_plot <- fviz_silhouette(sil) +
  ggtitle("Silhouette analysis (K = 2)")

ggsave("Silhouette_plot.png", sil_plot, width = 6, height = 4)

# =============================
# Risk score distribution
# =============================
density_plot <- ggplot(data, aes(x = risk_score, fill = Phenotype)) +
  geom_density(alpha = 0.5) +
  labs(title = "Risk score distribution by phenotype",
       x = "Risk score",
       y = "Density") +
  theme_minimal()

ggsave("RiskScore_Distribution.png", density_plot, width = 6, height = 4)

# =============================
# Export results
# =============================
centroids <- data %>%
  group_by(Phenotype) %>%
  summarise(across(PC1:PC4, mean), .groups = "drop")

write_xlsx(
  list(
    Cluster_Summary = cluster_summary,
    Centroids = centroids
  ),
  "Kmeans_Phenotype_Results.xlsx"
)

cat("Analysis completed.\n")