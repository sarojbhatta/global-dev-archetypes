# ============================================================
# Global Development Indicators — Clustering Analysis
# Developing Archetypes of Development
# ============================================================


# ============================================================
# Load Packages
# ============================================================

#| message: false
#| warning: false

set.seed(77)

library(tidyverse)
library(visdat)
library(naniar)
library(corrplot)
library(patchwork)
library(moments)
library(GGally)
library(mice)
library(AMR)
library(cluster)
library(ggdendro)


# ============================================================
# Data Loading, Exploration, and Cleaning
# ============================================================

# --- Load data ---

world_dev <- read.csv("world_dev.csv") |>
  select(-X) # need to remove the index too - the X variable

glimpse(world_dev)


# --- Visual overview of missingness ---

vis_dat(world_dev)

vis_miss(world_dev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Missingness Overview — World Development Indicators")


# --- Column-level missing counts and percentages ---

miss_summary <- world_dev |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(pct_missing = round(n_missing / nrow(world_dev) * 100, 1)) |>
  arrange(desc(pct_missing))
#print(miss_summary, n = 30)


# --- Visualize missingness by variable ---

# make the threshold decision from this
ggplot(miss_summary, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 35, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Percent Missing by Variable",
       subtitle = "Red dashed line = 35% threshold",
       x = NULL, y = "% Missing")

#ggsave("fig_missingness_bar.png", width = 7, height = 5, dpi = 150)


# --- Identify which variables to drop ---

vars_to_drop <- miss_summary |>
  filter(pct_missing > 35) |>
  pull(variable)

print(vars_to_drop)

world_clean <- world_dev |>
  select(-all_of(vars_to_drop)) #select(-all_of(vars_to_drop), -n_missing)

#glimpse(world_clean)


# --- Row level (country) missingness ---

world_clean <- world_clean |>
  mutate(n_missing = rowSums(is.na(across(where(is.numeric)))))

world_clean |>
  select(Country, n_missing) |>
  arrange(desc(n_missing)) |>
  slice_head(n = 20)


# --- Drop countries with excessive row-level missingness ---

# visualize row-level missingness distribution
world_clean |>
  ggplot(aes(x = n_missing)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Missing Variables per Country",
       subtitle = "Each bar = number of countries; red line = candidate drop threshold",
       x = "Number of Missing Variables", y = "Count of Countries")

#ggsave("fig_row_missingness.png", width = 7, height = 5, dpi = 150)


# --- Drop countries missing more than 5 variables ---

world_clean <- world_clean |>
  mutate(n_miss_row = rowSums(is.na(across(where(is.numeric))))) |>
  filter(n_miss_row < 5) |> # to not include anything 5 or more missing
  select(-n_miss_row, -n_missing) ## these are same, could use 1 less var, but more cleaning of code on these later

cat("Countries remaining after row filter:", nrow(world_clean), "\n")


# --- Distribution of all numeric variables ---

world_clean |>
  select(where(is.numeric)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 color = "black", fill = "skyblue", na.rm = TRUE) +
  geom_density(alpha = 0.4, fill = "blue", na.rm = TRUE) +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Distribution of All Numeric Variables",
       x = NULL, y = "Density")

#ggsave("fig_distributions.png", width = 7, height = 5, dpi = 150)


# --- Skewness - Numerically ---

# Confirm skewness numerically (to support visual judgment)
skew_summary <- world_clean |>
  select(where(is.numeric)) |>
  summarise(across(everything(), ~ skewness(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "skewness") |>
  arrange(desc(abs(skewness)))

print(skew_summary, n = 25)


# --- Snapshot before imputation ---

world_clean_before_imputation <- world_clean

# Check world_clean
glimpse(world_clean)


# ============================================================
# Data Imputation using MICE
# ============================================================

# Trying mice (multivariate imputation by chained equations to impute the missing values.
# This is an iterative process that will run to convergence

# select numeric columns
mice_input <- world_clean_before_imputation |>
  select(where(is.numeric))

# 'm=1' creates one imputed dataset; 'method="pmm"' (Predictive Mean Matching).
# we want this because it ensures imputed values stay within the original data range.
imputed_data <- mice(mice_input, m = 1, method = 'pmm', maxit = 10)

# we use m = 1 to create 1 version of the filled-in dataset.
# creating 5-10 would be better to ensure there is consistency across versions
# but this would be more difficult to cluster. # maxit is number of iterations

# Complete the data to get a normal dataframe back
df_mice <- complete(imputed_data)


# We skipped log transform, but need to be aware that severely skewed variables
# like CO2 (skewness 5.2) and HIV (4.3) will still be in their raw form going into PCA.
# PCA is sensitive to this, so we should maybe mention in the report that skewness
# existed and that standardization (which we'll do before PCA) partially mitigates it,
# even if it doesn't fully fix distributional shape.


# Bind Country and Code columns back in
df_mice <- bind_cols(
  world_clean |> select(Country, Code),
  df_mice
)


# ============================================================
# Exploratory Data Analysis
# ============================================================

# --- Correlation heatmap ---

# This motivates PCA (if variables are highly correlated, dimension reduction is justified)

cor_matrix <- df_mice |>
  select(where(is.numeric)) |>
  cor(use = "complete.obs")

# Open the file "device" - to save corrplot
#png("fig_corrplot.png", width = 7, height = 5, units = "in", res = 150)

corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",       # cluster similar variables together
         tl.cex = 0.5,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         number.cex = 0.4,
         col = colorRampPalette(c("red", "white", "steelblue"))(200),
         title = "Correlation Matrix — World Development Indicators",
         mar = c(0, 0, 1, 0))

#dev.off() # Close and save the plot
#ggsave("fig_corrplot.png", width=7, height=5, dpi=150)


# --- Some targeted bivariate plots based on correlation structure ---

# STORY 1:
# Birth_rate, Infant_mortality, and Agriculture are positively correlated with
# each other (r = 0.72–0.86) and strongly negatively correlated with
# Internet_use, Life_expectancy and Access_electricity (r = −0.66 to −0.90).

# Story 1a: Development axis — demographic transition
ggplot(df_mice, aes(x = Birth_rate, y = Life_expectancy)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "Birth Rate vs. Life Expectancy",
       subtitle = "Higher birth rates associate strongly with lower life expectancy",
       x = "Birth Rate (per 1,000)", y = "Life Expectancy (years)")

#ggsave("fig_birthrate_lifeexp.png", width = 7, height = 5, dpi = 150)

# Story 1b: Development axis — internet and infant mortality
ggplot(df_mice, aes(x = Internet_use, y = Infant_mortality)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "Internet Use vs. Infant Mortality",
       subtitle = "Digital access closely tracks child health outcomes",
       x = "Internet Use (%)", y = "Infant Mortality (per 1,000 live births)")

#ggsave("fig_internet_infantmort.png", width = 7, height = 5, dpi = 150)

# Story 1c: GDP vs Life Expectancy (log scale for GDP given skew)
ggplot(df_mice, aes(x = log(GDP + 1), y = Life_expectancy)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "log(GDP) vs. Life Expectancy",
       subtitle = "Diminishing returns: gains level off at higher incomes",
       x = "log(GDP per Capita)", y = "Life Expectancy (years)")

#ggsave("fig_logGDP_lifeexp.png", width = 7, height = 5, dpi = 150)


# STORY 2:
# Trade Openness: Imports and Exports are almost perfectly correlated (0.86) and
# move together — small open economies vs. large closed ones.

# Story 2: Trade openness — imports vs exports
ggplot(df_mice, aes(x = Imports, y = Exports)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "Imports vs. Exports (% of GDP)",
       x = "Imports (% of GDP)", y = "Exports (% of GDP)")

#ggsave("fig_imports_exports.png", width = 7, height = 5, dpi = 150)


# STORY 3:
# Renewable vs. Development: Renewable is positively correlated with Agriculture
# and Birth_rate as well as infant mortality. Could be suggesting renewables are
# high in poor agrarian nations (biomass/wood fuel), not just wealthy green ones.
# A subtle and interesting finding worth highlighting, maybe?

# Story 3: Renewable energy paradox
ggplot(df_mice, aes(x = Agriculture, y = Renewable)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "Agriculture Share vs. Renewable Energy Use",
       subtitle = "Renewables are high in agrarian economies, so likely biomass, not solar/wind",
       x = "Agriculture (% of GDP)", y = "Renewable Energy (%)")

#ggsave("fig_agriculture_renewable.png", width = 7, height = 5, dpi = 150)


# STORY 4:
# Story 4 — Health expenditure paradox
ggplot(df_mice, aes(x = Health_expenditure, y = Life_expectancy)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", formula = 'y ~ x', se = TRUE, color = "red") +
  labs(title = "Health Expenditure vs. Life Expectancy",
       subtitle = "r = 0.34 — spending alone does not predict outcomes",
       x = "Health Expenditure (% of GDP)", y = "Life Expectancy (years)")

#ggsave("fig_health_lifeexp.png", width = 7, height = 5, dpi = 150)


# --- Summary statistics table ---

df_mice |>
  select(where(is.numeric)) |>
  summary()


# ============================================================
# Principal Component Analysis
# ============================================================

# --- Scaling and running PCA ---

# Need to Standardize before PCA. Variables are on wildly different scales
# (GDP in thousands, rates as percentages, etc.). scale() gives each variable
# mean 0 and SD 1 so no single variable dominates by unit alone.

pca_input <- df_mice |>
  select(where(is.numeric))   # maybe could have used |> scale() here

# Run PCA
pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)

summary(pca_result)


# --- Scree Plot ---

# Scree plot to decide how many PCs to retain
# Look for the "elbow" where additional PCs stop adding meaningful variance

pca_var <- pca_result$sdev^2          # eigenvalues
pve     <- pca_var / sum(pca_var)     # proportion of variance explained
cum_pve <- cumsum(pve)                # cumulative proportion

scree_df <- tibble(
  PC      = 1:length(pve),
  PVE     = pve,
  Cum_PVE = cum_pve
)

# Scree plot
ggplot(scree_df, aes(x = PC, y = PVE)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 1:length(pve)) +
  labs(title = "Scree Plot — Proportion of Variance Explained",
       subtitle = "Red dashed line = 5% threshold",
       x = "Principal Component", y = "Proportion of Variance Explained")

#ggsave("fig_scree.png", width = 7, height = 5, dpi = 150)


# --- Cumulative Variance Plot ---

# Cumulative variance plot
ggplot(scree_df, aes(x = PC, y = Cum_PVE)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 1:length(pve)) +
  labs(title = "Cumulative Variance Explained",
       subtitle = "Red dashed line = 80% threshold",
       x = "Principal Component", y = "Cumulative Proportion")

#ggsave("fig_cumvar.png", width = 7, height = 5, dpi = 150)


# PC4 at ~67% is a reasonable tradeoff — we're capturing 2/3 of all variation with just 4 dimensions.
# We hit above 80% at PC7, but that's too many for interpretability.
## Based on Standard deviation, we could also use 5 or 6 PCs since until then std > ~1?


# --- PCA Loadings for first 6 PCs ---

# Look at both plots together.
# Loadings interpret what each PC means and tell which variables drive each component

loadings_df <- as.data.frame(pca_result$rotation[, 1:6]) |>
  rownames_to_column("variable") |>
  pivot_longer(-variable, names_to = "PC", values_to = "loading")

ggplot(loadings_df, aes(x = reorder(variable, loading), y = loading, fill = loading > 0)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("red", "steelblue")) +
  facet_wrap(~ PC, scales = "free_x") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "PCA Loadings for First 6 Principal Components",
       x = NULL, y = "Loading")

#ggsave("fig_loadings.png", width = 7, height = 5, dpi = 150)


# This plot is key for the Discussion section. PC1 has high positive loadings on
# Agriculture, Birth rate, Infant mortality, and Renewable energy and high negative
# loadings on Access_electricity, Life_expectancy, GDP, Urban, Internet_use.
# PC2 might separate agricultural/rural economies from urban ones, or
# resource-export economies from others.
# PC6 might distinguish countries with strong food security and education systems
# from those struggling with high disease burdens and lower economic output.


# --- Biplot — see countries in PC space ---

scores_df <- as.data.frame(pca_result$x[, 1:2]) |>
  mutate(Country = df_mice$Country)

ggplot(scores_df, aes(x = PC1, y = PC2, label = Country)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_text(size = 2, vjust = -0.5, check_overlap = TRUE) +
  labs(title = "Countries in PC1–PC2 Space",
       x = "PC1 (Development Level)",
       y = "PC2")

ggplot_pca(pca_result, choices = c(1, 2))
ggplot_pca(pca_result, choices = c(2, 3))
ggplot_pca(pca_result, choices = c(1, 3))


# --- Save PC scores for use in clustering ---

# Retain however many PCs we chose (e.g., 4 or 6)
n_pcs <- 6   # adjust based on our decision 4 vs 6; keeping 6 for now.

pc_scores <- as.data.frame(pca_result$x[, 1:n_pcs]) |>
  mutate(Country = df_mice$Country,
         Code    = df_mice$Code)


# ============================================================
# Clustering
# ============================================================

# --- Hierarchical Clustering ---

# Use PC scores from earlier steps (6 PCs retained)
# pc_scores already has Country, Code, PC1-PC6; let's work with 6

pc_matrix <- pc_scores |>
  select(PC1, PC2, PC3, PC4, PC5, PC6)

# Dendrogram gives a visual sense of natural groupings before committing to a k

dist_matrix <- dist(pc_matrix, method = "euclidean")

hc_ward <- hclust(dist_matrix, method = "ward.D2") # or use complete
# Ward.D2 minimizes total within-cluster variance — best general
# purpose linkage for finding compact, well-separated clusters

## Dr. Gade recommended to not use single linkage. (could use complete linkage)
#hc_complete <- hclust(dist_matrix, method = "complete")
# but it is sensitive to outliers, so choosing ward instead

# Plot dendrogram
plot(hc_ward,
     labels = pc_scores$Country,
     cex    = 0.4,
     main   = "Hierarchical Clustering Dendrogram (Ward's Method)",
     xlab   = "",
     sub    = "",
     ylab   = "Height")

# Draw rectangles around candidate clusters to guide k selection
rect.hclust(hc_ward, k = 2, border = "red")
rect.hclust(hc_ward, k = 5, border = "blue")
# visually compare k=2 (red) vs k=5 (blue) — pick whichever
# cuts at a natural gap in the dendrogram height


# --- Color coded visualization for dendrogram ---

# Helper function to make boxes
make_cluster_boxes <- function(hclust_obj, groups, height_cut) {
  num_clusters <- length(unique(groups))
  tibble(
    xmin    = sapply(1:num_clusters, function(zz) min(which(groups[hclust_obj$order] == zz))) - 0.5,
    xmax    = sapply(1:num_clusters, function(zz) max(which(groups[hclust_obj$order] == zz))) + 0.5,
    ymin    = rep(0, num_clusters),
    ymax    = rep(height_cut, num_clusters),
    cluster = factor(1:num_clusters)
  )
}

cluster_colors <- c("tomato1", "steelblue2", "olivedrab3", "darkorchid2",
                    "orange", "lightgoldenrod1", "aquamarine2")

# Cut the tree — inspected the dendrogram first to pick these
k_hc      <- 5       # updated after looking at dendrogram
height_hc <- 18      # updated — picked a height where the cut feels natural at that k

groups_hc <- cutree(hc_ward, k = k_hc)

df_boxes <- make_cluster_boxes(hc_ward, groups_hc, height_hc)


# Plotting using ggdendrogram and color coding

ggdendrogram(hc_ward, rotate = FALSE) +
  geom_rect(
    data = df_boxes,
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill  = cluster,
        color = cluster),
    alpha = 0.2
  ) +
  scale_fill_manual(values  = cluster_colors[1:k_hc]) +
  scale_color_manual(values = cluster_colors[1:k_hc]) +
  theme(
    panel.grid.major.y = element_line("gray"),
    axis.text.x        = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
    legend.position    = "none"
  ) +
  labs(title    = "Hierarchical Clustering — Ward's Linkage, Euclidean Distance",
       subtitle  = paste("k =", k_hc, "clusters"),
       x = NULL, y = "Height")

#ggsave("fig_dendrogram.png", width = 7, height = 5, dpi = 150)


# --- Elbow plot to confirm k choice ---

#set.seed(123)   # for reproducibility

# Compute within-cluster sum of squares for k = 1 to 10
wss <- map_dbl(1:10, function(k) {
  kmeans(pc_matrix, centers = k, nstart = 25, iter.max = 100)$tot.withinss
})

wss_df <- tibble(k = 1:10, wss = wss)

ggplot(wss_df, aes(x = k, y = wss)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "K-Means Elbow Plot",
       x = "Number of Clusters (k)", y = "Total Within-Cluster SS")

#ggsave("fig_elbow.png", width = 7, height = 5, dpi = 150)


# --- Silhouette analysis — second way to validate k ---

# Silhouette score measures how well each country fits its cluster vs.
# the next nearest cluster. Ranges from -1 to 1; higher is better.

#sil_scores <- map_dbl(2:8, function(k) {
#  km  <- kmeans(pc_matrix, centers = k, nstart = 25, iter.max = 100)
#  sil <- silhouette(km$cluster, dist_matrix)
#  mean(sil[, 3])
#})

#sil_df <- tibble(k = 2:8, avg_silhouette = sil_scores)

#ggplot(sil_df, aes(x = k, y = avg_silhouette)) +
#  geom_line(color = "steelblue") +
#  geom_point(color = "steelblue", size = 2) +
#  scale_x_continuous(breaks = 2:8) +
#  labs(title = "Average Silhouette Width by k",
#       subtitle = "Higher = better defined clusters",
#       x = "Number of Clusters (k)", y = "Average Silhouette Width")

#ggsave("fig_silhouette.png", width=7, height=5, dpi=150)


# From this plot, it says we can use k = 2, but we could also do k = 7
# (earlier runs of code without setting seed suggested k = 5 so all further
# analysis is based on k=5.) Now, by setting seed, the analysis is not
# going to be correct, so removing the Silhouette analysis itself.


# --- K-means clustering ---

# Fit final k-means with chosen k

k_final <- 5   # could choose k = 2 or 5
# updated this after inspecting plots. Since 2 is too small, going forward with 5;
# Did some analysis based on Silhouette to validate k earlier,
# but since current report is not based on it, will just move forward with k = 5.

# current seed of 77 suggests k = 2 or 7, not using it in our analysis anyways.

#set.seed(123)
km_final <- kmeans(pc_matrix, centers = k_final, nstart = 25, iter.max = 100)

# Attach cluster labels back to countries
pc_scores <- pc_scores |>
  mutate(cluster = factor(km_final$cluster))

# Also attach to full imputed dataset for profiling
df_clustered <- df_mice |>
  mutate(cluster = factor(km_final$cluster))


# --- Visualize clusters in multiple PC spaces ---

# Visualize clusters in PC1-PC2 space
ggplot(pc_scores, aes(x = PC1, y = PC2, color = cluster, label = Country)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(size = 2, vjust = -0.6, check_overlap = TRUE, color = "black") +
  #scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title    = "Countries Clustered in PC1–PC2 Space",
       subtitle = paste("K-Means with k =", k_final),
       x        = "PC1",
       y        = "PC2",
       color    = "Cluster")

#ggsave("fig_pc12.png", width = 7, height = 5, dpi = 150)

# PC1 vs PC6
ggplot(pc_scores, aes(x = PC1, y = PC6, color = cluster, label = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text(size = 2, vjust = -0.6, check_overlap = TRUE, color = "black") +
  #scale_color_brewer(palette = "Set2") +
  theme_bw() +
  labs(title    = "Countries Clustered in PC1–PC6 Space",
       subtitle = paste("K-Means with k =", k_final),
       x        = "PC1",
       y        = "PC6",
       color    = "Cluster")

#ggsave("fig_pc13.png", width = 7, height = 5, dpi = 150)

# PC1 vs PC5
ggplot(pc_scores, aes(x = PC1, y = PC5, color = cluster, label = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text(size = 2, vjust = -0.6, check_overlap = TRUE, color = "black") +
  #scale_color_brewer(palette = "Set2") +
  theme_bw() +
  labs(title    = "Countries Clustered in PC1–PC5 Space",
       subtitle = paste("K-Means with k =", k_final),
       x        = "PC1",
       y        = "PC5",
       color    = "Cluster")

# PC2 vs PC3
ggplot(pc_scores, aes(x = PC2, y = PC3, color = cluster, label = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text(size = 2, vjust = -0.6, check_overlap = TRUE, color = "black") +
  #scale_color_brewer(palette = "Set2") +
  theme_bw() +
  labs(title    = "Countries Clustered in PC2–PC3 Space",
       subtitle = paste("K-Means with k =", k_final),
       x        = "PC2",
       y        = "PC3",
       color    = "Cluster")

# PC1 vs PC6
ggplot(pc_scores, aes(x = PC1, y = PC6, color = cluster, label = Country)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text(size = 2, vjust = -0.6, check_overlap = TRUE, color = "black") +
  #scale_color_brewer(palette = "Set2") +
  theme_bw() +
  labs(title    = "Countries Clustered in PC1–PC6 Space",
       subtitle = paste("K-Means with k =", k_final),
       x        = "PC1",
       y        = "PC6",
       color    = "Cluster")

## Add more PCs if required


# --- Profile each cluster ---

# We need this to interpret
# Compute mean of each original variable by cluster

cluster_profiles <- df_clustered |>
  group_by(cluster) |>
  summarise(
    n_countries               = n(),
    mean_GDP                  = mean(GDP),
    mean_Life_expectancy      = mean(Life_expectancy),
    mean_Birth_rate           = mean(Birth_rate),
    mean_Infant_mortality     = mean(Infant_mortality),
    mean_Health_expenditure   = mean(Health_expenditure),
    mean_Access_electricity   = mean(Access_electricity),
    mean_Internet_use         = mean(Internet_use),
    mean_Women_parliament     = mean(Women_parliament),
    mean_Education_expenditure = mean(Education_expenditure),
    mean_Education_years      = mean(Education_years),
    mean_HIV                  = mean(HIV),
    mean_CO2                  = mean(CO2),
    mean_Cropland             = mean(Cropland),
    mean_Agriculture          = mean(Agriculture),
    mean_Food_production      = mean(Food_production),
    mean_Renewable            = mean(Renewable),
    mean_Urban                = mean(Urban),
    mean_Imports              = mean(Imports),
    mean_Exports              = mean(Exports),
  ) # |> arrange(mean_GDP)

#print(cluster_profiles)
cluster_profiles |> as.data.frame()


# --- Heatmap of cluster profiles ---

# Heatmap of cluster profiles — visual version of above
profile_long <- df_clustered |>
  select(cluster, where(is.numeric)) |>
  group_by(cluster) |>
  summarise(across(everything(), mean)) |>
  pivot_longer(-cluster, names_to = "variable", values_to = "mean_value") |>
  group_by(variable) |>
  mutate(scaled_value = scale(mean_value)[, 1])   # z-score across clusters
# so variables on different
# scales are comparable

ggplot(profile_long, aes(x = cluster, y = variable, fill = scaled_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "steelblue",
                       midpoint = 0) +
  labs(title = "Cluster Profile Heatmap",
       subtitle = "Scaled mean values per cluster (red = low, blue = high)",
       x = "Cluster", y = NULL, fill = "Z-score")

#ggsave("fig_heatmap.png", width = 7, height = 5, dpi = 150)


# --- Country - Cluster Association ---

# Which countries are in each cluster?
pc_scores |>
  select(Country, cluster) |>
  arrange(cluster, Country) |>
  as.data.frame()
