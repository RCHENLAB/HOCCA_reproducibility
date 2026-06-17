library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)


nhood <- read.csv("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/flat_df.csv")

heatmap_data <- nhood %>%
  group_by(row, col) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')
ggplot(heatmap_data, aes(x = factor(col), y = factor(row), fill = mean_value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean Value") +
  theme_minimal() +
  labs(x = "Column", y = "Row", title = "Heatmap of Mean Values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

wide_matrix <- heatmap_data %>%
  mutate(row = as.character(row), col = as.character(col)) %>%
  pivot_wider(names_from = col, values_from = mean_value) %>%
  arrange(row)
rownames_matrix <- wide_matrix$row
wide_matrix_mat <- as.matrix(wide_matrix[ , -1])
rownames(wide_matrix_mat) <- rownames_matrix
pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/nhood_enrichment_matrix/nhood_enrichment_no_cluster.pdf", width = 7, height = 6.5)
p <- pheatmap(wide_matrix_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE)
print(p)
dev.off()
