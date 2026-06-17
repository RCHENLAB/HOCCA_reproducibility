library(Seurat)
library(dplyr)
library(ggplot2)
library(princurve)

sets_epi_major <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_major.rds")
sets_epi_sub <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_sub.rds")
sets_others <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_others.rds")

one_donor <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/one_donor.rds")

slide_list <- paste0("sec0", c(1:6))
position_list <- list()
pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/epi_cell_proportion_across_bins.pdf", width = 6, height = 4)
for (i in c(1:6)){
  df <- one_donor@meta.data[one_donor@meta.data$slide_id == slide_list[i],]
  
  xy <- as.matrix(df[, c("x", "y")])
  
  set.seed(1)
  fit <- principal_curve(xy)  # fits a smooth curve through the cloud
  
  # 1D projection (arc-length-like parameter along the curve)
  lambda <- fit$lambda              # numeric vector, one value per point
  lambda_scaled <- (lambda - min(lambda)) / diff(range(lambda))  # 0..1
  
  # add back to your data
  if (i %in% c(2, 6)){
    df$pcurve_lambda <- lambda
    df$pcurve_pos01 <- 1 - lambda_scaled
  } else {
    df$pcurve_lambda <- lambda
    df$pcurve_pos01 <- lambda_scaled
  }
  
  
  # quick visualization
  p <- ggplot(df, aes(x, y, color = pcurve_pos01)) +
    geom_point(size = 1) +
   # geom_path(data = as.data.frame(fit$s), aes(x = x, y = y),
  #            color = "black", linewidth = 1) +
    scale_color_viridis_c(name = "position (0..1)") +
    coord_equal() + theme_minimal()
  print(p)
  

  # Example: df with pcurve_pos01 ∈ [0,1] and annot categories
  df_binned <- df %>%
    mutate(
      bin = cut(pcurve_pos01, breaks = 30, include.lowest = TRUE)  # 50 bins
    ) %>%
    group_by(bin, epi_annot_detail) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(bin) %>%
    mutate(proportion = n / sum(n))
  
  # Inspect results
  p <- ggplot(df_binned, aes(x = factor(bin), y = proportion, fill = epi_annot_detail)) +
    geom_bar(stat = "identity", width = 0.9) +
    theme_classic() +
    labs(x = "Cornea -> Limbus", y = "Proportion", fill = "", title = slide_list[i]) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = sets_epi_sub)
  
  print(p)
  
  
  # order points along the curve
  #df_ordered <- df[order(df$pcurve_lambda), ] 
  position_list[[i]] <- df
}

position_df <- do.call(rbind, position_list)
df_binned <- position_df %>%
  mutate(
    bin = cut(pcurve_pos01, breaks = 30, include.lowest = TRUE)  # 50 bins
  ) %>%
  group_by(bin, epi_annot_detail) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin) %>%
  mutate(proportion = n / sum(n))

p <- ggplot(df_binned, aes(x = factor(bin), y = proportion, fill = epi_annot_detail)) +
  geom_bar(stat = "identity", width = 0.9) +
  theme_classic() +
  labs(x = "Cornea -> Limbus", y = "Proportion", fill = "", title = "Epithelial Cell Proportion") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = sets_epi_sub)

print(p)

dev.off()
