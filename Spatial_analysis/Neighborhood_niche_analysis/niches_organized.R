library(dplyr)
library(purrr)
library(RANN)
library(Seurat)
library(ggplot2)

sets_epi_major <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_major.rds")
sets_epi_sub <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_sub.rds")
sets_others <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_others.rds")

xenium_data_list <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/basic_analysis/Surface/Surface_xenium_data_list_cutoff20.rds")
comb_meta <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/comb_meta_annot.rds")
table(comb_meta$slide_id)
sum(is.na(comb_meta$x))

metadata_list <- list()
for (sec in c(1:6)){
  xenium_sec01 <- xenium_data_list[[sec]]
  xenium_sec01@meta.data <- xenium_sec01@meta.data %>% left_join(comb_meta[, c(11, 12, 21:39)])
  rownames(xenium_sec01@meta.data) <- xenium_sec01@meta.data$barcode
  
  xenium_sec01@meta.data$x <- xenium_sec01@images$fov$centroids@coords[, 1]
  xenium_sec01@meta.data$y <- xenium_sec01@images$fov$centroids@coords[, 2]
  sum(is.na(xenium_sec01@meta.data$x))
  
  metadata_list[[sec]] <- xenium_sec01@meta.data
}

metadata <- do.call(rbind, metadata_list)
table(metadata$slide_id)
sum(is.na(metadata$x))

colnames(metadata)

# Step 1: Function to compute neighbor class counts for one section
compute_neighbor_matrix <- function(sec_df, k = 30) {
  coords <- sec_df[, c("x", "y")]
  nn <- RANN::nn2(coords, k = k + 1)  # +1 because the cell itself is included
  
  # Create a consistent final_annot set across sections
  all_classes <- sort(unique(metadata$epi_annot_detail))
  neighbor_matrix <- matrix(0, nrow = nrow(sec_df), ncol = length(all_classes))
  colnames(neighbor_matrix) <- all_classes
  rownames(neighbor_matrix) <- sec_df$barcode
  
  for (i in seq_len(nrow(sec_df))) {
    neighbors <- nn$nn.idx[i, ]
    neighbor_classes <- sec_df$epi_annot_detail[neighbors]
    neighbor_counts <- table(factor(neighbor_classes, levels = all_classes))
    neighbor_matrix[i, ] <- as.numeric(neighbor_counts)
  }
  
  return(as.data.frame(neighbor_matrix) %>%
           mutate(barcode = sec_df$barcode,
                  slide_id = sec_df$slide_id,
                  celltype = sec_df$celltype,
                  epi_annot = sec_df$epi_annot,
                  epi_annot_detail = sec_df$epi_annot_detail))
}

# Step 2: Run for each section, combine results
neighborhood_df <- metadata %>%
  group_split(slide_id) %>%
  map_df(compute_neighbor_matrix)

# Step 3: Prepare for Seurat
neighbor_matrix <- neighborhood_df %>%
  select(-barcode, -slide_id, -celltype, -epi_annot, -epi_annot_detail) %>%
  as.matrix()
rownames(neighbor_matrix) <- neighborhood_df$barcode

meta_data <- neighborhood_df %>%
  select(barcode, slide_id, celltype, epi_annot, epi_annot_detail) 
rownames(meta_data) <- meta_data$barcode


#
scaled_mat <- scale(neighbor_matrix)
set.seed(123)


km <- kmeans(scaled_mat, centers = 5)
clusters <- km$cluster

str(clusters)

### find k for k-means
library(mclust)

# find k
cluster_list <- list()
i <- 0
for (k in c(3:20)){
  i <- i+1
  km <- kmeans(scaled_mat, centers = k)
  cluster_list[[i]] <- km$cluster
}

ARI_list <- as.numeric()
for (i in c(1:17)){
  ARI_list[i] <- adjustedRandIndex(cluster_list[[i]], cluster_list[[i+1]])
}

plot(x = c(3:19), y = ARI_list, type = "l")

index_list <- as.numeric()
for (i in c(1:16)){
  index_list[i] <- (ARI_list[i]+ARI_list[i+1])/2
}

plot(x = c(4:19), y = index_list, type = "l") # k=8 if k=20 for neighbor; k=10 if k=10, 15 for neighbor; k=10 if k=30 for neighbor

###
clusters <- cluster_list[[13]]
meta_data$cluster_k15 <- clusters
clusters <- cluster_list[[8]]
meta_data$cluster_k10 <- clusters


write.csv(meta_data, "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_meta.csv")
saveRDS(cluster_list, "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_list.rds")

#table(meta_data$cluster_k10, meta_data$celltype)
#table(meta_data$cluster_k10, meta_data$epi_annot)
#table(meta_data$cluster_k10, meta_data$epi_annot_detail)
table(meta_data$epi_annot_detail, meta_data$cluster_k15)
table(meta_data$epi_annot_detail, meta_data$cluster_k10)

###
meta_data$niche <- factor(paste0("Niche ", meta_data$cluster_k10), levels = paste0("Niche ", c(1:10)))

celltype_df <- meta_data %>%
  count(niche) %>%
  mutate(proportion = n / sum(n))

# Plot
library(RColorBrewer)

p <- ggplot(celltype_df, aes(x = niche, y = proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),
            hjust = -0.1,  # adjust position (works with coord_flip)
            size = 4) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Niches", y = "Proportion") +
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip() +
  ylim(0, 0.37) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Niche Proportion")
  

pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_proportion_k10.pdf", width = 6, height = 4)
print(p)
dev.off()

## neighborhood of the niches
table(meta_data$cluster_k10, meta_data$epi_annot_detail)

# Average neighborhood composition per niche
niche_profile <- neighbor_matrix %>%
  as.data.frame() %>%
  mutate(niche = meta_data$niche) %>%
  group_by(niche) %>%
  summarise(across(everything(), mean))

niche_name <- niche_profile$niche
niche_profile <- niche_profile[, -1]
niche_profile <- apply(niche_profile, 2, as.numeric)
rownames(niche_profile) <- niche_name

# Heatmap of niche compositions
library(pheatmap)
library(RColorBrewer)
library(viridis)
pheatmap(as.matrix(niche_profile),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = data.frame(niche = rownames(niche_profile)),
         scale = "none")
pheatmap(niche_profile, scale = "column", cluster_cols = F, cluster_rows = F)

pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_composition_heatmap.pdf", width = 5, height = 4)
p <- pheatmap(niche_profile, scale = "column", cluster_cols = F, cluster_rows = F, color = viridis(n = 10, alpha = 1, 
                                                                                                   begin = 0, end = 1, option = "viridis"), border_color = "black")
print(p)
dev.off()

#========
# plot niches
niche_label <- paste0("Niche ", c(1:10))
niche_color <- brewer.pal(12, "Paired")[c(1:10)]
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A"
sets_niche <- setNames(niche_color, niche_label)

for (sec in c(1:6)){
  xenium_sec01 <- xenium_data_list[[sec]]
  
  xenium_sec01@meta.data <- xenium_sec01@meta.data %>% left_join(meta_data)
  rownames(xenium_sec01@meta.data) <- xenium_sec01@meta.data$barcode
  
  xenium_sec01@meta.data <- xenium_sec01@meta.data %>% left_join(comb_meta[, c(11, 12, 21:39)])
  rownames(xenium_sec01@meta.data) <- xenium_sec01@meta.data$barcode
  
  xenium_sec01@meta.data$x <- xenium_sec01@images$fov$centroids@coords[, 1]
  xenium_sec01@meta.data$y <- xenium_sec01@images$fov$centroids@coords[, 2]
  
  colnames(xenium_sec01@meta.data)
  
  DefaultBoundary(xenium_sec01[["fov"]]) <- "segmentation"
  
  pdf(paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_sec0", sec, ".pdf"), width = 20, height = 16)
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "niche", border.size = 0.05, cols = sets_niche))
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "niche", border.size = 0.05, cols = sets_niche, dark.background = F, border.color = "grey"))
  dev.off()
  
  saveRDS(xenium_sec01, paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/sec0", sec, "_object.rds"))

}

#=======
library(tidyverse)
sets_epi_major <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_major.rds")
sets_epi_sub <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_sub.rds")
sets_others <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_others.rds")


df <- as.data.frame(niche_profile)
df$Niche <- rownames(df)

df_long <- df %>%
  pivot_longer(
    cols = -Niche,
    names_to = "CellType",
    values_to = "Fraction"
  )

df_long$Niche <- factor(df_long$Niche, levels = rownames(niche_profile))

p <- ggplot(df_long, aes(x = Niche, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(sets_epi_sub, sets_others)) +
  theme_bw() +
  xlab("Niche") +
  ylab("Cell composition") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_composition_stacked_bar.pdf", height = 5, width = 7)
print(p)
dev.off()

p <- ggplot(df_long, aes(x = Niche, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(sets_epi_sub, sets_others)) +
  coord_flip() +
  theme_bw() +
  xlab("Niche") +
  ylab("Cell composition") +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(size = 10)
  )

p <- p + guides(fill = guide_legend(nrow = 8, byrow = TRUE))

pdf("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/niche_composition_stacked_bar_rotate.pdf", height = 9, width = 9)
print(p)
dev.off()
