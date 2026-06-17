library(dplyr)
library(purrr)
library(RANN)
library(Seurat)
library(ggplot2)

sets_epi_major <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_major.rds")
sets_epi_sub <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_sub.rds")
sets_others <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_others.rds")

outdir <- "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/"

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

library(dplyr)
library(RANN)

cell_id_col <- "barcode"
slide_col <- "slide_id"
annot_col <- "epi_annot_detail"
x_col <- "x"
y_col <- "y"

find_closest_melanocyte <- function(df,
                                    cell_id_col = "cell_id",
                                    annot_col = "epi_annot_detail",
                                    x_col = "x",
                                    y_col = "y") {

    mel_df <- df %>% filter(.data[[annot_col]] == "Melanocytes")
  
  if (nrow(mel_df) == 0) {
    df$nearest_melanocyte_id <- NA
    df$nearest_melanocyte_dist <- NA
    return(df)
  }
  
  mel_mat <- as.matrix(mel_df[, c(x_col, y_col)])
  query_mat <- as.matrix(df[, c(x_col, y_col)])
  
  nn <- RANN::nn2(data = mel_mat, query = query_mat, k = 1)
  
  df$nearest_melanocyte_id <- mel_df[[cell_id_col]][nn$nn.idx[, 1]]
  df$nearest_melanocyte_dist <- nn$nn.dists[, 1]
  
  return(df)
}

comb_meta_out <- metadata %>%
  group_split(.data[[slide_col]]) %>%
  lapply(function(df) {
    find_closest_melanocyte(
      df,
      cell_id_col = cell_id_col,
      annot_col = annot_col,
      x_col = x_col,
      y_col = y_col
    )
  }) %>%
  bind_rows()

library(ggplot2)
library(dplyr)

df_plot <- comb_meta_out %>%
  filter(!is.na(nearest_melanocyte_dist),
         is.finite(nearest_melanocyte_dist))

p <- ggplot(df_plot, aes(x = epi_annot_detail, y = nearest_melanocyte_dist)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  xlab("") +
  ylab("Distance") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Distance to closest Melanocyte")

pdf(paste0(outdir, "boxplot_dist_to_closest_melanocytes.pdf"), width = 8, height = 4)
print(p)
dev.off()

df_plot <- df_plot[df_plot$celltype == "Epithelium",]

p <- ggplot(df_plot, aes(x = epi_annot_detail, y = nearest_melanocyte_dist)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  xlab("") +
  ylab("Distance") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Distance to closest Melanocyte")

pdf(paste0(outdir, "boxplot_dist_to_closest_melanocytes_epi.pdf"), width = 4, height = 4)
print(p)
dev.off()

mean_df <- df_plot %>%
  group_by(epi_annot_detail) %>%
  summarise(
    mean_dist = mean(nearest_melanocyte_dist),
    median_dist = median(nearest_melanocyte_dist),
    n = n(),
    .groups = "drop"
  )

mean_df


#======
df_plot <- df_plot[df_plot$epi_annot_detail %in% c("Limbus_basal_1", "Limbus_basal_2", "Limbus_suprabasal", "Limbus_superficial"),]
library(forcats)
sum_df <- df_plot %>%
  group_by(epi_annot_detail) %>%
  summarise(
    mean_dist = mean(nearest_melanocyte_dist),
    sem_dist = sd(nearest_melanocyte_dist) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p <- ggplot(sum_df, aes(x = epi_annot_detail, y = mean_dist, fill = epi_annot_detail)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_dist - sem_dist, ymax = mean_dist + sem_dist),
    width = 0.2
  ) +
  scale_fill_manual(values = sets_epi_sub) +
  theme_bw() +
  xlab("") +
  ylab("Distance") +
  ggtitle("Distance to closest Melanocyte") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

pdf(paste0(outdir, "boxplot_dist_to_closest_melanocytes_epi_limbus.pdf"), width = 4, height = 4)
print(p)
dev.off()
