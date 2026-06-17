library(Seurat)
library(dplyr)

xenium_data_list <- readRDS(paste0(output_PATH, target_tissue, "_xenium_data_list_cutoff20.rds"))

xenium_sec01 <- xenium_data_list[[1]]
ImageDimPlot(xenium_sec01, fov = "fov", molecules = c("KRT5"), nmols = 20000)
ImageFeaturePlot(xenium_sec01, features = c("KRT5"), max.cutoff = c(25), size = 0.75, cols = c("white", "red"))

xenium_sec01@meta.data <- xenium_sec01@meta.data %>% left_join(xenium.combined@meta.data)
rownames(xenium_sec01@meta.data) <- xenium_sec01@meta.data$barcode
ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_celltype", size = 1.5)
ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_majorclass")

ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_celltype", size = 1.5, dark.background = F)


library(fastDummies)
df_dummy <- dummy_cols(xenium_sec01@meta.data, select_columns = "xenium_clusters_p3_celltype", remove_first_dummy = FALSE)
rownames(df_dummy) <- df_dummy$barcode
colnames(df_dummy)
colnames(df_dummy) <- gsub(" ", "_", colnames(df_dummy))
colnames(df_dummy) <- gsub("/", "_", colnames(df_dummy))
colnames(df_dummy)
xenium_sec01@meta.data <- df_dummy
ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_celltype_Corneal_Stromal_Fibroblasts", size = 1.5, dark.background = F)
DefaultBoundary(xenium_sec01[["fov"]]) <- "segmentation"

ir_list <- colnames(xenium_sec01@meta.data)[20:32]
pdf("basic_analysis/Surface/annot_sep_sec01.pdf", width = 20, height = 16)
for (i in ir_list){
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = i, size = 1.5, dark.background = T, border.size = 0.1, border.color = "white") + 
          scale_fill_manual(values = c("0" = "grey", "1" = "red")))
}
dev.off()

pdf("basic_analysis/Surface/annot_sec01.pdf", width = 20, height = 16)
print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_celltype", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p5", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_1", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
dev.off()

for (sec in c(2:5)){
  xenium_sec01 <- xenium_data_list[[sec]]
  xenium_sec01@meta.data <- xenium_sec01@meta.data %>% left_join(xenium.combined@meta.data)
  rownames(xenium_sec01@meta.data) <- xenium_sec01@meta.data$barcode
  
  df_dummy <- dummy_cols(xenium_sec01@meta.data, select_columns = "xenium_clusters_p3_celltype", remove_first_dummy = FALSE)
  rownames(df_dummy) <- df_dummy$barcode
  colnames(df_dummy)
  colnames(df_dummy) <- gsub(" ", "_", colnames(df_dummy))
  colnames(df_dummy) <- gsub("/", "_", colnames(df_dummy))
  colnames(df_dummy)
  xenium_sec01@meta.data <- df_dummy
  DefaultBoundary(xenium_sec01[["fov"]]) <- "segmentation"
  
  ir_list <- colnames(xenium_sec01@meta.data)[20:32]
  pdf(paste0("basic_analysis/Surface/annot_sep_sec0", sec, ".pdf"), width = 20, height = 16)
  for (i in ir_list){
    print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = i, size = 1.5, dark.background = T, border.size = 0.1, border.color = "white") + 
            scale_fill_manual(values = c("0" = "grey", "1" = "red")))
  }
  dev.off()
  
  pdf(paste0("basic_analysis/Surface/annot_sec0", sec, ".pdf"), width = 20, height = 16)
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3_celltype", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p3", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_p5", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
  print(ImageDimPlot(xenium_sec01, fov = "fov", group.by = "xenium_clusters_1", size = 1.5, dark.background = T, border.size = 0.1, border.color = "white"))
  dev.off()
  
}

