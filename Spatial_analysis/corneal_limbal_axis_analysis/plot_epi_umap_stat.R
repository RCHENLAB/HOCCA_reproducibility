library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
source("/dfs3b/ruic20_lab/tingty7/codes/Xenium/subset_obj_seurat_v2.R")


sets_epi_major <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_major.rds")
sets_epi_sub <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_epi_sub.rds")
sets_others <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/colors/sets_others.rds")

outdir <- "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/"

xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/coembed/epi_coembed.combined_umap_annot.rds")
comb_meta <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/comb_meta_annot.rds")

colnames(xenium.combined@meta.data)
colnames(comb_meta)

xenium.combined@meta.data$comb_barcode <- rownames(xenium.combined@meta.data)
xenium.combined@meta.data <- xenium.combined@meta.data %>% left_join(comb_meta[, c(11, 12, 37:39)])
rownames(xenium.combined@meta.data) <- xenium.combined@meta.data$comb_barcode

table(xenium.combined$celltype)
table(xenium.combined$epi_annot)
table(xenium.combined$epi_annot_detail)


epi_clean <- subset(xenium.combined, subset = celltype %in% c("Epithelium"))
dim(epi_clean)

pdf(paste0(outdir, "epi_coembed.combined_umap_annot_clean.pdf"), width = 8, height = 5)
DimPlot(epi_clean, group.by = "epi_annot", label = T, repel = T, cols = sets_epi_major) + ggtitle("Epithelial Cells")
DimPlot(epi_clean, group.by = "epi_annot_detail", label = T, repel = T, cols = sets_epi_sub) + ggtitle("Epithelial Cells - subtypes")
DimPlot(epi_clean, group.by = "epi_annot", label = F, repel = T, cols = sets_epi_major) + ggtitle("Epithelial Cells")
DimPlot(epi_clean, group.by = "epi_annot_detail", label = F, repel = T, cols = sets_epi_sub) + ggtitle("Epithelial Cells - subtypes")
dev.off()

pdf(paste0(outdir, "epi_coembed.combined_umap_annot_clean_split.pdf"), width = 40)
DimPlot(epi_clean, group.by = "epi_annot", split.by = "slide_id", cols = sets_epi_major) + ggtitle("Epithelial Cells")
DimPlot(epi_clean, group.by = "epi_annot_detail", split.by = "slide_id", cols = sets_epi_sub) + ggtitle("Epithelial Cells - subtypes")
dev.off()

table(xenium.combined$epi_annot_detail)
table(epi_clean$epi_annot_detail)

saveRDS(epi_clean, paste0(outdir, "epi_clean.rds"))

table(epi_clean$slide_id)
one_donor <- subset_opt(epi_clean, subset = slide_id %in% c(paste0("sec0", c(1:6))))
table(one_donor$epi_annot, one_donor$epi_annot_detail)
saveRDS(one_donor, paste0(outdir, "one_donor.rds"))

##### plot cell proportion

# Calculate proportions
majorclass_df <- one_donor@meta.data %>%
  count(epi_annot) %>%
  mutate(proportion = n / sum(n))

# Plot
p <- ggplot(majorclass_df, aes(x = epi_annot, y = proportion, fill = epi_annot)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),
            hjust = -0.1,  # adjust position (works with coord_flip)
            size = 4) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Epithelial Cells", y = "Proportion") +
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip() +
  ylim(0, 0.38) +
  scale_fill_manual(values = sets_epi_major) +
  ggtitle("Epithelial cell type proportion")

pdf(paste0(outdir, "epi_one_donor_proportion.pdf"), width = 6, height = 3)
print(p)
dev.off()

#
celltype_df <- one_donor@meta.data %>%
  count(epi_annot_detail) %>%
  mutate(proportion = n / sum(n))

# Plot
p <- ggplot(celltype_df, aes(x = epi_annot_detail, y = proportion, fill = epi_annot_detail)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)),
            hjust = -0.1,  # adjust position (works with coord_flip)
            size = 4) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Epithelial Cells", y = "Proportion") +
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip() +
  ylim(0, 0.17) +
  scale_fill_manual(values = sets_epi_sub) +
  ggtitle("Epithelial cell states proportion")

pdf(paste0(outdir, "epi_state_one_donor_proportion.pdf"), width = 6, height = 5)
print(p)
dev.off()

unique(celltype_df$epi_clusters_2_annot_detail)


#=======
Idents(epi_clean) <- "epi_annot"
markers <- FindAllMarkers(epi_clean, only.pos = TRUE, group.by = "epi_annot")
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/epi_subtype_DEGs.csv")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  dplyr::filter(pct.1 > 0.4) %>%
  ungroup() -> top5

p <- DoHeatmap(epi_clean, features = top5$gene, group.colors = sets_epi_major)
pdf(paste0(outdir, "epi_DEGs_heatmap.pdf"), width = 7, height = 5)
print(p)
dev.off()

#-======
Idents(epi_clean) <- "epi_annot_detail"
markers <- FindAllMarkers(epi_clean, only.pos = TRUE, group.by = "epi_annot_detail")
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/surface_combine/spatial_plots/epi_states_DEGs.csv")

FeaturePlot(epi_clean, features = c("VIT", "GJB6"), min.cutoff = "q10", max.cutoff = "q90")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  dplyr::filter(pct.1 > 0.4) %>%
  ungroup() -> top5

p <- DoHeatmap(epi_clean, features = top5$gene, group.colors = sets_epi_sub)
pdf(paste0(outdir, "epi_state_DEGs_heatmap.pdf"), width = 9, height = 9)
print(p)
dev.off()
