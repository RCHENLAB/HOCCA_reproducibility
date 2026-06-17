library(hypeR)
library(dplyr)
library(ggplot2)

epi_module <- read.csv("/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA_wo_fetal/RNA_velocity/hotspot/mod_reordered.csv")
table(epi_module$Module)

module_list <- list()
for (i in c(1:7)){
  module_list[[i]] <- epi_module$X[epi_module$Module == i]
}
names(module_list) <- paste0("Module", c(1:7))

GO_BP <- msigdb_gsets(species="Homo sapiens", collection="C5", subcollection = "BP")

mhyp <- hypeR(module_list, GO_BP, test = "hypergeometric", background = 36601, plotting = T, fdr = 0.05)

saveRDS(mhyp, "/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA_wo_fetal/RNA_velocity/module_GO_enrichment.rds")

p <- hyp_dots(mhyp)
p

hyp_dots(mhyp, merge=TRUE)

result_list <- list()
for (i in c(1:7)){
  result <-  mhyp[["data"]][[i]][["data"]]
  result$module <- paste0("Module", i)
  result_list[[i]] <- result
}
final_results <- do.call(rbind, result_list)
final_results$gene_ratio <- final_results$overlap/final_results$signature
write.csv(final_results, "/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA_wo_fetal/RNA_velocity/module_GO_result_full.csv")

sub_final_results <- final_results[final_results$gene_ratio > 0.08,]
table(sub_final_results$module)
#sub_final_results <- sub_final_results[sub_final_results$module != "Module7",]


for (i in c(1:7)){
  print(m <- paste0("Module", c(1:7))[i])
  print(head(sub_final_results$label[sub_final_results$module == m]))
}

sub_final_results %>%
  group_by(module) %>%
  slice_max(order_by = gene_ratio, n = 5) %>%  # top 5 per module
  ungroup() %>%
  mutate(label = factor(label, levels = rev(unique(label)))) %>%  # reverse for horizontal bars
  ggplot(aes(x = gene_ratio, y = label, fill = -log10(fdr))) +
  geom_col() +
  facet_wrap(~ module, scales = "free_y") +
  scale_fill_gradient(low = "lightblue", high = "red") +
  labs(x = "Gene Ratio", y = NULL, fill = "-log10(FDR)") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

sub_final_results <- sub_final_results %>%
  mutate(log10_fdr_capped = pmin(-log10(fdr), 10))

top_go <- sub_final_results %>%
  group_by(module) %>%
  slice_head(n = 5) %>%
  ungroup()

# Plot horizontal barplot, one facet per module
p <- ggplot(top_go, aes(x = gene_ratio, y = reorder(label, gene_ratio), fill = log10_fdr_capped)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) + 
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(FDR)") +
  labs(x = "Gene Ratio", y = "GO Term") +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_blank(),     # Hide y-axis text
    axis.ticks.y = element_blank(),    # Hide y-axis ticks
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  ) +
  xlim(0, max(top_go$gene_ratio) * 2)

pdf("/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA_wo_fetal/RNA_velocity/module_GO_barplot.pdf", width = 10, height = 15)
print(p)
dev.off()
