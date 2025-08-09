# ─────────────────────────────────────────────────────────────────────────────
# Setup for Micheli dataset (fresh start)
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(crayon)
  library(stringr)
  library(harmony)
  library(patchwork)
  library(scico)
})

# Updated base directory
base_dir     <- "D:/Masters/Processing Codes/Micheli_Cosgrove_Annotation"
seurat_dir   <- file.path(base_dir, "Seurat_Object")
micheli_rds <- file.path(seurat_dir, "ref_DeMicheli.rds")

# Output folders specific to this version of Micheli
micheli_dir              <- base_dir  # Keep for consistency
umap_dir                  <- file.path(base_dir, "Plots", "UMAPs")
barplot_dir               <- file.path(base_dir, "Plots", "Barplots")
prop_dir                  <- file.path(base_dir, "Proportions")
log_dir                   <- file.path(base_dir, "Logs")
annotation_base_dir       <- file.path(base_dir, "Annotation")
renad_plot_dir            <- file.path(annotation_base_dir, "Renad Markers")
top_marker_dir            <- file.path(annotation_base_dir, "Top Markers")

# Ensure all directories are created
dir.create(seurat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(renad_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(top_marker_dir, recursive = TRUE, showWarnings = FALSE)

# Load object
cat(blue$bold("\n Loading Micheli Seurat object...\n"))
micheli <- readRDS(micheli_rds)
DefaultAssay(micheli) <- "RNA"

# ─────────────────────────────────────────────────────────────────────────────
# Run PCA + Harmony + fresh UMAP specific to this dataset
# ─────────────────────────────────────────────────────────────────────────────

cat(blue$bold(" Running PCA + Harmony + UMAP...\n"))

# Run standard scaling and PCA
micheli <- ScaleData(micheli, verbose = TRUE)
micheli <- RunPCA(micheli, npcs = 20, verbose = TRUE)

# Run Harmony batch correction using 'sample' metadata

micheli <- RunHarmony(
  object = micheli,
  group.by.vars = "sample",        
  reduction.use = "pca",
  assay.use = "RNA",
  verbose = TRUE
)

# Generate UMAP from Harmony-reduced space
micheli <- RunUMAP(
  object = micheli,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap_harmony_clean",
  reduction.key = "umapHarmonyClean_"
)

cat(green$bold("Harmony-corrected UMAP complete.\n"))

final_umap_rds <- file.path(annotation_base_dir, "Final UMAP Annotation", "micheli_harmony_umap.rds")
dir.create(dirname(final_umap_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(micheli, final_umap_rds)
cat(green$bold("\n Saved Harmony+UMAP Seurat object to:\n"))
cat(green(final_umap_rds), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Logging key dataset info (metadata, dimensions, clustering, assay, etc.)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\nLogging dataset information...\n"))

log_file <- file.path(log_dir, paste0("Micheli_log_", format(Sys.Date(), "%Y%m%d"), ".txt"))
log_con  <- file(log_file, open = "wt")

writeLines("========= Micheli Dataset Diagnostics =========\n", con = log_con)

# Basic dimensions
writeLines(paste0(" Cells: ", ncol(micheli)), con = log_con)
writeLines(paste0(" Genes: ", nrow(micheli)), con = log_con)

# Assay info
writeLines(paste0("\n Default assay: ", DefaultAssay(micheli)), con = log_con)
writeLines(paste0("Available assays: ", paste(Assays(micheli), collapse = ", ")), con = log_con)

# Reductions present
writeLines(paste0("\n Reductions: ", paste(Reductions(micheli), collapse = ", ")), con = log_con)

# Metadata columns
writeLines("\n Metadata columns:\n", con = log_con)
writeLines(paste(colnames(micheli@meta.data), collapse = "\n"), con = log_con)

# Cluster info
if ("seurat_clusters" %in% colnames(micheli@meta.data)) {
  cluster_counts <- table(micheli$seurat_clusters)
  writeLines("\n Cluster levels:\n", con = log_con)
  writeLines(paste0(levels(factor(micheli$seurat_clusters)), collapse = ", "), con = log_con)
  writeLines("\n Cells per cluster:\n", con = log_con)
  writeLines(capture.output(cluster_counts), con = log_con)
}

# Annotation info
if ("harmony_factorIDs" %in% colnames(micheli@meta.data)) {
  annotation_counts <- table(micheli$harmony_factorIDs)
  writeLines("\n️ Annotation labels (harmony_factorIDs):\n", con = log_con)
  writeLines(paste0(levels(factor(micheli$harmony_factorIDs)), collapse = ", "), con = log_con)
  writeLines("\n Cells per annotation:\n", con = log_con)
  writeLines(capture.output(annotation_counts), con = log_con)
}

# Timepoint + sample info
if ("injury.days" %in% colnames(micheli@meta.data)) {
  day_counts <- table(micheli$injury.days)
  writeLines("\n Timepoints (injury.days):\n", con = log_con)
  writeLines(paste0(levels(factor(micheli$injury.days)), collapse = ", "), con = log_con)
  writeLines("\n Cells per timepoint:\n", con = log_con)
  writeLines(capture.output(day_counts), con = log_con)
}
if ("sample" %in% colnames(micheli@meta.data)) {
  sample_counts <- table(micheli$sample)
  writeLines("\n Samples:\n", con = log_con)
  writeLines(paste0(levels(factor(micheli$sample)), collapse = ", "), con = log_con)
  writeLines("\n Cells per sample:\n", con = log_con)
  writeLines(capture.output(sample_counts), con = log_con)
}

# Top expressed genes
top_genes <- head(sort(Matrix::rowSums(GetAssayData(micheli, slot = "counts")), decreasing = TRUE), 20)
writeLines("\n Top 20 most expressed genes:\n", con = log_con)
writeLines(capture.output(top_genes), con = log_con)

# Preview UMAP embedding (if available)
if ("umap_fresh" %in% Reductions(micheli)) {
  writeLines("\n UMAP (umap_fresh) coordinate preview:\n", con = log_con)
  writeLines(capture.output(head(Embeddings(micheli, "umap_fresh"))), con = log_con)
}

# R Session info
writeLines("\n R Session Info:\n", con = log_con)
writeLines(capture.output(sessionInfo()), con = log_con)

# Close log file
close(log_con)
cat(green$bold(paste("✅ Log saved to:", log_file, "\n")))

# ─────────────────────────────────────────────────────────────────────────────
#                          UMAP by Clusters Only (Harmony-Corrected)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Generating UMAP by cluster IDs only...\n"))

# Set identity to cluster assignments
Idents(micheli) <- "seurat_clusters"

# Plot UMAP based on Harmony embedding
p_cluster <- DimPlot(
  micheli,
  reduction = "umap_harmony_clean",   
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4,
  pt.size = 0.3,
  repel = TRUE
) +
  ggtitle("Micheli – Harmony-Corrected UMAP by Clusters") +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

# Print and save
print(p_cluster)
ggsave(file.path(umap_dir, "Micheli_UMAP_clusters_only_Harmony.png"), plot = p_cluster, width = 9, height = 7, dpi = 300)

cat(green$bold("Saved Harmony-corrected UMAP by clusters only.\n"))

# ─────────────────────────────────────────────────────────────────────────────
#                   UMAP BY EXISTING ANNOTATION SECTION (Harmony-Corrected)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Generating annotated UMAP with single-column legend...\n"))

# Control order of annotation labels
micheli$harmony_factorIDs <- factor(micheli$harmony_factorIDs,
                                     levels = sort(unique(micheli$harmony_factorIDs)))

# UMAP by Cosgrove annotation
p_annot <- DimPlot(
  micheli,
  reduction = "umap_harmony_clean",      
  group.by = "harmony_factorIDs",
  label = TRUE,
  label.size = 3.5,
  pt.size = 0.3,
  repel = TRUE
) +
  ggtitle("Micheli – Harmony-Corrected UMAP by Cosgrove Annotation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))  

# UMAP by sample
p_sample <- DimPlot(
  micheli,
  reduction = "umap_harmony_clean",      
  group.by = "sample",
  pt.size = 0.3
) +
  ggtitle("Micheli – Harmony-Corrected UMAP by Sample") +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

print(p_annot)
print(p_sample)

# Save updated UMAPs
ggsave(file.path(umap_dir, "Micheli_UMAP_annotation_Harmony.png"), plot = p_annot, width = 9, height = 7, dpi = 300)
ggsave(file.path(umap_dir, "Micheli_UMAP_sample_Harmony.png"), plot = p_sample, width = 9, height = 7, dpi = 300)

cat(green$bold("Saved Harmony-corrected UMAPs by annotation and sample.\n"))

# ─────────────────────────────────────────────────────────────────────────────
#               UMAPs split by day (Harmony-Corrected)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Generating UMAPs split by injury.days...\n"))

unique_days <- sort(unique(micheli$injury.days))

for (day in unique_days) {
  cells_day <- WhichCells(micheli, expression = injury.days == day)
  if (length(cells_day) > 0) {
    subset_obj <- subset(micheli, cells = cells_day)
    
    p_day <- DimPlot(
      subset_obj,
      reduction = "umap_harmony_clean",  
      group.by = "harmony_factorIDs",
      label = TRUE,
      label.size = 3.5,
      pt.size = 0.3,
      repel = TRUE
    ) +
      ggtitle(paste0("Micheli – Day ", day, " (Harmony-Corrected UMAP)")) +  
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"),
        panel.grid = element_blank()
      ) +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
    
    fname <- file.path(umap_dir, paste0("Micheli_UMAP_day_", day, "_Harmony.png"))  
    ggsave(fname, plot = p_day, width = 9, height = 7, dpi = 300)
  }
}

cat(green$bold("Saved Harmony-corrected UMAPs split by day.\n"))

# ─────────────────────────────────────────────────────────────────────────────
# Combine UMAPs split by injury.days (Harmony-Corrected)
# ─────────────────────────────────────────────────────────────────────────────

cat(blue$bold("Generating combined UMAP grid by injury.days...\n"))

unique_days <- sort(unique(micheli$injury.days))
umap_list <- list()

for (day in unique_days) {
  cells_day <- WhichCells(micheli, expression = injury.days == day)
  if (length(cells_day) > 0) {
    subset_obj <- subset(micheli, cells = cells_day)
    
    p_day <- DimPlot(
      subset_obj,
      reduction = "umap_harmony_clean",  
      group.by = "harmony_factorIDs",
      label = TRUE,
      label.size = 2.8,
      pt.size = 0.3,
      repel = TRUE
    ) +
      ggtitle(paste("Day", day)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank()
      ) +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
    
    umap_list[[as.character(day)]] <- p_day
  }
}

# Combine plots with shared legend
combined_umap <- wrap_plots(umap_list, ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# Save the combined UMAP grid
ggsave(
  filename = file.path(umap_dir, "Micheli_UMAP_combined_by_day_Harmony.png"),  
  plot = combined_umap,
  width = 12,
  height = 10,
  dpi = 300
)

cat(green$bold("Saved combined Harmony-corrected UMAP grid by day (single legend).\n"))

# ─────────────────────────────────────────────────────────────────────────────
# Save cluster-to-annotation mapping (original-based)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Saving Cosgrove-based cluster-to-annotation table...\n"))

# Make sure the identities are set
Idents(micheli) <- "seurat_clusters"

# Extract mapping: take the *most common* annotation for each cluster
cluster_annotation_map <- micheli@meta.data %>%
  group_by(seurat_clusters) %>%
  count(harmony_factorIDs) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  select(seurat_clusters, harmony_factorIDs) %>%
  arrange(as.numeric(as.character(seurat_clusters)))  

# Save to CSV
write.csv(cluster_annotation_map,
          file.path(prop_dir, "Micheli_cluster_to_annotation_Cosgrove.csv"),
          row.names = FALSE)

cat(green$bold("Saved cluster-to-annotation table (Cosgrove-based).\n"))

# ─────────────────────────────────────────────────────────────────────────────
#                         PROPORTIONS – Original  ANNOTATION
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Calculating overall cell type proportions (Cosgrove annotation)...\n"))

df <- micheli@meta.data

# Overall proportions
overall_props <- df %>%
  count(harmony_factorIDs) %>%
  mutate(Proportion = n / sum(n))

# Save both .csv and .rds
write.csv(overall_props, file.path(prop_dir, "Micheli_overall_proportions_Cosgrove.csv"), row.names = FALSE)
saveRDS(overall_props, file.path(prop_dir, "Micheli_overall_proportions_Cosgrove_raw.rds"))

# Set global factor levels for consistency
all_levels <- levels(factor(overall_props$harmony_factorIDs))
overall_props$harmony_factorIDs <- factor(overall_props$harmony_factorIDs, levels = all_levels)

# fixed color map for harmony_factorIDs (consistent across all plots)
color_map <- setNames(scico(length(all_levels), palette = "batlow"), all_levels)

# Barplot: Overall
p_overall <- ggplot(overall_props, aes(x = harmony_factorIDs, y = Proportion, fill = harmony_factorIDs)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_map) +
  theme_minimal() +
  labs(
    title = "Micheli – Overall Cell Type Proportions (Cosgrove Annotation)",
    x = "Cell Type", y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(barplot_dir, "Micheli_barplot_overall_Cosgrove.png"), plot = p_overall, width = 10, height = 6, dpi = 300)
cat(green$bold("Saved overall barplot (Cosgrove annotation).\n"))

# ─────────────────────────────────────────────────────────────────────────────
# PER SAMPLE barplots 
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Generating per-sample barplots (Cosgrove annotation)...\n"))

sample_props <- df %>%
  count(sample, harmony_factorIDs) %>%
  group_by(sample) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

write.csv(sample_props, file.path(prop_dir, "Micheli_per_sample_proportions_Cosgrove.csv"), row.names = FALSE)

samples <- unique(sample_props$sample)

for (s in samples) {
  sample_df <- filter(sample_props, sample == s)
  sample_df$harmony_factorIDs <- factor(sample_df$harmony_factorIDs, levels = all_levels)
  
  p <- ggplot(sample_df, aes(x = harmony_factorIDs, y = Proportion, fill = harmony_factorIDs)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_map) +
    theme_minimal() +
    labs(
      title = paste("Micheli –", s, "(Cosgrove Annotation)"),
      x = "Cell Type", y = "Proportion"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(barplot_dir, paste0("Micheli_barplot_sample_", s, "_Cosgrove.png")), plot = p, width = 10, height = 6, dpi = 300)
}
cat(green$bold("Saved barplots for each sample (Cosgrove annotation).\n"))

# ─────────────────────────────────────────────────────────────────────────────
# PER DAY barplots 
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("Generating per-day barplots (Cosgrove annotation)...\n"))

day_props <- df %>%
  count(injury.days, harmony_factorIDs) %>%
  group_by(injury.days) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

write.csv(day_props, file.path(prop_dir, "Micheli_per_day_proportions_Cosgrove.csv"), row.names = FALSE)

# Combined facet barplot (days)
p_combined <- ggplot(day_props, aes(x = harmony_factorIDs, y = Proportion, fill = harmony_factorIDs)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_map) +
  facet_wrap(~ injury.days, scales = "free_x") +
  theme_bw(base_size = 9) +
  labs(title = "Micheli – Cell Type Proportions by Day (Cosgrove)", x = "Cell Type", y = "Proportion") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(barplot_dir, "Micheli_barplot_by_day_combined_Cosgrove.png"), plot = p_combined, width = 12, height = 8, dpi = 300)

# Individual plots by day
for (d in unique(day_props$injury.days)) {
  day_df <- filter(day_props, injury.days == d)
  day_df$harmony_factorIDs <- factor(day_df$harmony_factorIDs, levels = all_levels)
  
  p <- ggplot(day_df, aes(x = harmony_factorIDs, y = Proportion, fill = harmony_factorIDs)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_map) +
    theme_minimal() +
    labs(title = paste("Micheli – Day", d, "(Cosgrove Annotation)"), x = "Cell Type", y = "Proportion") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(barplot_dir, paste0("Micheli_barplot_day_", d, "_Cosgrove.png")), plot = p, width = 10, height = 6, dpi = 300)
}

cat(green$bold("Saved barplots for each day (Cosgrove annotation).\n"))


######################################################### SETTING UP THE MANIUAL ANNOTATION PIPELINE #########################################################

# ─────────────────────────────────────────────────────────────────────────────
# PLACEHOLDER — DEFINE Renad_markers BEFORE RUNNING
# ─────────────────────────────────────────────────────────────────────────────

Renad_markers <- list(
  "Endothelial" = c("PECAM1", "CD31", "CDH5", "VEGFR2", "KDR"),
  "Smooth Muscle" = c("ACTA2", "MYL9", "MYH11"),
  "Pericytes" = c("MCAM", "PDGFRB", "NG2", "CSPG4"),
  "Tenocytes" = c("TNMD", "SCX", "COL1A1"),
  "Mature" = c("MYH1", "MYH2", "MYH4", "ACTA1"),
  "FAPs (Adipogenic)" = c("PPARG", "CD36"),
  "FAPs (Stem)" = c("PDGFRA", "SCA1"),
  "FAPs (Pro-remodeling)" = c("TGFB", "MMP14"),
  "FAPs (Cxc14+)" = c("CXCL14"),
  "Neutrophils" = c("S100A8", "S100A9", "CXCR2", "CCL3"),
  "M0M Monocytes" = c("CTSA", "CTSL", "CTSB", "CCR2", "CD14", "ITGAM", "LYZ2", "S100A4", "CCL2"),
  "Ly6C+ Monocytes" = c("CD177", "LY6C1", "LY6C2C", "IRF5"), 
  "Classical M1 Macrophages" = c("CLL2", "CCL6", "CCL9", "TNF", "IL10", "NOS2"),
  "M1a Macrophages" = c("SPP1", "IL17RA", "CXCL3", "CXCL16"),
  "M1b Macrophages" = c("CCR2", "CXCL16", "C1QC", "IL6RA", "CSF1R", "AIF1"),
  "CD14+ Macrophages" = c("CD14", "SLFN4", "CXCL2", "FOSL2", "ATF4", "MDM2"),
  "CX3CR1+ Macrophages" = c("CX3CR1", "LY86", "LPXN", "DHRS3", "DNMT3A"),
  "M2 Macrophages" = c("ARG1", "CD163", "MRC1", "CCL8", "FOLR2", "APOE", "TREM2"),
  "T cells" = c("CD3E", "CD4", "CCL5", "CD8A"),
  "NK cells" = c("NKG7", "KLRA7", "NCAM1", "CD56"),
  "B cells" = c("CD19", "MS4A1", "CD20"),
  "Dendritic cells" = c("CD1C", "HLA-DR", "HLADR", "CD11C", "ITGAX"),
  "Dendritic_own1" = c( "ITGAX", "H2-AB1", "H2-AA", "FLT3", "ZBTB46", "CD74", "CD209A"), 
  "Dendritic_own2" = c( "XCR1", "CLEC9A", "CD8A", "BATF3", "IRF8", "SIRPA", "IRF4", "ITGAM"), 
  "Dendritic_own3" = c( "SIGLECH", "BST2", "IRF7", "TCF4", "LY6D", "SPIB" ), 
  "QSCs" = c("CHODL", "HEYL", "PAX7", "SDC3"),
  "Self-Renewing QSCs" = c("NOTCH2", "NOTCH3", "DPT", "CALCR", "COL15A1", "SPRY1", "NOTCH2/3"),
  "ASCs (early activated)" = c("MYOD1", "SDC4", "TCF7", "CARM1", "MAPK", "FOS", "MEST"),
  "IMB" = c("CTSB", "IFIT1", "IFIT3", "ISG15", "CXCL"),
  "Proliferating ASCs" = c("KI67", "RAC1", "TOP2A", "EZH2", "BIRC5", "CDK1"),
  "Myocytes & Differentiating" = c("ACTA1", "MYOG", "TTN", "MYH", "MYL4")
)

Renad_markers_titlecase <- lapply(Renad_markers, stringr::str_to_title)

# ─────────────────────────────────────────────────────────────────────────────
#                         Marker Gene Presence Check
# ─────────────────────────────────────────────────────────────────────────────

cat(blue$bold("\nRunning marker gene presence check...\n"))

check_marker_coverage <- function(marker_list, seurat_obj, dataset_name) {
  cat(blue$bold(paste0("\n===== Checking ", dataset_name, " =====\n")))
  
  results <- data.frame(Cell_Type = character(), Total_Markers = numeric(),
                        Found = numeric(), Missing = numeric(), stringsAsFactors = FALSE)
  
  found_marker_list <- list()
  
  for (cell_type in names(marker_list)) {
    genes <- marker_list[[cell_type]]
    found_genes <- genes[genes %in% rownames(seurat_obj)]
    
    results <- rbind(results, data.frame(
      Cell_Type = cell_type,
      Total_Markers = length(genes),
      Found = length(found_genes),
      Missing = length(genes) - length(found_genes)
    ))
    
    found_marker_list[[cell_type]] <- found_genes
    
    cat(green$bold(sprintf("%-30s", cell_type)),
        sprintf("– Found: %2d / %2d", length(found_genes), length(genes)),
        if (length(found_genes) == 0) red(" ❌")
        else if (length(found_genes) < length(genes)) yellow(" ⚠️")
        else green(" ✅"), "\n")
  }
  
  return(list(summary_table = results, found_markers = found_marker_list))
}

# ── Run check and save output ───────────────────────────────────────────────

marker_output     <- check_marker_coverage(Renad_markers_titlecase, micheli, "Micheli – RNA assay")
marker_summary    <- marker_output$summary_table
found_marker_list <- marker_output$found_markers

# ── Save to log file ────────────────────────────────────────────────────────

marker_log_path <- file.path(log_dir, "marker_presence_summary_RNA.txt")
sink(marker_log_path)
cat("=== Marker Gene Presence Summary (RNA assay) ===\n\n")
cat("Micheli:\n\n")
print(marker_summary)
sink()

cat(magenta$bold("\n Marker presence summary saved to:\n"))
cat(magenta(marker_log_path), "\n")

# ── Updated Plotting Function ───────────────────────────────────────────────

generate_combined_plots <- function(
    seurat_obj,
    marker_list,
    plots_dir,
    target_cell_type = NULL,
    target_marker = NULL
) {
  # Create base output directory
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  cat(blue$bold("\n Base plot folder created at:\n"))
  cat(blue(plots_dir), "\n")
  
  # Determine which cell types to plot
  cell_types_to_plot <- if (!is.null(target_cell_type)) {
    if (!(target_cell_type %in% names(marker_list))) {
      stop(red$bold(sprintf("Cell type '%s' not found in marker list.\n", target_cell_type)))
    }
    target_cell_type
  } else {
    names(marker_list)
  }
  
  for (cell_type in cell_types_to_plot) {
    all_markers <- marker_list[[cell_type]]
    
    markers <- if (!is.null(target_marker)) {
      if (!(target_marker %in% all_markers)) {
        stop(red$bold(sprintf("Marker '%s' not found in '%s' marker list.\n", target_marker, cell_type)))
      }
      target_marker
    } else {
      all_markers
    }
    
    if (length(markers) == 0) {
      cat(red$bold(sprintf("\n  Skipping '%s' — No valid markers found.\n", cell_type)))
      next
    }
    
    # Create output folder
    cell_dir <- file.path(plots_dir, gsub("[ /]", "_", cell_type))
    if (!dir.exists(cell_dir)) dir.create(cell_dir)
    
    if (!is.null(target_marker)) {
      cell_dir <- file.path(cell_dir, target_marker)
      dir.create(cell_dir, showWarnings = FALSE)
    }
    
    cat(yellow$bold(sprintf("\n Generating plots for: %s\n", cell_type)))
    cat(yellow(sprintf("   Markers: %s\n", paste(markers, collapse = ", "))))
    
    DefaultAssay(seurat_obj) <- "RNA"
    
    # ─ Feature Plot ─
    p_feat_list <- lapply(markers, function(gene) {
      tryCatch({
        FeaturePlot(seurat_obj, features = gene, reduction = "umap_harmony_clean", pt.size = 0.3)+
          ggtitle(gene)
      }, error = function(e) {
        ggplot() + ggtitle(paste(gene, "(failed)"))
      })
    })
    p_feat <- wrap_plots(p_feat_list) + plot_annotation(title = paste(cell_type, "– Feature Plots"))
    ggsave(file.path(cell_dir, "FeaturePlot.png"), p_feat, width = 12, height = 6)
    
    # ─ Violin Plot ─
    p_vln <- wrap_plots(lapply(markers, function(gene) {
      tryCatch(
        VlnPlot(seurat_obj, features = gene, pt.size = 0.1, group.by = "seurat_clusters") +
          ggtitle(gene) + NoLegend(),
        error = function(e) ggplot() + ggtitle(paste(gene, "(failed)"))
      )
    })) + plot_annotation(title = paste(cell_type, "– Violin Plots"))
    ggsave(file.path(cell_dir, "ViolinPlot.png"), p_vln, width = 12, height = 6)
    
    # ─ Dot Plot ─
    p_dot <- tryCatch({
      DotPlot(seurat_obj, features = markers, group.by = "seurat_clusters") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste(cell_type, "– Dot Plot"))
    }, error = function(e) {
      ggplot() + ggtitle("DotPlot failed")
    })
    ggsave(file.path(cell_dir, "DotPlot.png"), p_dot, width = 10, height = 6)
    
    cat(green$bold(sprintf("Saved all 3 plots for '%s'%s\n", cell_type,
                           if (!is.null(target_marker)) paste0(" (", target_marker, ")") else "")))
  }
  
  cat(magenta$bold("\n Plot generation complete!\n"))
}

# ── Loop over all marker sets ───────────────────────────────────────────────

# All plots for all cell types
generate_combined_plots(seurat_obj = micheli,marker_list = Renad_markers_titlecase,plots_dir = renad_plot_dir)

# Only for "Dendritic_Own"
generate_combined_plots(seurat_obj = micheli,marker_list = Renad_markers_titlecase,plots_dir = renad_plot_dir,target_cell_type = "Dendritic_own3")

# ─────────────────────────────────────────────────────────────────────────────
#                 Top Marker Detection + Visualization per Cluster
# ─────────────────────────────────────────────────────────────────────────────

cat(blue$bold("\n Identifying top 6 markers per cluster...\n"))

# Ensure RNA is the default assay for marker detection
DefaultAssay(micheli) <- "RNA"

# Find markers (adjust if needed to avoid very large gene sets)
all_markers <- FindAllMarkers(
  object = micheli,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "RNA"
)

# Save full marker table to CSV
write.csv(all_markers, file = file.path(top_marker_dir, "All_Top_Markers_By_Cluster.csv"), row.names = FALSE)
cat(green$bold("Full marker CSV saved to:\n"))
cat(green(file.path(top_marker_dir, "All_Top_Markers_By_Cluster.csv")), "\n")

# Get top 6 markers per cluster
top6_by_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC)

# Convert to marker list format compatible with your plotting function
top_marker_list <- split(top6_by_cluster$gene, top6_by_cluster$cluster)
names(top_marker_list) <- paste0("Cluster_", names(top_marker_list))  

# Save top 6 markers as a second CSV
top6_path <- file.path(top_marker_dir, "Top6_Markers_Per_Cluster.csv")
write.csv(top6_by_cluster, top6_path, row.names = FALSE)
cat(green$bold("Top 6 markers CSV saved to:\n"))
cat(green(top6_path), "\n")

# ─────────────────────────────────────────────────────────────────────────────
#             Generate Feature/Vln/Dot plots for top 6 markers per cluster
# ─────────────────────────────────────────────────────────────────────────────

cat(blue$bold("\n Generating plots for top 6 markers per cluster...\n"))

generate_combined_plots(
  seurat_obj = micheli,
  marker_list = top_marker_list,
  plots_dir = top_marker_dir
)

cat(magenta$bold("\n Top marker visualization complete.\n"))


