# ==================================================
# Code for annotating and displaying the 
#          Myeloid subcluster   
# ==================================================

# ====================================================================== Libraries and Dataset  ======================================================================

suppressPackageStartupMessages({
  library(pryr)
  library(dplyr)
  library(tidyr)
  library(Seurat)
  library(plotly)
  library(readxl)
  library(harmony)
  library(stringr)
  library(ggplot2)
  library(celldex)
  library(SingleR)
  library(viridis)
  library(ggrastr)
  library(openxlsx)
  library(clustree)
  library(pheatmap)
  library(patchwork)
  library(effectsize)
  library(RColorBrewer)
  library(DoubletFinder)
  library(SingleCellExperiment)
})

# ====================================================================== Annotation Setup   ======================================================================

#  ==  Load the Final Harmony-Processed Seurat Object
cat("Loading the Harmony-processed Seurat object...\n")
final_harmony_file <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony/Final_Seurat_Harmony_Clustered.rds"
seurat_obj <- readRDS(final_harmony_file)
cat("Successfully loaded the Seurat object.\n")

# ====================================================================== Overall Seurat Object Annotation   ======================================================================

#  == manual annotations for the clusters directly with updated labels
manual_annotations <- c(
  "Myeloid Cells",               # Cluster 0
  "FAPs (Cxc14+)",               # Cluster 1
  "Myeloid Cells",               # Cluster 2
  "B Cells",                     # Cluster 3
  "QSCs",                        # Cluster 4
  "FAPs (Cxc14+)",               # Cluster 5
  "T + NK Cells",                # Cluster 6
  "Endothelial Cells",           # Cluster 7
  "Myeloid Cells",               # Cluster 8
  "Prlf.ASCs",                   # Cluster 9
  "Neutrophils",                 # Cluster 10
  "Myeloid Cells",               # Cluster 11
  "FAPs (Stem)",                 # Cluster 12
  "Keratinocytes",               # Cluster 13
  "Prlf.ASCs",                   # Cluster 14
  "Myocytes",                    # Cluster 15
  "Tenocytes",                   # Cluster 16 
  "ASCs",                        # Cluster 17
  "Osteoclasts",                 # Cluster 18 
  "Prlf.ASCs",                   # Cluster 19
  "Myeloid Cells",               # Cluster 20
  "Myeloid Cells",               # Cluster 21
  "Smooth Muscle",               # Cluster 22
  "FAPs (Adipogenic)",           # Cluster 23
  "QSCs",                        # Cluster 24
  "B Cells",                     # Cluster 25
  "FAPs (Cxc14+)"                # Cluster 26
)

#  == Convert `seurat_clusters` to numeric and add manual annotations
seurat_obj$cell_annotation <- factor(
  manual_annotations[as.numeric(as.character(seurat_obj$seurat_clusters)) + 1], 
  levels = unique(manual_annotations)
)

#  == Generate the UMAP plot with the new annotations
umap_plot <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_annotation",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("UMAP with Manual Cell Annotations") +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#  == Define UMAP directory
umap_output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Myeloid Typing/Main UMAP"

#  == Ensure the directory exists
if (!dir.exists(umap_output_dir)) {
  dir.create(umap_output_dir, recursive = TRUE)
}

#  == Save the UMAP 
umap_file <- file.path(umap_output_dir, "Manual_Annotation_UMAP.png")
ggsave(filename = umap_file, plot = umap_plot, width = 12, height = 10)

print(umap_plot)

cat("UMAP with manual annotations saved to:", umap_file, "\n")

# ====================================================================== Subclustering Myeloid Subc cluster  ======================================================================

#  == Subset the Data for Myeloid Cells
cat("Subsetting data for Myeloid Cells...\n")
myeloid_cells <- subset(seurat_obj, idents = c(0, 2, 8, 11, 20, 21))

#  == Re-Normalize and Scale the Data
myeloid_cells <- NormalizeData(myeloid_cells)
myeloid_cells <- FindVariableFeatures(myeloid_cells)
myeloid_cells <- ScaleData(myeloid_cells)

#  == Perform PCA
myeloid_cells <- RunPCA(myeloid_cells, npcs = 20)

#  == Run Harmony for Batch Correction
myeloid_cells <- RunHarmony(
  object = myeloid_cells,
  group.by.vars = c("strain", "time_point", "replicate")
)

#  == Perform Clustering and UMAP at Resolution 0.3
resolution <- 0.3
output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Myeloid Typing/Subcluster UMAP_0.3"

#  == Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#  == Find Neighbors and Clusters
myeloid_cells <- FindNeighbors(myeloid_cells, reduction = "harmony", dims = 1:20)
myeloid_cells <- FindClusters(myeloid_cells, resolution = resolution)
myeloid_cells <- RunUMAP(myeloid_cells, reduction = "harmony", dims = 1:20)

#  == Save UMAP 
cat("Saving UMAP plot to the output directory...\n")
umap_plot <- DimPlot(myeloid_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("Myeloid Cells Subclusters UMAP")

png_filename <- file.path(output_dir, "myeloid_cells_Subclusters_UMAP.png")
ggsave(filename = png_filename, plot = umap_plot, width = 8, height = 6, dpi = 300)
cat("UMAP plot saved successfully at:", png_filename, "\n")

# ====================================================================== Define Myeloid Marker Sets with Normalization   ======================================================================

#  ==  myeloid marker genes (From Renad's list)
myeloid_markers <- list(
  "M1_Macrophages" = c("CCL6", "CCL9", "TNF", "NOS2", "SPP1", "IL10", "CXCL3", "CXCL16", "CCR2", "SLFN4", "CXCL2", "FOSL2", "ATF4", "MDM2", "CX3CR1", "LY86", "LPXN", "DHRS3", "DNMT3A"),
  "M2_Macrophages" = c("ARG1", "CD163", "MRC1", "CCL8", "FOLR2", "APOE", "TREM2"),
  "Monocytes" = c("CTSA", "CTSL", "CTSB", "CCR2", "ITGAM", "LYZ2", "S100A4","CD177", "LY6C1", "LY6C2C", "IRF5"),
  "IMB_Cells" = c("CTSB", "IFIT1", "ISG15", "CXCL"),
  "Neutrophils" = c("S100A8", "S100A9", "CXCR2", "CCL3"),
  "T_Cells" = c("CD3E", "CD4", "CCL5", "CD8A"),
  "Dendritic_Cells" = c("Xcr1", "Batf3", "Clec9a", "Irf8", "Tlr3", "Cd24a", "Sirpa", "Irf4", "Siglech", "Tcf4", "Ly6d", "Tlr7", "Zbtb46", "Itgax")
)

#  == Title Case
title_case_markers <- lapply(myeloid_markers, stringr::str_to_title)

#  == Define directory for saving The OG Renad Markers plots 
renad_markers_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Renad Markers"

#  == Ensure the directory exists
if (!dir.exists(renad_markers_dir)) {
  dir.create(renad_markers_dir, recursive = TRUE)
}

#   ==  Function to generate the annotation verification plots 
generate_combined_plots <- function(seurat_object, cell_type, group.by = "seurat_clusters", reduction = "umap", output_dir) {
  #  == Check if cell type exists in the marker list
  if (!cell_type %in% names(Renad_markers)) {
    stop("Cell type not found in marker list ;(")
  }
  
  #  == Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #  == Get the markers for the selected cell type
  markers <- Renad_markers[[cell_type]]
  
  #  == Generate violin plot
  vln_plot <- VlnPlot(
    object = seurat_object,
    features = markers,
    group.by = group.by,
    pt.size = 0.5,
    ncol = 3
  ) +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
  
  vln_plot <- vln_plot +
    plot_annotation(title = paste(cell_type, "Marker Expression"))  # Add global title using patchwork
  
  #  == Generate feature plot
  feature_plot <- FeaturePlot(
    object = seurat_object,
    features = markers,
    reduction = reduction,
    cols = c("lightgrey", "blue")
  ) +
    theme(
      text = element_text(size = 10)
    )
  
  feature_plot <- feature_plot +
    plot_annotation(title = paste(cell_type, "Marker Expression"))  # Add global title using patchwork
  
  #  == Save plots to specified directory
  violin_file <- file.path(output_dir, paste0(cell_type, "_ViolinPlot.png"))
  feature_file <- file.path(output_dir, paste0(cell_type, "_FeaturePlot.png"))
  
  ggsave(filename = violin_file, plot = vln_plot, width = 12, height = 10)
  ggsave(filename = feature_file, plot = feature_plot, width = 12, height = 10)
  
  #  ==  immediate visualization
  return(list(violin_plot = vln_plot, feature_plot = feature_plot))
}

# ====================================================================== Examine Multiple Cell Types (OG and Combined) ======================================================================

#  == automate generating plots for all cell types in the marker list
for (cell_type in names(myeloid_markers)) {
  print(paste("Generating plots for:", cell_type))
  
  #  == Generate and save plots
  tryCatch({
    plots <- generate_combined_plots(
      seurat_object = myeloid_cells,
      cell_type = cell_type,
      group.by = "seurat_clusters",
      reduction = "umap",
      renad_markers_dir 
    )
    
    # Display plots interactively
    print(plots$violin_plot)
    print(plots$feature_plot)
  }, error = function(e) {
    message(paste("Error generating plots for cell type:", cell_type, "-", e$message))
  })
}

# ====================================================================== Compute Enrichment Scores ======================================================================

#  == Filter markers to include only genes in the dataset
filtered_markers <- lapply(title_case_markers, function(markers) {
  markers[markers %in% rownames(myeloid_cells)]
})

#  == Compute enrichment scores with a reduced control gene adjustment
for (cell_type in names(filtered_markers)) {
  valid_markers <- filtered_markers[[cell_type]]
  
  if (length(valid_markers) > 0) {
    myeloid_cells <- AddModuleScore(
      object = myeloid_cells,
      features = list(valid_markers),
      name = cell_type,
      ctrl = 2 * length(valid_markers)  # Reduced control gene adjustment
    )
  } else {
    cat(paste("Skipping", cell_type, "- No valid markers found.\n"))
  }
}

# ====================================================================== Z-score Normalization and Compute Average Enrichment Scores  ======================================================================

normalize_zscore <- function(scores) {
  return(scale(scores))  
}

for (cell_type in names(filtered_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(myeloid_cells@meta.data)) {
    myeloid_cells@meta.data[[score_col]] <- normalize_zscore(myeloid_cells@meta.data[[score_col]])
  }
}

cat("z-score normalization applied.\n")

#  == Compute average enrichment scores after normalization
enrichment_scores_z <- data.frame(
  Cell_Barcode = rownames(myeloid_cells@meta.data),
  Cluster = myeloid_cells$seurat_clusters
)

for (cell_type in names(filtered_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(myeloid_cells@meta.data)) {
    enrichment_scores_z[[cell_type]] <- myeloid_cells@meta.data[[score_col]]
  }
}

average_enrichment_scores_z <- enrichment_scores_z %>%
  group_by(Cluster) %>%
  summarise(across(-Cell_Barcode, \(x) mean(x, na.rm = TRUE)))

#  == Save normalized averages
z_output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Myeloid Typing/Z_Score_Analysis"
if (!dir.exists(z_output_dir)) dir.create(z_output_dir, recursive = TRUE)

z_scores_file <- file.path(z_output_dir, "ZScore_Average_Enrichment_Scores_By_Cluster.csv")
write.csv(average_enrichment_scores_z, file = z_scores_file, row.names = FALSE)

# ======================================================================  Min-Max Normalization and Compute Average Enrichment Scores  ======================================================================

normalize_minmax <- function(scores) {
  return((scores - min(scores, na.rm = TRUE)) / (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE)))
}

for (cell_type in names(filtered_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(myeloid_cells@meta.data)) {
    myeloid_cells@meta.data[[score_col]] <- normalize_minmax(myeloid_cells@meta.data[[score_col]])
  }
}

#  ==  average enrichment scores after Min-Max normalization
enrichment_scores_mm <- data.frame(
  Cell_Barcode = rownames(myeloid_cells@meta.data),
  Cluster = myeloid_cells$seurat_clusters
)

for (cell_type in names(filtered_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(myeloid_cells@meta.data)) {
    enrichment_scores_mm[[cell_type]] <- myeloid_cells@meta.data[[score_col]]
  }
}

average_enrichment_scores_mm <- enrichment_scores_mm %>%
  group_by(Cluster) %>%
  summarise(across(-Cell_Barcode, \(x) mean(x, na.rm = TRUE)))

#  == Save Min-Max normalized averages
mm_output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Myeloid Typing/Min_Max_Analysis"
if (!dir.exists(mm_output_dir)) dir.create(mm_output_dir, recursive = TRUE)

mm_scores_file <- file.path(mm_output_dir, "MinMax_Average_Enrichment_Scores_By_Cluster.csv")
write.csv(average_enrichment_scores_mm, file = mm_scores_file, row.names = FALSE)

# ====================================================================== Categorize Enrichment Scores into High, Medium, Low Groups ======================================================================

cat("Categorizing enrichment scores into quartiles...\n")

categorize_quartiles <- function(scores) {
  quantiles <- quantile(scores, probs = c(0.25, 0.75), na.rm = TRUE)
  return(ifelse(scores < quantiles[1], "Low", ifelse(scores > quantiles[2], "High", "Medium")))
}

for (cell_type in names(filtered_markers)) {
  score_col <- paste0(cell_type, "1")
  if (score_col %in% colnames(myeloid_cells@meta.data)) {
    myeloid_cells@meta.data[[paste0(cell_type, "_Category")]] <- categorize_quartiles(myeloid_cells@meta.data[[score_col]])
  }
}

#  == Save categorized enrichment scores 
category_output_file <- file.path(output_dir, "Categorized_Enrichment_Scores.csv")
categorized_scores <- myeloid_cells@meta.data[, grepl("_Category$", colnames(myeloid_cells@meta.data))]
categorized_scores <- cbind(Cell_Barcode = rownames(myeloid_cells@meta.data), categorized_scores)
write.csv(categorized_scores, file = category_output_file, row.names = FALSE)

# ====================================================================== Visualize Categorized Enrichment Scores ======================================================================

#  == Convert categorized scores into long format 
categorized_long <- myeloid_cells@meta.data %>%
  select(matches("_Category$"), seurat_clusters) %>%
  pivot_longer(-seurat_clusters, names_to = "Cell_Type", values_to = "Category")

#  == Clean cell type names
categorized_long$Cell_Type <- gsub("_Category", "", categorized_long$Cell_Type)

#  == Generate stacked bar plot
ggplot(categorized_long, aes(x = as.factor(seurat_clusters), fill = Category)) +
  geom_bar(position = "fill") +  # "fill" makes it a percentage stacked bar plot
  facet_wrap(~ Cell_Type, scales = "free_x") +
  labs(title = "Proportion of Enrichment Categories Per Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Enrichment Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

category_summary <- categorized_long %>%
  group_by(seurat_clusters, Cell_Type, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)

#  == Save the summarized table
summary_file <- file.path(output_dir, "Cluster_Wise_Category_Distribution.csv")
write.csv(category_summary, file = summary_file, row.names = FALSE)

# ====================================================================== Subcluster Annotation After inspection  ======================================================================

# Assign manual annotations to clusters
manual_annotations <- c(
  "Monocytes",        # Cluster 0
  "M2 Macrophages",    # Cluster 1
  "M1 Macrophages",    # Cluster 2
  "Monocytes",         # Cluster 3
  "IMB Cells",         # Cluster 4
  "M1 Macrophages",    # Cluster 5
  "Neutrophils",       # Cluster 6
  "M2 Macrophages",    # Cluster 7
  "M2 Macrophages",    # Cluster 8
  "M2 Macrophages",    # Cluster 9
  "M2 Macrophages",    # Cluster 10
  "Dendritic Cells",   # Cluster 11
  "T Cells",           # Cluster 12
  "Dendritic Cells",   # Cluster 13
  "Osteoclasts",       # Cluster 14
  "M2 Macrophages"     # Cluster 15
)

#  == Convert `seurat_clusters` to numeric and assign manual annotations
myeloid_cells$final_annotation <- factor(
  manual_annotations[as.numeric(as.character(myeloid_cells$seurat_clusters)) + 1], 
  levels = unique(manual_annotations)
)

#  == Generate the UMAP plot with final cell annotations
umap_plot <- DimPlot(
  object = myeloid_cells,
  reduction = "umap",
  group.by = "final_annotation",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("Final Annotation of Myeloid Cell Subclusters") +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )


#  == Re-define output directory
output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Myeloid Typing/Final Annotation"

#  == Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#  == Define the output path for the UMAP
umap_file <- file.path(output_dir, "Myeloid_Cells_Final_Annotation_UMAP.png")

#  == Save the UMAP plot
ggsave(filename = umap_file, plot = umap_plot, width = 12, height = 10, dpi = 300)

#  == Save the updated Seurat object with final annotations
final_seurat_file <- file.path(output_dir, "Myeloid_Cells_Final_Annotated.rds")
saveRDS(myeloid_cells, file = final_seurat_file)
cat("Final annotated Seurat object saved at:", final_seurat_file, "\n")




