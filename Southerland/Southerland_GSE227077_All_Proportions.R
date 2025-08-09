# ==================================================
# Code shows the Southerland annotation process  
# And computes the proportions that will be used 
# eventually for the comparison with the model 
# ==================================================

# ====================================================================== Load any potential Libraries ======================================================================  

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

# ====================================================================== Setup the directories ====================================================================== 

#  == Define base directory
base_dir <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing"

#  == Existing directories
seurat_dir <- file.path(base_dir, "Seurat_Object")
overall_plots_dir <- file.path(base_dir, "Annotation", "Overall_Annotation", "Plots")
myeloid_plots_dir <- file.path(base_dir, "Annotation", "Myeloid_Annotation", "Plots")

#  == Proportion directories
prop_main_data_dir  <- file.path(base_dir, "Proportions", "Main Annotation", "Data")
prop_main_plots_dir <- file.path(base_dir, "Proportions", "Main Annotation", "Plots")
prop_myel_data_dir  <- file.path(base_dir, "Proportions", "Myeloid Annotation", "Data")
prop_myel_plots_dir <- file.path(base_dir, "Proportions", "Myeloid Annotation", "Plots")

#  == Create all required directories
dir.create(seurat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(overall_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(myeloid_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_main_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_main_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_myel_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_myel_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ====================================================================== Load and Annotate the dataset (Inherited) ====================================================================== 

cat("Loading the Harmony-processed Seurat object...\n")
final_harmony_file <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony/Final_Seurat_Harmony_Clustered.rds"
seurat_obj <- readRDS(final_harmony_file)
cat("Successfully loaded the Seurat object.\n") 

#  == Define manual annotations for the clusters directly
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

#  == Convert seurat_clusters to numeric and add manual annotations
seurat_obj$cell_annotation <- factor(
  manual_annotations[as.numeric(as.character(seurat_obj$seurat_clusters)) + 1], 
  levels = unique(manual_annotations)
)  

#  == Save the updated Seurat object
annotated_rds_path <- file.path(seurat_dir, "Southerland_Overall_Seurat_Annotated.rds")
saveRDS(seurat_obj, annotated_rds_path)
cat("Annotated Seurat object saved to:", annotated_rds_path, "\n")

#  == Generate UMAP
umap_plot <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_annotation",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("UMAP Overall Cell Annotations") +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#  == Save UMAP (PDF and PNG)
umap_pdf_path <- file.path(overall_plots_dir, "UMAP_Manual_Annotation.pdf")
umap_png_path <- file.path(overall_plots_dir, "UMAP_Manual_Annotation.png")

ggsave(filename = umap_pdf_path, plot = umap_plot, width = 8, height = 6)
ggsave(filename = umap_png_path, plot = umap_plot, width = 8, height = 6, dpi = 300)

# ====================================================================== Subcluster and Plot Myeloid (Inherited) ====================================================================== 

#  == Subset the Data for Myeloid Cells
cat("Subsetting data for Myeloid Cells...\n")
myeloid_cells <- subset(seurat_obj, idents = c(0, 2, 8, 11, 20, 21))

#  == Re-Normalize and Scale the Data
myeloid_cells <- NormalizeData(myeloid_cells)
myeloid_cells <- FindVariableFeatures(myeloid_cells)
myeloid_cells <- ScaleData(myeloid_cells)

#  == PCA and Harmony Batch Correction
myeloid_cells <- RunPCA(myeloid_cells, npcs = 20)
myeloid_cells <- RunHarmony(
  object = myeloid_cells,
  group.by.vars = c("strain", "time_point", "replicate")
)

#  == Clustering and UMAP
resolution <- 0.3
myeloid_cells <- FindNeighbors(myeloid_cells, reduction = "harmony", dims = 1:20)
myeloid_cells <- FindClusters(myeloid_cells, resolution = resolution)
myeloid_cells <- RunUMAP(myeloid_cells, reduction = "harmony", dims = 1:20)

#  == UMAP of Clusters
umap_plot <- DimPlot(
  myeloid_cells,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("Myeloid Cells Subclusters UMAP")

#  == Save UMAP (pre-annotation)
umap_path <- file.path(myeloid_plots_dir, "Myeloid_Subclusters_UMAP.png")
ggsave(filename = umap_path, plot = umap_plot, width = 8, height = 6, dpi = 300)
cat("UMAP plot saved successfully at:", umap_path, "\n")

# ====================================================================== Subcluster Annotation Myeloid ====================================================================== 

#  == Assign manual annotations to clusters
manual_annotations <- c(
  "Monocytes",         # Cluster 0
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

#  == Apply Manual Annotations
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

#  == Save final annotation UMAP
final_umap_path <- file.path(myeloid_plots_dir, "Myeloid_Cells_Final_Annotation_UMAP.png")
ggsave(filename = final_umap_path, plot = umap_plot, width = 12, height = 10, dpi = 300)
cat("Final annotation UMAP saved at:", final_umap_path, "\n")

#  == Save final annotated Seurat object
final_rds_path <- file.path(seurat_dir, "Southerland_Myeloid_Seurat_Annotated.rds")
saveRDS(myeloid_cells, file = final_rds_path)
cat("Final annotated Seurat object saved at:", final_rds_path, "\n")

# ====================================================================== Function to plot UMAPs over time ====================================================================== 

generate_timepoint_umaps <- function(seurat_obj, output_dir, label) {
  cat(paste("\n Generating UMAPs across time points for:", label, "...\n"))
  
  #  == Define subdirectory for saving UMAPs
  umap_dir <- file.path(output_dir, paste0("Timepoint_UMAPs_", label))
  if (!dir.exists(umap_dir)) {
    dir.create(umap_dir, recursive = TRUE)
  }
  
  #  == Standardize time_point values for correct ordering
  seurat_obj$time_point <- factor(
    seurat_obj$time_point,
    levels = c("sham", "day1", "day3", "day7"),
    labels = c("Sham", "Day 1", "Day 3", "Day 7"),
    ordered = TRUE
  )
  
  #  == Loop through strains and time points to generate UMAPs
  unique_strains <- unique(seurat_obj$strain)
  unique_timepoints <- levels(seurat_obj$time_point)
  
  for (strain in unique_strains) {
    for (tp in unique_timepoints) {
      #  == Subset data for strain + time point
      subset_obj <- subset(seurat_obj, strain == strain & time_point == tp)
      
      if (ncol(subset_obj) == 0) {
        cat(paste0(" Skipping: No cells for ", strain, " - ", tp, "\n"))
        next
      }
      
      #  == Generate UMAP 
      umap_plot <- DimPlot(
        object = subset_obj, reduction = "umap", group.by = "cell_annotation",
        label = TRUE, label.size = 4, repel = TRUE
      ) +
        ggtitle(paste("UMAP -", strain, "-", tp)) +
        theme_minimal() +
        theme(
          text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10)
        )
      
      #  == file name
      plot_file <- file.path(umap_dir, paste0("UMAP_", strain, "_", tp, "_", label, ".png"))
      
      #  == Save plot
      ggsave(filename = plot_file, plot = umap_plot, width = 10, height = 8, dpi = 300)
      
      cat("Saved UMAP for:", strain, "-", tp, "to", plot_file, "\n")
    }
  }
  
  cat(" All UMAPs saved in:", umap_dir, "\n")
}

#  == Generate for the overall annotations == 
generate_timepoint_umaps(
  seurat_obj = seurat_obj,
  output_dir = overall_plots_dir,
  label = "Overall"
)

#  == Generate for the Myeloid Annotations ==
generate_timepoint_umaps(
  seurat_obj = myeloid_cells,
  output_dir = myeloid_plots_dir,
  label = "Myeloid"
)

# ====================================================================== Compute Cell Type Counts & Proportions (By replicate and Strain) ====================================================================== 

compute_cell_counts_percentages <- function(seurat_obj, output_dir, label) {
  cat(paste("\n Computing cell type counts & percentages for:", label, "...\n"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  metadata_df <- seurat_obj@meta.data
  
  #  == Fix replicate naming issue at the source
  metadata_df$replicate <- sub("REP", "Rep", metadata_df$replicate, ignore.case = TRUE)
  
  #  == Count the number of cells per type
  cell_type_counts <- metadata_df %>%
    group_by(strain, time_point, replicate, cell_annotation) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  #  == Convert to wide format
  cell_type_counts_wide <- cell_type_counts %>%
    pivot_wider(names_from = cell_annotation, values_from = cell_count, values_fill = 0)
  
  #  == Save counts
  counts_file <- file.path(output_dir, paste0("Cell_Type_Counts_", label, ".csv"))
  write.csv(cell_type_counts_wide, file = counts_file, row.names = FALSE)
  cat("Counts saved:", counts_file, "\n")
  
  #  == Compute percentages instead of proportions
  cell_type_percentages <- cell_type_counts %>%
    group_by(strain, time_point, replicate) %>%
    mutate(percentage = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  
  #  == Save percentages
  percentages_file <- file.path(output_dir, paste0("Cell_Type_Percentages_", label, ".csv"))
  write.csv(cell_type_percentages, file = percentages_file, row.names = FALSE)
  cat("Percentages saved:", percentages_file, "\n")

  return(cell_type_percentages)
}

# ====================================================================== Calculate Percentages for both Main and Sub cluster ====================================================================== 

#  == Compute percentages for overall annotations ==
percentages_overall <- compute_cell_counts_percentages(
  seurat_obj = seurat_obj,
  output_dir = prop_main_data_dir,
  label = "Overall"
)

#  == Compute percentages for myeloid subcluster annotations ==
percentages_myeloid <- compute_cell_counts_percentages(
  seurat_obj = myeloid_cells,
  output_dir = prop_myel_data_dir,
  label = "Myeloid"
)

# ====================================================================== Compute percentage for the entire thing ====================================================================== 

compute_collapsed_percentages <- function(seurat_obj, output_dir, label, annotation_col = "cell_annotation") {
  cat(paste("\n Computing collapsed cell type counts & percentages for:", label, "using", annotation_col, "...\n"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  metadata_df <- seurat_obj@meta.data
  
  #  == Standardize time_point format
  metadata_df$time_point <- as.character(metadata_df$time_point)
  metadata_df$time_point <- gsub("day1", "day 1", metadata_df$time_point)
  metadata_df$time_point <- gsub("day3", "day 3", metadata_df$time_point)
  metadata_df$time_point <- gsub("day7", "day 7", metadata_df$time_point)
  
  metadata_df$time_point <- factor(
    metadata_df$time_point,
    levels = c("sham", "day 1", "day 3", "day 7"),
    ordered = TRUE
  )
  
  #  == Count total cells per annotation per time point
  cell_type_counts <- metadata_df %>%
    group_by(time_point, .data[[annotation_col]]) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    mutate(cell_annotation = .data[[annotation_col]])  # Use mutate instead of rename
  
  #  == Convert to wide format
  cell_type_counts_wide <- cell_type_counts %>%
    select(-!!annotation_col) %>%
    pivot_wider(names_from = cell_annotation, values_from = cell_count, values_fill = 0)
  
  #  == Save counts
  counts_file <- file.path(output_dir, paste0("Cell_Type_Counts_Collapsed_", label, ".csv"))
  write.csv(cell_type_counts_wide, file = counts_file, row.names = FALSE)
  cat("✅ Collapsed counts saved:", counts_file, "\n")
  
  #  == Compute percentages
  cell_type_percentages <- cell_type_counts %>%
    group_by(time_point) %>%
    mutate(percentage = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  
  #  == Save percentages
  percentages_file <- file.path(output_dir, paste0("Cell_Type_Percentages_Collapsed_", label, ".csv"))
  write.csv(cell_type_percentages, file = percentages_file, row.names = FALSE)
  cat("✅ Collapsed percentages saved:", percentages_file, "\n")
  
  return(cell_type_percentages)
}

# == Run the function ==

collapsed_overall <- compute_collapsed_percentages(
  seurat_obj = seurat_obj,
  output_dir = prop_main_data_dir,
  label = "Overall"
)

collapsed_myeloid <- compute_collapsed_percentages(
  seurat_obj = myeloid_cells,
  output_dir = prop_myel_data_dir,
  label = "Myeloid",
  annotation_col = "final_annotation"
)

# ====================================================================== Combined plotting function ====================================================================== 

generate_full_combined_barplot <- function(cell_type_percentages, output_dir, label) {
  cat("\n Generating a single full combined stacked bar plot...\n")
  
  #  == Ensure time_point is a character before processing
  cell_type_percentages$time_point <- as.character(cell_type_percentages$time_point)
  
  #  == Standardize time_point format
  cell_type_percentages$time_point <- gsub("day1", "day 1", cell_type_percentages$time_point)
  cell_type_percentages$time_point <- gsub("day3", "day 3", cell_type_percentages$time_point)
  cell_type_percentages$time_point <- gsub("day7", "day 7", cell_type_percentages$time_point)
  
  #  == Convert time_point to ordered factor
  cell_type_percentages$time_point <- factor(
    cell_type_percentages$time_point, 
    levels = c("sham", "day 1", "day 3", "day 7"),
    ordered = TRUE
  )
  
  #  == Create sample_group column for combined x-axis
  cell_type_percentages$sample_group <- paste(
    cell_type_percentages$strain, cell_type_percentages$replicate, cell_type_percentages$time_point,
    sep = "_"
  )

  #  == Define desired sample order
  correct_order <- c(
    "Balbc_Rep1_sham", "Balbc_Rep1_day 1", "Balbc_Rep1_day 3", "Balbc_Rep1_day 7",
    "Balbc_Rep2_sham", "Balbc_Rep2_day 1", "Balbc_Rep2_day 3", "Balbc_Rep2_day 7",
    "Bl6_Rep1_sham", "Bl6_Rep1_day 1", "Bl6_Rep1_day 3", "Bl6_Rep1_day 7",
    "Bl6_Rep2_sham", "Bl6_Rep2_day 1", "Bl6_Rep2_day 3", "Bl6_Rep2_day 7"
  )
  
  #  == Convert sample_group to factor
  cell_type_percentages$sample_group <- factor(
    cell_type_percentages$sample_group,
    levels = correct_order,
    ordered = TRUE
  )

  #  == Safety check
  if (nrow(cell_type_percentages) == 0) {
    stop("\n ERROR: The dataset is empty after ordering! Please check your input data.\n")
  }
  
  #  == Generate the plot
  plot <- ggplot(cell_type_percentages, aes(x = sample_group, y = percentage, fill = cell_annotation)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = paste("Full Cell Type Percentages -", label),
      x = "Sample Group (Strain + Replicate + Time Point)",
      y = "Percentage (%)",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  #  == Save plot
  plot_file <- file.path(output_dir, paste0("Full_Combined_Cell_Type_Barplot_", label, ".png"))
  ggsave(plot_file, plot = plot, width = 14, height = 8)
  cat("Full combined bar plot saved:", plot_file, "\n")
}

# == Run the function ==

# For Overall Data 
generate_full_combined_barplot(
  cell_type_percentages = percentages_overall,
  output_dir = prop_main_plots_dir,
  label = "Overall"
)

# For Myeloid Data 
generate_full_combined_barplot(
  cell_type_percentages = percentages_myeloid,
  output_dir = prop_myel_plots_dir,
  label = "Myeloid"
)

# ====================================================================== Final Plotting function ====================================================================== 

generate_collapsed_barplot <- function(cell_type_percentages, output_dir, label) {
  cat("\n Generating collapsed bar plot (already collapsed input)...\n")
  
  #  == Standardize time_point values
  cell_type_percentages$time_point <- as.character(cell_type_percentages$time_point)
  cell_type_percentages$time_point <- gsub("day1", "day 1", cell_type_percentages$time_point)
  cell_type_percentages$time_point <- gsub("day3", "day 3", cell_type_percentages$time_point)
  cell_type_percentages$time_point <- gsub("day7", "day 7", cell_type_percentages$time_point)
  
  #  == Convert to ordered factor
  cell_type_percentages$time_point <- factor(
    cell_type_percentages$time_point,
    levels = c("sham", "day 1", "day 3", "day 7"),
    ordered = TRUE
  )
  
  #  == Generate the plot directly from collapsed percentages
  plot <- ggplot(cell_type_percentages, aes(x = time_point, y = percentage, fill = cell_annotation)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = paste("Collapsed Cell Type Percentages -", label),
      x = "Time Point",
      y = "Percentage (%)",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"))
  
  #  == Save the plot
  plot_file <- file.path(output_dir, paste0("Collapsed_Cell_Type_Barplot_", label, ".png"))
  ggsave(plot_file, plot = plot, width = 10, height = 6)
  cat("Collapsed bar plot saved:", plot_file, "\n")
}

# == Call the Function ==

# For Overall Data
generate_collapsed_barplot(
  cell_type_percentages = collapsed_overall,
  output_dir = prop_main_plots_dir,
  label = "Overall"
)

# For Myeloid Data
generate_collapsed_barplot(
  cell_type_percentages = collapsed_myeloid,
  output_dir = prop_myel_plots_dir,
  label = "Myeloid"
)





