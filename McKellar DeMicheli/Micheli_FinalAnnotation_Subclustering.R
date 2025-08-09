# ===========================================================================
# Apply final cluster annotations and calculate proportions + plots
# ===========================================================================

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

# =======================
#      Directories
# =======================
base_dir              <- "D:/Masters/Processing Codes/Micheli_Cosgrove_Annotation"
annotation_base_dir   <- file.path(base_dir, "Annotation")
annot_dir             <- file.path(annotation_base_dir, "Final UMAP Annotation")
umap_dir              <- file.path(annot_dir, "UMAPs")
barplot_dir           <- file.path(annot_dir, "Barplots")
prop_dir              <- file.path(annot_dir, "Proportions")

dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)

# =======================
#     Load Seurat Object
# =======================
cat(blue$bold("\n Loading Harmony+UMAP Seurat object...\n"))
micheli <- readRDS(file.path(annot_dir, "micheli_harmony_umap.rds"))
DefaultAssay(micheli) <- "RNA"

Reductions(micheli)

# =======================
#  Final Cluster Annotations
# =======================
final_annotations <- c(
  "M2 Macrophages",           # 0
  "M1b Macrophages",          # 1
  "FAPs (Pro-remodeling)",    # 2
  "FAPs (Pro-remodeling)",    # 3
  "Endothelial",              # 4
  "Ly6c Monocytes",           # 5
  "M1a Macrophages",          # 6
  NA,                         # 7
  NA,                         # 8
  "Endothelial",              # 9
  "MuSCs",                    # 10
  "Endothelial",              # 11
  "T & NK cells",             # 12
  "MuSCs",                    # 13
  "M2 Macrophages",           # 14
  NA,                         # 15
  "Dendritic",                # 16
  "FAPs (Pro-remodeling)",    # 17
  "FAPs (Stem)",              # 18
  "MuSCs",                    # 19
  "Smooth Muscle & Pericytes",# 20
  "Tenocytes",                # 21
  "Neutrophils",              # 22
  "B Cells",                  # 23
  "Mature Muscle",            # 24
  NA,                         # 25
  NA,                         # 26
  "Neural",                   # 27
  NA,                         # 28
  "Mo/M Monocytes",           # 29
  NA,                         # 30
  NA,                         # 31
  NA,                         # 32
  "Mature Muscle",            # 33
  NA,                         # 34
  "Endothelial",              # 35
  "Dendritic",                # 36
  "M1b Macrophages"           # 37
)

Idents(micheli) <- "seurat_clusters"
micheli$final_annotation <- factor(
  final_annotations[as.numeric(as.character(micheli$seurat_clusters)) + 1],
  levels = unique(na.omit(final_annotations))
)

# ==========================
#  Generate Annotated UMAP
# ==========================
cat(blue$bold("\n️Generating annotated UMAP plot...\n"))
p_annot <- DimPlot(
  micheli,
  reduction = "umap_harmony_clean",
  group.by = "final_annotation",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) + ggtitle("De Micheli UMAP – Final Annotation")

print(p_annot)

ggsave(file.path(umap_dir, "UMAP_final_annotation.png"), p_annot, width = 12, height = 10, dpi = 300)
ggsave(file.path(umap_dir, "UMAP_final_annotation.pdf"), p_annot, width = 12, height = 10)

cat(green$bold("\n✅ Annotated UMAP saved.\n")) 

# =======================
#  Save Annotated Object
# =======================
annotated_path <- file.path(annot_dir, "micheli_annotated_final.rds")
saveRDS(micheli, annotated_path)
cat(green$bold("\n✅ Annotated object saved to:\n"))
cat(green(annotated_path), "\n")

# ============= UMAP COORDINATES =============================

# ===================================================
#  Export UMAP Coordinates with Annotation (Micheli)
# ===================================================
cat(blue$bold("\n Exporting De Micheli UMAP coordinates with cell types...\n"))

# Extract UMAP coordinates
umap_coords <- as.data.frame(Embeddings(micheli[["umap_harmony_clean"]]))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Add cell type annotations
umap_coords$Cell_Type <- micheli$final_annotation

# Add cell barcodes as a proper column
umap_coords$Cell_Barcode <- rownames(umap_coords)

# Remove rows with NA annotations
umap_coords <- umap_coords[!is.na(umap_coords$Cell_Type), ]

# Reorder columns for clean export
umap_coords <- umap_coords[, c("Cell_Barcode", "UMAP_1", "UMAP_2", "Cell_Type")]

# Save to Excel
excel_path <- file.path(annot_dir, "UMAP_Coordinates_DeMicheli.xlsx")
write_xlsx(umap_coords, excel_path)

##############################################################################
# Section: Calculate and Plot Cell Type Proportions
##############################################################################

# =====================================
#     Proportion Calculation Function
# =====================================

calculate_proportions <- function(seurat_obj, group_by = "injury.days") {
  total_cells <- seurat_obj@meta.data %>%
    group_by_at(group_by) %>%
    summarise(total_cells = n(), .groups = "drop")
  
  cell_counts <- seurat_obj@meta.data %>%
    group_by_at(c(group_by, "final_annotation")) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  proportions <- left_join(cell_counts, total_cells, by = group_by) %>%
    mutate(proportion = cell_count / total_cells) %>%
    select_at(c(group_by, "final_annotation", "proportion"))
  
  return(proportions)
}

# ===========================
#     Calculate Proportions
# ==========================

cat(blue$bold("\n Calculating proportions...\n"))
proportions_by_day <- calculate_proportions(micheli, "injury.days")
proportions_by_sample <- calculate_proportions(micheli, "sample")

write.csv(proportions_by_day, file.path(prop_dir, "Proportions_by_injury_days.csv"), row.names = FALSE)
write.csv(proportions_by_sample, file.path(prop_dir, "Proportions_by_sample.csv"), row.names = FALSE)

cat(green$bold("Proportions CSVs saved.\n"))

# =======================
#     Plotting Functions
# =======================

generate_bar_plot <- function(data, group_by, output_file, title) {
  p <- ggplot(data, aes_string(x = "final_annotation", y = "proportion", fill = "final_annotation")) +
    geom_bar(stat = "identity") +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title) +
    theme(legend.position = "none")
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

generate_line_plot <- function(data, output_file, title) {
  p <- ggplot(data, aes(x = `injury.days`, y = proportion, color = final_annotation, group = final_annotation)) +
    geom_line() + geom_point() +
    theme_minimal() +
    ggtitle(title) +
    labs(x = "Injury Day", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

generate_horizontal_bar_plot <- function(data, group_by, output_file, title) {
  p <- ggplot(data, aes(x = proportion, y = reorder(final_annotation, proportion), fill = final_annotation)) +
    geom_bar(stat = "identity") +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 12),
      legend.position = "none"
    ) +
    ggtitle(title) +
    xlab("Proportion") +
    ylab("Cell Type")
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

# ============================
#     Generate and Save Plots
# ============================
cat(blue$bold("\n Generating all proportion plots...\n"))

generate_bar_plot(proportions_by_day, "injury.days", file.path(barplot_dir, "BarPlot_by_injury_days.png"), "Cell Type Proportions by Injury Day")
generate_bar_plot(proportions_by_sample, "sample", file.path(barplot_dir, "BarPlot_by_sample.png"), "Cell Type Proportions by Sample")
generate_line_plot(proportions_by_day, file.path(barplot_dir, "LinePlot_by_injury_days.png"), "Cell Type Proportions Over Time")

for (cell_type in unique(proportions_by_day$final_annotation)) {
  cell_data <- filter(proportions_by_day, final_annotation == cell_type)
  
  # Skip if no data for this cell type
  if (nrow(cell_data) == 0) next
  
  out_file <- file.path(
    barplot_dir,
    paste0("BarPlot_", str_replace_all(cell_type, "[^A-Za-z0-9]", "_"), "_by_injury_days.png")
  )
  
  generate_bar_plot(cell_data, "injury.days", out_file, paste(cell_type, "- Proportion by Injury Day"))
}

generate_horizontal_bar_plot(proportions_by_day, "injury.days", file.path(barplot_dir, "Horizontal_BarPlot_by_injury_days.png"), "Proportions by Injury Day (Horizontal)")
generate_horizontal_bar_plot(proportions_by_sample, "sample", file.path(barplot_dir, "Horizontal_BarPlot_by_sample.png"), "Proportions by Sample (Horizontal)")

cat(magenta$bold("\n All annotation and proportion analysis complete.\n"))

# ====================================
# Check before subclustering 
# ====================================

# Organization of metadata columns 
head(micheli@meta.data)

# What assays are active 
Assays(micheli)
DefaultAssay(micheli)

# How many cells here? 
table(micheli$seurat_clusters %in% c("10", "13", "19"))

##############################################################################
# Subcluster MuSCs (clusters "10", "13", "19") from Micheli dataset
##############################################################################

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

# =======================
#      Setup
# =======================
base_dir            <- "D:/Masters/Processing Codes/Micheli_Cosgrove_Annotation"
annotation_dir      <- file.path(base_dir, "Annotation", "Final UMAP Annotation")
subcluster_dir      <- file.path(annotation_dir, "MuSC_Subclustering")
dir.create(subcluster_dir, recursive = TRUE, showWarnings = FALSE)

cat(blue$bold("\n Loading annotated Seurat object...\n"))
micheli <- readRDS(file.path(annotation_dir, "micheli_annotated_final.rds"))
DefaultAssay(micheli) <- "RNA"

# Log folder
log_dir <- file.path(subcluster_dir, "Log")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# General plots folder under Annotation
plots_dir <- file.path(subcluster_dir, "Annotation", "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Renad marker-specific plots
renad_plot_dir <- file.path(plots_dir, "Renad Markers")
dir.create(renad_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Top marker-specific plots
top_marker_plot_dir <- file.path(plots_dir, "Top Markers")
dir.create(top_marker_plot_dir, recursive = TRUE, showWarnings = FALSE)

# =======================
#   Subset MuSC Clusters
# =======================
cat(blue$bold("\n Subsetting MuSC clusters (10, 13, 19)...\n"))
musc_obj <- subset(micheli, idents = c("10", "13", "19"))

cat(green$bold("Cells in subset: "), ncol(musc_obj), "\n")  # There's  3221 cells 

# =============================================================
#   SKIP normalization, variable features, scaling, Harmony
# =============================================================

cat(blue$bold("\n Finding subclusters directly (no re-Harmony)...\n"))
musc_obj <- FindNeighbors(musc_obj, dims = 1:20)
musc_obj <- FindClusters(musc_obj, resolution = 0.3) 

cat(blue$bold("\n Running UMAP for MuSC subset...\n"))
musc_obj <- RunUMAP(
  musc_obj,
  dims = 1:20,
  reduction = "pca",
  reduction.name = "umap_harmony_musc",
  reduction.key = "UMAPMuSC_"
)

# =======================
#     UMAP Plot
# =======================
cat(blue$bold("\n️ Generating MuSC UMAP...\n"))
p_umap <- DimPlot(
  musc_obj,
  reduction = "umap_harmony_musc",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("MuSC Subclusters – 0.3 Res")

print(p_umap)

ggsave(
  filename = file.path(subcluster_dir, "MuSCs_UMAP.png"),
  plot = p_umap,
  width = 10, height = 8, dpi = 300
)

# ===================================
#   Day-wise UMAPs (Pre-Annotation)
# ===================================

day_umap_dir <- file.path(subcluster_dir, "UMAPs_By_Day")
dir.create(day_umap_dir, recursive = TRUE, showWarnings = FALSE)

unique_days <- sort(unique(musc_obj$injury.days))

cat(blue$bold("\n Generating day-wise UMAPs (pre-annotation)...\n"))
for (day in unique_days) {
  cat(cyan(sprintf("  → Day %s\n", day)))
  
  day_obj <- subset(musc_obj, subset = injury.days == day)
  
  p_day <- DimPlot(
    day_obj,
    reduction = "umap_harmony_musc",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste("MuSC Subclusters – Day", day, "(Resolution 0.3)"))
  
  ggsave(
    filename = file.path(day_umap_dir, paste0("MuSCs_UMAP_Day_", day, "_Res_0.3.png")),
    plot = p_day,
    width = 10, height = 8, dpi = 300
  )
}

cat(green$bold("\n Day-wise pre-annotation UMAPs saved to:\n"))
cat(green(day_umap_dir), "\n")

# =======================
#     Save Object
# =======================
saveRDS(musc_obj, file = file.path(subcluster_dir, "micheli_MuSCs_subclustered.rds"))
cat(green$bold("\n MuSC subclustering complete and saved.\n"))
cat(cyan(" -> UMAP saved to: "), file.path(subcluster_dir, "MuSCs_UMAP.png"), "\n")
cat(cyan(" -> Object saved to: "), file.path(subcluster_dir, "micheli_MuSCs_subclustered.rds"), "\n")

# ===================================================================== Setting things up for annotation of the subcluster =====================================================================

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

marker_output     <- check_marker_coverage(Renad_markers_titlecase, musc_obj, "Micheli MuSCs – RNA assay")
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
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  cat(blue$bold("\n Base plot folder created at:\n"))
  cat(blue(plots_dir), "\n")
  
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
    
    cell_dir <- file.path(plots_dir, gsub("[ /]", "_", cell_type))
    if (!dir.exists(cell_dir)) dir.create(cell_dir)
    
    if (!is.null(target_marker)) {
      cell_dir <- file.path(cell_dir, target_marker)
      dir.create(cell_dir, showWarnings = FALSE)
    }
    
    cat(yellow$bold(sprintf("\n Generating plots for: %s\n", cell_type)))
    cat(yellow(sprintf("   ➤ Markers: %s\n", paste(markers, collapse = ", "))))
    
    DefaultAssay(seurat_obj) <- "RNA"
    
    # ─ Feature Plot ─
    p_feat_list <- lapply(markers, function(gene) {
      tryCatch({
        FeaturePlot(seurat_obj, features = gene, reduction = "umap_harmony_musc", pt.size = 0.3) +
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
generate_combined_plots(
  seurat_obj = musc_obj,
  marker_list = found_marker_list,  # From presence check
  plots_dir = renad_plot_dir
)

# Only for "Dendritic_Own"
generate_combined_plots(
  seurat_obj = musc_obj,
  marker_list = found_marker_list,
  plots_dir = renad_plot_dir,
  target_cell_type = "Dendritic_own3"
)

# ======================= Subcluster Annotation  =======================

cat(blue$bold("\n Assigning final annotations to MuSC subclusters...\n"))

# Define annotations
muscs_final_annotations <- c(
  "ASCs",                     # Cluster 0
  "QSCs",                     # Cluster 1
  "QSCs",                     # Cluster 2
  "Differentiating Myocytes", # Cluster 2
  "Differentiating Myocytes"  # Cluster 2
)

# Apply to object
Idents(musc_obj) <- "seurat_clusters"

# Add annotation column to metadata
musc_obj$muscs_annotation <- factor(
  muscs_final_annotations[as.numeric(as.character(musc_obj$seurat_clusters)) + 1],
  levels = unique(muscs_final_annotations)
)

# =======================
#  Annotated UMAP Plot
# =======================
cat(blue$bold("\n Generating final annotated MuSC UMAP...\n"))

p_musc_annot <- DimPlot(
  musc_obj,
  reduction = "umap_harmony_musc",
  group.by = "muscs_annotation",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) + ggtitle("De Micheli MuSC Subclusters – Final Annotation")

print(p_musc_annot)

ggsave(
  filename = file.path(subcluster_dir, "MuSCs_UMAP_final_annotation.png"),
  plot = p_musc_annot,
  width = 10, height = 8, dpi = 300
)

ggsave(
  filename = file.path(subcluster_dir, "MuSCs_UMAP_final_annotation.pdf"),
  plot = p_musc_annot,
  width = 10, height = 8
)

cat(green$bold("\n Annotated MuSC UMAP saved to:\n"))
cat(green(file.path(subcluster_dir, "MuSCs_UMAP_final_annotation.png")), "\n")

# ====================================
#  Generate day-wise annotated UMAPs
# ====================================

day_umap_dir <- file.path(plots_dir, "Day_UMAPs")
dir.create(day_umap_dir, recursive = TRUE, showWarnings = FALSE)

unique_days <- sort(unique(musc_obj$injury.days))

cat(blue$bold("\n Generating day-wise annotated MuSC UMAPs...\n"))
for (day in unique_days) {
  cat(cyan(sprintf("  → Day %s\n", day)))
  
  day_obj <- subset(musc_obj, subset = injury.days == day)
  
  p_day <- DimPlot(
    day_obj,
    reduction = "umap_harmony_musc",
    group.by = "muscs_annotation",
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste("MuSC Subclusters –", day))
  
  ggsave(
    filename = file.path(day_umap_dir, paste0("MuSCs_UMAP_Day_", day, ".png")),
    plot = p_day,
    width = 10, height = 8, dpi = 300
  )
}

cat(green$bold("\n Day-wise UMAPs saved to:\n"))
cat(green(day_umap_dir), "\n")

# =======================
#  Save Annotated Object
# =======================
annotated_musc_path <- file.path(subcluster_dir, "micheli_MuSCs_annotated_final.rds")
saveRDS(musc_obj, annotated_musc_path)

cat(green$bold("\n✅ Annotated MuSC object saved to:\n"))
cat(green(annotated_musc_path), "\n")

# ======================== UMAP COORDINATES =============================================

# ===============================================
#  Export MuSC UMAP Coordinates with Annotation
# ===============================================
cat(blue$bold("\n Exporting MuSC UMAP coordinates with cell types...\n"))

# Extract UMAP coordinates from subclustered object
umap_musc <- as.data.frame(Embeddings(musc_obj[["umap_harmony_musc"]]))
colnames(umap_musc) <- c("UMAP_1", "UMAP_2")

# Add annotation
umap_musc$Cell_Type <- musc_obj$muscs_annotation

# Add cell barcode as column
umap_musc$Cell_Barcode <- rownames(umap_musc)

# Remove unannotated cells
umap_musc <- umap_musc[!is.na(umap_musc$Cell_Type), ]

# Reorder columns
umap_musc <- umap_musc[, c("Cell_Barcode", "UMAP_1", "UMAP_2", "Cell_Type")]

# Save to Excel
musc_excel_path <- file.path(subcluster_dir, "Micheli_UMAP_Coordinates_MuSCs.xlsx")
write_xlsx(umap_musc, musc_excel_path)

##############################  Subcluster proprotions ##############################

barplot_dir <- file.path(subcluster_dir, "Proportions", "Barplots")
dir.create(barplot_dir, recursive = TRUE, showWarnings = FALSE)

# ==================================
#     MuSC Proportion Calculation
# ==================================
calculate_musc_proportions <- function(seurat_obj, group_by = "injury.days") {
  total_cells <- seurat_obj@meta.data %>%
    group_by_at(group_by) %>%
    summarise(total_cells = n(), .groups = "drop")
  
  cell_counts <- seurat_obj@meta.data %>%
    group_by_at(c(group_by, "muscs_annotation")) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  proportions <- left_join(cell_counts, total_cells, by = group_by) %>%
    mutate(proportion = cell_count / total_cells) %>%
    select_at(c(group_by, "muscs_annotation", "proportion"))
  
  return(proportions)
}

generate_bar_plot <- function(data, group_by, output_file, title, annotation_col = "final_annotation") {
  p <- ggplot(data, aes(x = .data[[annotation_col]], y = proportion, fill = .data[[annotation_col]])) +
    geom_bar(stat = "identity") +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title) +
    theme(legend.position = "none")
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

generate_line_plot <- function(data, output_file, title, annotation_col = "final_annotation") {
  p <- ggplot(data, aes(x = `injury.days`, y = proportion, color = .data[[annotation_col]], group = .data[[annotation_col]])) +
    geom_line() + geom_point() +
    theme_minimal() +
    ggtitle(title) +
    labs(x = "Injury Day", y = "Proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

generate_horizontal_bar_plot <- function(data, group_by, output_file, title, annotation_col = "final_annotation") {
  p <- ggplot(data, aes(x = proportion, y = reorder(.data[[annotation_col]], proportion), fill = .data[[annotation_col]])) +
    geom_bar(stat = "identity") +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 12),
      legend.position = "none"
    ) +
    ggtitle(title) +
    xlab("Proportion") +
    ylab("Cell Type")
  
  ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300)
}

# =================================
#     Calculate MuSC Proportions
# =================================
cat(blue$bold("\n Calculating MuSC subcluster proportions...\n"))
musc_prop_day <- calculate_musc_proportions(musc_obj, "injury.days")
musc_prop_sample <- calculate_musc_proportions(musc_obj, "sample")

write.csv(musc_prop_day, file.path(barplot_dir, "MuSCs_Proportions_by_injury_days.csv"), row.names = FALSE)
write.csv(musc_prop_sample, file.path(barplot_dir, "MuSCs_Proportions_by_sample.csv"), row.names = FALSE)

cat(green$bold("MuSC proportion CSVs saved.\n"))

# =================================
#     Generate and Save MuSC Plots
# ==================================
cat(blue$bold("\n Generating MuSC proportion plots...\n"))

generate_bar_plot(musc_prop_day, "injury.days", file.path(barplot_dir, "MuSCs_BarPlot_by_injury_days.png"), "MuSC Subcluster Proportions by Injury Day", annotation_col = "muscs_annotation")

generate_bar_plot(musc_prop_sample, "sample", file.path(barplot_dir, "MuSCs_BarPlot_by_sample.png"), "MuSC Subcluster Proportions by Sample", annotation_col = "muscs_annotation")

generate_line_plot(musc_prop_day, file.path(barplot_dir, "MuSCs_LinePlot_by_injury_days.png"), "MuSC Subcluster Proportions Over Time", annotation_col = "muscs_annotation")

for (label in unique(musc_prop_day$muscs_annotation)) {
  subset_data <- filter(musc_prop_day, muscs_annotation == label)
  if (nrow(subset_data) == 0) next
  
  out_file <- file.path(
    barplot_dir,
    paste0("MuSCs_BarPlot_", str_replace_all(label, "[^A-Za-z0-9]", "_"), "_by_injury_days.png")
  )
  
  generate_bar_plot(subset_data, "injury.days", out_file, paste(label, "- Proportion by Injury Day"), annotation_col = "muscs_annotation")
}

generate_horizontal_bar_plot(musc_prop_day, "injury.days", file.path(barplot_dir, "MuSCs_Horizontal_BarPlot_by_injury_days.png"), "MuSC Proportions (Horizontal by Injury Day)", annotation_col = "muscs_annotation")

generate_horizontal_bar_plot(musc_prop_sample, "sample", file.path(barplot_dir, "MuSCs_Horizontal_BarPlot_by_sample.png"), "MuSC Proportions (Horizontal by Sample)", annotation_col = "muscs_annotation")

cat(magenta$bold("\n MMuSC annotation and proportion plots complete.\n"))
