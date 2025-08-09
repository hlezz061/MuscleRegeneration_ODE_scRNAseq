# ==================================================
# Code for the entire Southerland processing  
# Along with annotation for the overall UMAP   
# ==================================================

# ====================================================================== Libraries and Dataset  ======================================================================

# == Load required libraries
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

# == Directory for Saving Dataset Summaries
dataset_summary_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Data/Dataset_Summary"
if (!dir.exists(dataset_summary_dir)) {
  dir.create(dataset_summary_dir, recursive = TRUE)
}

# == Initialize metadata log file
metadata_log_file <- file.path(dataset_summary_dir, "Metadata_Log.txt")
writeLines("Metadata Log\n===============\n", con = metadata_log_file)  # Start the log file

# == Load the Data
counts <- tryCatch({
  Read10X(data.dir = "D:/Masters/Datasets/GSE227077/GSE227075_RAW/supplementary_file/mm10_matrix")
}, error = function(e) {
  stop("Error loading data: ", e$message)
})

# == Create a Seurat Object
seurat_obj <- tryCatch({
  CreateSeuratObject(counts = counts, project = "Mouse_Ischemia", min.cells = 3, min.features = 200)
}, error = function(e) {
  stop("Error creating Seurat object: ", e$message)
})

# == Extract strain, condition, replicate, and time point information
cell_metadata <- data.frame(
  strain = sub("_.*", "", colnames(seurat_obj)),                                   # Extract strain
  condition = sub("^[^_]*_", "", sub("_rep.*", "", colnames(seurat_obj))),         # Extract condition
  replicate = sub(".*_(rep\\d+)_.*", "\\1", colnames(seurat_obj)),                 # Extract replicate
  row.names = colnames(seurat_obj)                                                # Ensure row names match cell barcodes
)
cell_metadata$time_point <- ifelse(grepl("day", cell_metadata$condition), cell_metadata$condition, "sham")  # Derive time point

# == Validate metadata parsing
metadata_output <- capture.output(head(cell_metadata))
cat("\nMetadata Parsing Validation:\n", 
    paste(metadata_output, collapse = "\n"), 
    file = metadata_log_file, 
    append = TRUE)

# == Add the metadata back to the Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_metadata)

# == Save Seurat Object and Metadata
saveRDS(seurat_obj, file.path(dataset_summary_dir, "seurat_obj_with_metadata.rds"))
saveRDS(cell_metadata, file.path(dataset_summary_dir, "parsed_metadata.rds"))

# == Inspect Metadata Columns
metadata_inspection <- capture.output({
  print("Unique Strains:")
  print(unique(seurat_obj@meta.data$strain))
  print("Unique Conditions:")
  print(unique(seurat_obj@meta.data$condition))
  print("Unique Replicates:")
  print(unique(seurat_obj@meta.data$replicate))
  print("Unique Time Points:")
  print(unique(seurat_obj@meta.data$time_point))
})

# == Write metadata inspection to the log file
cat("\nMetadata Inspection:\n", 
    paste(metadata_inspection, collapse = "\n"), 
    file = metadata_log_file, 
    append = TRUE)

# == Explore the Raw Count Data
count_summary <- capture.output({
  print("Dimensions of Count Data (Genes x Cells):")
  print(dim(counts))
  print("Subset of Count Data:")
  print(counts[1:10, 1:10])
})

# == Write raw count data inspection to the log file
cat("\nRaw Count Data Inspection:\n", 
    paste(count_summary, collapse = "\n"), 
    file = metadata_log_file, 
    append = TRUE)

# == Summarize the Raw Data
raw_data_summary <- capture.output({
  print("Summary of Counts per Cell:")
  print(summary(colSums(counts)))
  print("Summary of Genes per Cell:")
  print(summary(rowSums(counts)))
})

# == Write raw data summary to the log file
cat("\nRaw Data Summary:\n", 
    paste(raw_data_summary, collapse = "\n"), 
    file = metadata_log_file, 
    append = TRUE)

#  == Save Basic Summary
summary_file <- file.path(dataset_summary_dir, "Seurat_Object_Summary.txt")
writeLines(capture.output(print(seurat_obj)), con = summary_file)

# ====================================================================== Dataset Division  ======================================================================

#  == Complete Division of the Dataset by Strain, Time Point, and Replicate

#  == Divide by Strain
Balbc <- subset(seurat_obj, subset = strain == "Balbc")
Bl6 <- subset(seurat_obj, subset = strain == "Bl6")

#  == Divide by Time Point for Balbc
Balbc_sham <- subset(Balbc, subset = time_point == "sham")
Balbc_day1 <- subset(Balbc, subset = time_point == "day1")
Balbc_day3 <- subset(Balbc, subset = time_point == "day3")
Balbc_day7 <- subset(Balbc, subset = time_point == "day7")

#  == Divide by Time Point for Bl6
Bl6_sham <- subset(Bl6, subset = time_point == "sham")
Bl6_day1 <- subset(Bl6, subset = time_point == "day1")
Bl6_day3 <- subset(Bl6, subset = time_point == "day3")
Bl6_day7 <- subset(Bl6, subset = time_point == "day7")

#  == Divide by Replicate for Balbc
Balbc_sham_rep1 <- subset(Balbc_sham, subset = replicate == "rep1")
Balbc_sham_rep2 <- subset(Balbc_sham, subset = replicate == "rep2")
Balbc_day1_rep1 <- subset(Balbc_day1, subset = replicate == "rep1")
Balbc_day1_rep2 <- subset(Balbc_day1, subset = replicate == "rep2")
Balbc_day3_rep1 <- subset(Balbc_day3, subset = replicate == "rep1")
Balbc_day3_rep2 <- subset(Balbc_day3, subset = replicate == "rep2")
Balbc_day7_rep1 <- subset(Balbc_day7, subset = replicate == "rep1")
Balbc_day7_rep2 <- subset(Balbc_day7, subset = replicate == "rep2")

#  == Divide by Replicate for Bl6
Bl6_sham_rep1 <- subset(Bl6_sham, subset = replicate == "rep1")
Bl6_sham_rep2 <- subset(Bl6_sham, subset = replicate == "rep2")
Bl6_day1_rep1 <- subset(Bl6_day1, subset = replicate == "rep1")
Bl6_day1_rep2 <- subset(Bl6_day1, subset = replicate == "rep2")
Bl6_day3_rep1 <- subset(Bl6_day3, subset = replicate == "rep1")
Bl6_day3_rep2 <- subset(Bl6_day3, subset = replicate == "rep2")
Bl6_day7_rep1 <- subset(Bl6_day7, subset = replicate == "rep1")
Bl6_day7_rep2 <- subset(Bl6_day7, subset = replicate == "rep2")

#  == Save Subsets
dataset_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Data/Subsets"
if (!dir.exists(dataset_dir)) {
  dir.create(dataset_dir, recursive = TRUE)
}

save_subset <- function(obj, filename) {
  if (!is.null(obj) && ncol(obj) > 0) {
    saveRDS(obj, file.path(dataset_dir, filename))
    cat("Saved:", filename, "\n")
  } else {
    cat("Skipped (empty subset):", filename, "\n")
  }
}

save_subset(Balbc_sham_rep1, "Balbc_sham_rep1.rds")
save_subset(Balbc_sham_rep2, "Balbc_sham_rep2.rds")
save_subset(Balbc_day1_rep1, "Balbc_day1_rep1.rds")
save_subset(Balbc_day1_rep2, "Balbc_day1_rep2.rds")
save_subset(Balbc_day3_rep1, "Balbc_day3_rep1.rds")
save_subset(Balbc_day3_rep2, "Balbc_day3_rep2.rds")
save_subset(Balbc_day7_rep1, "Balbc_day7_rep1.rds")
save_subset(Balbc_day7_rep2, "Balbc_day7_rep2.rds")

save_subset(Bl6_sham_rep1, "Bl6_sham_rep1.rds")
save_subset(Bl6_sham_rep2, "Bl6_sham_rep2.rds")
save_subset(Bl6_day1_rep1, "Bl6_day1_rep1.rds")
save_subset(Bl6_day1_rep2, "Bl6_day1_rep2.rds")
save_subset(Bl6_day3_rep1, "Bl6_day3_rep1.rds")
save_subset(Bl6_day3_rep2, "Bl6_day3_rep2.rds")
save_subset(Bl6_day7_rep1, "Bl6_day7_rep1.rds")
save_subset(Bl6_day7_rep2, "Bl6_day7_rep2.rds")

#  == Function to print summary 
print_seurat_summary <- function(seurat_obj, name) {
  cat("\n\n### Summary for:", name, "###\n")
  if (is.null(seurat_obj) || ncol(seurat_obj) == 0) {
    cat("This subset is NULL or empty (no data available).\n")
  } else {
    cat("Dimensions (Genes x Cells):", dim(seurat_obj), "\n")
    cat("Unique Strains:", unique(seurat_obj@meta.data$strain), "\n")
    cat("Unique Conditions:", unique(seurat_obj@meta.data$condition), "\n")
    cat("Unique Replicates:", unique(seurat_obj@meta.data$replicate), "\n")
    cat("Preview of Metadata:\n")
    print(head(seurat_obj@meta.data, 5))  #  the first 5 rows of metadata
  }
}

#  == List all subsets for validation
seurat_subsets <- list(
  Balbc_sham_rep1, Balbc_sham_rep2,
  Balbc_day1_rep1, Balbc_day1_rep2,
  Balbc_day3_rep1, Balbc_day3_rep2,
  Balbc_day7_rep1, Balbc_day7_rep2,
  Bl6_sham_rep1, Bl6_sham_rep2,
  Bl6_day1_rep1, Bl6_day1_rep2,
  Bl6_day3_rep1, Bl6_day3_rep2,
  Bl6_day7_rep1, Bl6_day7_rep2
)

subset_names <- c(
  "Balbc_sham_rep1", "Balbc_sham_rep2",
  "Balbc_day1_rep1", "Balbc_day1_rep2",
  "Balbc_day3_rep1", "Balbc_day3_rep2",
  "Balbc_day7_rep1", "Balbc_day7_rep2",
  "Bl6_sham_rep1", "Bl6_sham_rep2",
  "Bl6_day1_rep1", "Bl6_day1_rep2",
  "Bl6_day3_rep1", "Bl6_day3_rep2",
  "Bl6_day7_rep1", "Bl6_day7_rep2"
)

#  == Iterate through subsets and print them
for (i in seq_along(seurat_subsets)) {
  print_seurat_summary(seurat_subsets[[i]], subset_names[i])
  # Remove subset object directly from memory
  rm(list = subset_names[i])
  gc()  # Trigger garbage collection
}

cat("\nDataset division and saving completed successfully.\n")

# Just a test for later
test_subset <- readRDS(file.path(dataset_dir, "Balbc_sham_rep1.rds"))
print(dim(test_subset))  # Check dimensions
print(head(test_subset@meta.data))  # Check metadata

# ====================================================================== Quality Control  ======================================================================

#  == Directory for saving Quality Control outputs
qc_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control"
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

qc_analysis <- function(seurat_obj, qc_dir, min_features = 1000, max_features = 4000, max_counts = 20000, max_percent_mt = 15) {
  #  == Create output directory if it doesn't exist
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE)
  }
  
  #  == Extract count data properly for Assay5
  counts <- GetAssayData(seurat_obj, layer = "counts")
  
  #  == Gene counts per cell
  gene_counts_per_cell <- Matrix::colSums(counts > 0)
  png(file.path(qc_dir, "Gene_Counts_per_Cell_Distribution.png"), width = 800, height = 600)
  hist(gene_counts_per_cell, breaks = 50, main = "Distribution of Gene Counts per Cell", xlab = "Gene Counts")
  dev.off()
  
  #  == Cell counts per gene
  cell_counts_per_gene <- Matrix::rowSums(counts > 0)
  png(file.path(qc_dir, "Cell_Counts_per_Gene_Distribution.png"), width = 800, height = 600)
  hist(cell_counts_per_gene, breaks = 50, main = "Distribution of Cell Counts per Gene", xlab = "Cell Counts")
  dev.off()
  
  #  == Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  #  == Visualize QC Metrics
  vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(qc_dir, "Violin_Plots_QC_Metrics.png"), plot = vln_plot, width = 10, height = 6)
  
  #  == Scatter plots
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter_plot <- plot1 + plot2
  ggsave(file.path(qc_dir, "Scatter_Plots_QC_Metrics.png"), plot = scatter_plot, width = 12, height = 6)
  
  #  == Filter Cells Based on QC Metrics
  seurat_obj_filtered <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_features & nFeature_RNA < max_features &
      nCount_RNA > 0 & nCount_RNA < max_counts &
      percent.mt < max_percent_mt
  )
  
  #  == Save filtered object summary
  filtered_summary_file <- file.path(qc_dir, "Filtered_Seurat_Object_Summary.txt")
  writeLines(capture.output(print(seurat_obj_filtered)), con = filtered_summary_file)
  
  #  == Save filtered object
  saveRDS(seurat_obj_filtered, file.path(qc_dir, "seurat_obj_filtered.rds"))
  
  return(seurat_obj_filtered)
}

#   ==  Directory to store the RDS files
og_dataset_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Data/Subsets"

#  == List all RDS files in the directory
rds_files <- list.files(og_dataset_dir, pattern = "\\.rds$", full.names = TRUE)

#  == Load each RDS file into the environment and name the object based on the filename
for (file in rds_files) {
  # Extract the base name (without the directory or file extension)
  object_name <- sub("\\.rds$", "", basename(file))
  
  # Load the RDS file and assign it to the environment with the extracted name
  assign(object_name, readRDS(file))
}

#  == Print loaded objects to confirm
loaded_objects <- ls()
print("The following objects were loaded:")
print(loaded_objects)

#  == Balbc strain - Sham
Balbc_sham_rep1_filtered <- qc_analysis(
  seurat_obj = Balbc_sham_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_sham_rep1_QC"
)

Balbc_sham_rep2_filtered <- qc_analysis(
  seurat_obj = Balbc_sham_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_sham_rep2_QC"
)

#  == Balbc strain - Day 1
Balbc_day1_rep1_filtered <- qc_analysis(
  seurat_obj = Balbc_day1_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day1_rep1_QC"
)
Balbc_day1_rep2_filtered <- qc_analysis(
  seurat_obj = Balbc_day1_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day1_rep2_QC"
)

#  == Balbc strain - Day 3
Balbc_day3_rep1_filtered <- qc_analysis(
  seurat_obj = Balbc_day3_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day3_rep1_QC"
)
Balbc_day3_rep2_filtered <- qc_analysis(
  seurat_obj = Balbc_day3_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day3_rep2_QC"
)

#  == Balbc strain - Day 7
Balbc_day7_rep1_filtered <- qc_analysis(
  seurat_obj = Balbc_day7_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day7_rep1_QC"
)
Balbc_day7_rep2_filtered <- qc_analysis(
  seurat_obj = Balbc_day7_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Balbc_day7_rep2_QC"
)

#  == Bl6 strain - Sham
Bl6_sham_rep1_filtered <- qc_analysis(
  seurat_obj = Bl6_sham_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_sham_rep1_QC"
)
Bl6_sham_rep2_filtered <- qc_analysis(
  seurat_obj = Bl6_sham_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_sham_rep2_QC"
)

#  == Bl6 strain - Day 1
Bl6_day1_rep1_filtered <- qc_analysis(
  seurat_obj = Bl6_day1_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day1_rep1_QC"
)
Bl6_day1_rep2_filtered <- qc_analysis(
  seurat_obj = Bl6_day1_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day1_rep2_QC"
)

#  == Bl6 strain - Day 3
Bl6_day3_rep1_filtered <- qc_analysis(
  seurat_obj = Bl6_day3_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day3_rep1_QC"
)
Bl6_day3_rep2_filtered <- qc_analysis(
  seurat_obj = Bl6_day3_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day3_rep2_QC"
)

#  == Bl6 strain - Day 7
Bl6_day7_rep1_filtered <- qc_analysis(
  seurat_obj = Bl6_day7_rep1,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day7_rep1_QC"
)
Bl6_day7_rep2_filtered <- qc_analysis(
  seurat_obj = Bl6_day7_rep2,
  qc_dir = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Quality Control/Bl6_day7_rep2_QC"
)

print(Balbc_sham_rep1_filtered)  # Check dimensions and metadata
summary(Balbc_sham_rep1_filtered@meta.data$percent.mt)  # Ensure no values > 15
summary(Balbc_sham_rep1_filtered@meta.data$nFeature_RNA)  # Ensure values in range 1000–4000
summary(Balbc_sham_rep1_filtered@meta.data$nCount_RNA)  # Ensure values < 20000

# ====================================================================== NSPCA + Doublet Detection  ======================================================================

#  == Define the combined preprocessing function
combined_preprocessing <- function(seurat_obj_filtered, dataset_name, nfeatures = 3000, ndims = 35, PCs = 1:20) {
  #  == base output directory
  base_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/NSPCA_Doub"
  output_dir <- file.path(base_dir, dataset_name)
  
  #  == Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #  == Normalize Data
  seurat_obj_filtered <- NormalizeData(
    seurat_obj_filtered, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
  
  #  == Find Variable Features
  seurat_obj_filtered <- FindVariableFeatures(
    seurat_obj_filtered, 
    selection.method = "vst", 
    nfeatures = nfeatures
  )
  
  #  == Scale Data
  seurat_obj_filtered <- ScaleData(
    seurat_obj_filtered, 
    features = VariableFeatures(seurat_obj_filtered),
    do.scale = TRUE,
    do.center = TRUE
  )
  
  #  == Run PCA
  seurat_obj_filtered <- RunPCA(
    seurat_obj_filtered, 
    features = VariableFeatures(seurat_obj_filtered)
  )
  
  #  == Save intermediate Seurat object post-PCA
  saveRDS(seurat_obj_filtered, file.path(output_dir, "seurat_obj_postPCA.rds"))
  
  #  == Save PCA Results
  pca_dir <- file.path(output_dir, "PCA")
  if (!dir.exists(pca_dir)) {
    dir.create(pca_dir)
  }
  elbow_plot <- ElbowPlot(seurat_obj_filtered, ndims = ndims)
  ggsave(file.path(pca_dir, "Elbow_Plot.png"), plot = elbow_plot, width = 8, height = 6)
  
  #  == Run Doublet Detection
  doublet_dir <- file.path(output_dir, "Doublet")
  if (!dir.exists(doublet_dir)) {
    dir.create(doublet_dir)
  }
  
  #  == Parameter Sweep
  sweep.res <- paramSweep(seurat_obj_filtered, PCs = PCs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  write.csv(sweep.stats, file.path(doublet_dir, "Sweep_Stats.csv"))
  
  #  == Determine Optimal pK
  bcmvn <- find.pK(sweep.stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  optimal_pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
  
  #  == Expected Doublet Rate
  nExp <- round(0.05 * ncol(seurat_obj_filtered))  # Default 5% doublet rate
  
  #  == Add Doublet Classifications
  seurat_obj_filtered <- doubletFinder(
    seurat_obj_filtered, 
    PCs = PCs, 
    pN = 0.25, 
    pK = optimal_pK, 
    nExp = nExp
  )
  
  #  == Save Seurat object post-doublet detection
  saveRDS(seurat_obj_filtered, file.path(doublet_dir, "seurat_obj_postDoubletDetection.rds"))
  
  #  == Save Doublet Classification Results
  DF_col <- grep("DF.classifications_.*", colnames(seurat_obj_filtered@meta.data), value = TRUE)
  seurat_obj_filtered@meta.data$DF.classifications <- seurat_obj_filtered@meta.data[[DF_col[1]]]
  doublet_table <- table(seurat_obj_filtered@meta.data$DF.classifications)
  write.csv(as.data.frame(doublet_table), file.path(doublet_dir, "Doublet_Classification_Results.csv"))
  
  #  == Visualize doublets and singlets on UMAP
  seurat_obj_filtered <- RunUMAP(seurat_obj_filtered, dims = 1:20)
  doublet_umap <- DimPlot(seurat_obj_filtered, group.by = "DF.classifications", label = TRUE, pt.size = 0.5) +
    labs(title = "UMAP Distribution of Doublets and Singlets")
  ggsave(file.path(doublet_dir, "UMAP_Doublets_vs_Singlets.png"), plot = doublet_umap, width = 10, height = 6)
  
  #  == Visualize QC metrics for doublets and singlets
  doublet_vln_plot <- VlnPlot(
    seurat_obj_filtered, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    group.by = "DF.classifications", 
    pt.size = 0.1
  ) +
    labs(title = "QC Metrics for Doublets vs. Singlets")
  ggsave(file.path(doublet_dir, "QC_Metrics_Doublets_vs_Singlets.png"), plot = doublet_vln_plot, width = 10, height = 6)
  
  #  == Filter to retain only singlets
  seurat_obj_singlets <- subset(seurat_obj_filtered, subset = DF.classifications == "Singlet")
  
  #  == Confirm the number of cells after filtering
  print("Number of cells after filtering Singlets:")
  print(dim(seurat_obj_singlets))
  
  #  == Visualize UMAP for singlets after filtering
  singlet_umap <- DimPlot(seurat_obj_singlets, reduction = "umap", label = TRUE, pt.size = 0.5) +
    labs(title = "UMAP Distribution of Singlets After Filtering")
  ggsave(file.path(doublet_dir, "UMAP_Singlets_After_Filtering.png"), plot = singlet_umap, width = 10, height = 6)
  
  #  == Save Seurat object after filtering singlets
  saveRDS(seurat_obj_singlets, file.path(doublet_dir, "seurat_obj_filtered_postSingletFiltering.rds"))
  
  return(seurat_obj_singlets)
}

#  == Balbc strain - Sham
Balbc_sham_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_sham_rep1_filtered,
  dataset_name = "Balbc_sham_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_sham_rep1_filtered, Balbc_sham_rep1_singlets)
gc()

Balbc_sham_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_sham_rep2_filtered,
  dataset_name = "Balbc_sham_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_sham_rep2_filtered, Balbc_sham_rep2_singlets)
gc()

#  == Balbc strain - Day 1
Balbc_day1_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day1_rep1_filtered,
  dataset_name = "Balbc_day1_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day1_rep1_filtered, Balbc_day1_rep1_singlets)
gc()

Balbc_day1_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day1_rep2_filtered,
  dataset_name = "Balbc_day1_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day1_rep2_filtered, Balbc_day1_rep2_singlets)
gc()

#  == Balbc strain - Day 3
Balbc_day3_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day3_rep1_filtered,
  dataset_name = "Balbc_day3_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day3_rep1_filtered, Balbc_day3_rep1_singlets)
gc()

Balbc_day3_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day3_rep2_filtered,
  dataset_name = "Balbc_day3_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day3_rep2_filtered, Balbc_day3_rep2_singlets)
gc()

#  == Balbc strain - Day 7
Balbc_day7_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day7_rep1_filtered,
  dataset_name = "Balbc_day7_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day7_rep1_filtered, Balbc_day7_rep1_singlets)
gc()

Balbc_day7_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Balbc_day7_rep2_filtered,
  dataset_name = "Balbc_day7_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Balbc_day7_rep2_filtered, Balbc_day7_rep2_singlets)
gc()

#  == Bl6 strain - Sham
Bl6_sham_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_sham_rep1_filtered,
  dataset_name = "Bl6_sham_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_sham_rep1_filtered, Bl6_sham_rep1_singlets)
gc()

Bl6_sham_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_sham_rep2_filtered,
  dataset_name = "Bl6_sham_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_sham_rep2_filtered, Bl6_sham_rep2_singlets)
gc()

#  == Bl6 strain - Day 1
Bl6_day1_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day1_rep1_filtered,
  dataset_name = "Bl6_day1_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day1_rep1_filtered, Bl6_day1_rep1_singlets)
gc()

Bl6_day1_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day1_rep2_filtered,
  dataset_name = "Bl6_day1_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day1_rep2_filtered, Bl6_day1_rep2_singlets)
gc()

#  == Bl6 strain - Day 3
Bl6_day3_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day3_rep1_filtered,
  dataset_name = "Bl6_day3_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day3_rep1_filtered, Bl6_day3_rep1_singlets)
gc()

Bl6_day3_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day3_rep2_filtered,
  dataset_name = "Bl6_day3_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day3_rep2_filtered, Bl6_day3_rep2_singlets)
gc()

#  == Bl6 strain - Day 7
Bl6_day7_rep1_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day7_rep1_filtered,
  dataset_name = "Bl6_day7_rep1",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day7_rep1_filtered, Bl6_day7_rep1_singlets)
gc()

Bl6_day7_rep2_singlets <- combined_preprocessing(
  seurat_obj_filtered = Bl6_day7_rep2_filtered,
  dataset_name = "Bl6_day7_rep2",
  nfeatures = 3000,
  ndims = 35,
  PCs = 1:20
)
rm(Bl6_day7_rep2_filtered, Bl6_day7_rep2_singlets)
gc()

# ====================================================================== Harmony Integration  ======================================================================

#  == Define Harmony Integration Directory
harmony_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony"
if (!dir.exists(harmony_dir)) {
  dir.create(harmony_dir, recursive = TRUE)
}

#  == Loadand List  Filtered Singlet Datasets
filtered_files <- list.files(
  path = "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/NSPCA_Doub",
  pattern = "seurat_obj_filtered_postSingletFiltering.rds",
  recursive = TRUE,
  full.names = TRUE
)

#  == Load datasets into a list
filtered_objects <- lapply(filtered_files, readRDS)

#  == Merge all filtered datasets into one Seurat object
merged_seurat <- merge(
  x = filtered_objects[[1]],
  y = filtered_objects[-1],
  add.cell.ids = gsub(".*/|\\..*", "", filtered_files),
  project = "Mouse_Ischemia_Merged"
)

#  == Save merged object for reference
saveRDS(merged_seurat, file.path(harmony_dir, "merged_seurat_before_harmony.rds"))

#  == Normalize and Find Variable Features
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 3000)

#  == Scale Data and Run PCA
merged_seurat <- ScaleData(merged_seurat, features = VariableFeatures(merged_seurat))
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(merged_seurat))

#  == Save PCA diagnostics
pca_dir <- file.path(harmony_dir, "PCA")
if (!dir.exists(pca_dir)) {
  dir.create(pca_dir)
}
elbow_plot <- ElbowPlot(merged_seurat, ndims = 50)
ggsave(file.path(pca_dir, "Elbow_Plot_Merged.png"), plot = elbow_plot, width = 8, height = 6)

#  == Run Harmony (metadata fields for batch correction)
merged_seurat <- RunHarmony(
  object = merged_seurat, 
  group.by.vars = c("strain", "condition", "replicate"), 
  dims.use = 1:20  # Updated to match NSPCA step
)

#  == Save Harmony-corrected object for reference
saveRDS(merged_seurat, file.path(harmony_dir, "merged_seurat_post_harmony.rds"))

#  == UMAP Visualization Pre-Harmony
merged_seurat <- RunUMAP(merged_seurat, reduction = "pca", dims = 1:20)  # Updated to use 20 PCs
umap_pre <- DimPlot(merged_seurat, group.by = "condition", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Before Harmony")
ggsave(file.path(harmony_dir, "UMAP_Before_Harmony.png"), plot = umap_pre, width = 10, height = 6)

#  == UMAP Visualization Post-Harmony
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:20)  # Updated to use 20 PCs
umap_post <- DimPlot(merged_seurat, group.by = "condition", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP After Harmony")
ggsave(file.path(harmony_dir, "UMAP_After_Harmony.png"), plot = umap_post, width = 10, height = 6)

#  == Find neighbors and clusters after Harmony integration
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:20)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.6)

#  == Save clustering results
clustering_dir <- file.path(harmony_dir, "Clustering")
if (!dir.exists(clustering_dir)) {
  dir.create(clustering_dir)
}

cluster_umap <- DimPlot(merged_seurat, group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP with Clustering (Post-Harmony)")
ggsave(file.path(clustering_dir, "UMAP_with_Clusters.png"), plot = cluster_umap, width = 10, height = 6)

#  == Save clustering results to CSV
cluster_table <- table(Idents(merged_seurat))
write.csv(as.data.frame(cluster_table), file.path(clustering_dir, "Cluster_Summary.csv"))

#  == Save the final Harmony-integrated and clustered Seurat object
final_seurat_file <- file.path(harmony_dir, "final_seurat_harmony_clustered.rds")
saveRDS(merged_seurat, final_seurat_file)

cat("Harmony integration and clustering complete. All outputs saved to:", harmony_dir, "\n")

# ====================================================================== Cell Typing Setup  ======================================================================

#  == Ensure 'final_seurat_harmony_clustered' is loaded
final_seurat_file <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony/final_seurat_harmony_clustered.rds"
if (!exists("final_seurat_harmony_clustered")) {
  cat("Loading Seurat object...\n")
  final_seurat_harmony_clustered <- readRDS(final_seurat_file)
}

#  == Renad's markers 
renad_markers <- c("PECAM1", "CD31" , "CDH5", "VEGFR2", "KDR", "ACTA2", "MYL9", "MYH11", "MCAM", "PDGFRB", "NG2", "CSPG4" ,
                   "TNMD", "SCX", "COL1A1", "MYH1", "MYH2", "MYH4", "ACTA1", "PPARG", "CD36",
                   "PDGFRA", "SCA1", "TGFB", "MMP14", "CXCL14", "S100A8", "S100A9", "CXCR2", "CCL3",
                   "CTSA", "CTSL" ,"CTSB", "CCR2", "CD14", "ITGAM", "LYZ2", "S100A4", "CCL2", "CD177",
                   "LY6C1", "LY6C2C", "IRF5", "CLL2", "CCL6", "CCL9", "TNF", "IL10", "NOS2",
                   "SPP1", "IL17RA", "CXCL3", "CXCL16", "CCR2", "CXCL16", "C1QC", "IL6RA", "CSF1R",
                   "AIF1", "CD14", "SLFN4", "CXCL2", "FOSL2", "ATF4", "MDM2", "CX3CR1", "LY86",
                   "LPXN", "DHRS3", "DNMT3A", "ARG1", "CD163", "MRC1", "CCL8", "FOLR2", "APOE",
                   "TREM2", "CD3E", "CD4", "CCL5", "CD8A", "NKG7", "KLRA7", "NCAM1", "CD56", "CD19",
                   "MS4A1", "CD20", "CD1C", "HLA-DR", "HLADR", "CD11C", "ITGAX", "CHODL", "HEYL", "PAX7", "SDC3", "NOTCH2/3",
                   "DPT", "CALCR", "COL15A1", "SPRY1", "MYOD1", "SDC4", "TCF7", "CARM1", "MAPK",
                   "FOS", "MEST", "CTSB", "IFIT1", "IFIT3", "ISG15", "CXCL", "KI67", "RAC1", "TOP2A", "EZH2",
                   "BIRC5", "CDK1", "ACTA1","MYOG", "TTN", "MYH", "MYL4", "NOTCH2","NOTCH3","HLADR")

#  == Convert markers to title case to match Seurat object rownames
renad_markers <- stringr::str_to_title(renad_markers)

present_markers <- c()
missing_markers <- c()

#  == Check each marker
for (marker in renad_markers) {
  if (marker %in% rownames(final_seurat_harmony_clustered)) {
    present_markers <- c(present_markers, marker)
  } else {
    missing_markers <- c(missing_markers, marker)
  }
}

# Print results
cat("Present markers:\n", paste(present_markers, collapse = ", "), "\n\n")
cat("Missing markers:\n", paste(missing_markers, collapse = ", "), "\n")

# Original Renad Markers 
Renad_markers <- list(
  "Endothelial" = c("PECAM1", "CD31", "CDH5", "VEGFR2", "KDR"),
  "Smooth Muscle" = c("ACTA2", "MYL9", "MYH11"),
  "Pericytes" = c("MCAM", "PDGFRB", "NG2", "CSPG4"),
  "Tenocytes" = c("TNMD", "SCX", "COL1A1"),
  "Mature" = c("MYH1", "MYH2", "MYH4", "ACTA1"),
  "FAPs COmbined" = c("PPARG", "CD36", "TGFB", "MMP14", "PDGFRA", "SCA1", "CXCL14"),
  "FAPs (Adipogenic)" = c("PPARG", "CD36"),
  "FAPs (Stem)" = c("PDGFRA", "SCA1"),
  "FAPs (Pro-remodeling)" = c("TGFB", "MMP14"),
  "FAPs (Cxc14+)" = c("CXCL14"),
  "Neutrophils" = c("S100A8", "S100A9", "CXCR2", "CCL3"),
  "M0M Monocytes" = c("CTSA", "CTSL", "CTSB", "CCR2", "CD14", "ITGAM", "LYZ2", "S100A4", "CCL2"),
  "Ly6C+ Monocytes" = c("CD177", "LY6C1", "LY6C2C", "IRF5"),
  "Monocytes Combined" = c("CTSA", "CTSL", "CTSB", "CCR2", "CD14", "ITGAM", "LYZ2", "S100A4", "CCL2","CD177", "LY6C1", "LY6C2C", "IRF5"),
  "M1 Macrophages COmbined" = c("CLL2", "CCL6", "CCL9", "TNF", "IL10", "NOS2", "SPP1", "IL17RA", "CXCL3", "CXCL16", "CCR2"),
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
  "QSCs" = c("CHODL", "HEYL", "PAX7", "SDC3"),
  "Self-Renewing QSCs" = c("NOTCH2", "NOTCH3", "DPT", "CALCR", "COL15A1", "SPRY1", "NOTCH2/3"),
  "QSCs" = c("CHODL", "HEYL", "PAX7", "SDC3", "NOTCH2", "NOTCH3", "DPT", "CALCR", "COL15A1", "SPRY1", "NOTCH2/3"),
  "ASCs (early activated)" = c("MYOD1", "SDC4", "TCF7", "CARM1", "MAPK", "FOS", "MEST"),
  "IMB" = c("CTSB", "IFIT1", "IFIT3", "ISG15", "CXCL"),
  "Proliferating ASCs" = c("KI67", "RAC1", "TOP2A", "EZH2", "BIRC5", "CDK1"),
  "Myocytes & Differentiating" = c("ACTA1", "MYOG", "TTN", "MYH", "MYL4")
)

#  == Convert each marker list to title case
Renad_markers <- lapply(Renad_markers, stringr::str_to_title)

# Define directory for saving The OG Renad Markers plots 
renad_markers_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Renad Markers"

# Ensure the directory exists
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
for (cell_type in names(Renad_markers)) {
  print(paste("Generating plots for:", cell_type))
  
  # Generate and save plots
  tryCatch({
    plots <- generate_combined_plots(
      seurat_object = final_seurat_harmony_clustered,
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

# ======================================================================Identify Top Markers and Generate Plots  ======================================================================

#  ==  parameters for dynamic thresholds and top marker selection
min_pct <- 0.25  # Minimum percentage of cells expressing the marker
logfc_threshold <- 0.25  # Minimum log-fold change threshold
top_n_markers <- 9  # Number of top markers per cluster

#  ==  directories for saving outputs
top_markers_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Top Markers"
plots_dir <- file.path(top_markers_dir, "Plots")
heatmap_dir <- file.path(top_markers_dir, "Heatmaps")

#  == make sure directories exist
if (!dir.exists(top_markers_dir)) dir.create(top_markers_dir, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)

#  == Identify markers for each cluster
cluster_markers <- FindAllMarkers(
  object = final_seurat_harmony_clustered,
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

#  == Save the complete marker list
write.csv(cluster_markers, file.path(top_markers_dir, "Cluster_Markers.csv"))

#  == Step 2: Extract top markers per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = top_n_markers, order_by = avg_log2FC)

#  == Save the top markers
write.csv(top_markers, file.path(top_markers_dir, "Top_Cluster_Markers.csv"))

#  == Generate plots for each cluster
for (cluster_id in unique(top_markers$cluster)) {
  print(paste("Processing Cluster:", cluster_id))
  
  #  == Subset top markers for this cluster
  cluster_top_markers <- top_markers %>% filter(cluster == cluster_id)
  genes <- cluster_top_markers$gene  # Extract top marker genes
  
  #  == Generate Combined Violin Plot
  vln_plot <- VlnPlot(
    object = final_seurat_harmony_clustered,
    features = genes,
    group.by = "seurat_clusters",
    pt.size = 0.5,
    ncol = min(3, length(genes))
  ) +
    ggtitle(paste("Violin Plot for Top Markers (Cluster", cluster_id, ")")) +
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  #  == Save Plot
  ggsave(
    filename = file.path(plots_dir, paste0("Cluster_", cluster_id, "_ViolinPlot.png")),
    plot = vln_plot,
    width = 15,
    height = 10
  )
  
  #  == Generate Combined Feature Plot
  feature_plot <- FeaturePlot(
    object = final_seurat_harmony_clustered,
    features = genes,
    reduction = "umap",
    cols = c("lightgrey", "blue"),
    ncol = min(3, length(genes))
  ) +
    ggtitle(paste("Feature Plot for Top Markers (Cluster", cluster_id, ")")) +
    theme(
      text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  #  == Save Plot
  ggsave(
    filename = file.path(plots_dir, paste0("Cluster_", cluster_id, "_FeaturePlot.png")),
    plot = feature_plot,
    width = 15,
    height = 10
  )
}

################################################################## Enrichment Scores Pipeline ##################################################################

#  == Define output directories
output_dir <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Cell Typing/Broad_Annotations"
enrichment_dir <- file.path(output_dir, "Enrichment_Scores")
plots_dir <- file.path(output_dir, "Validation_Plots")

#  == Ensure directories exist
if (!dir.exists(enrichment_dir)) dir.create(enrichment_dir, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

#  == Ensure 'seurat_obj_harmony' is loaded
final_seurat_file <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony/final_seurat_harmony_clustered.rds"
if (!exists("seurat_obj_harmony")) {
  cat("Loading Seurat object...\n")
  seurat_obj_harmony <- readRDS(final_seurat_file)
}

#  == Check if Seurat object exists and is valid
if (!inherits(seurat_obj_harmony, "Seurat")) {
  stop("The loaded object is not a valid Seurat object")
}

#  == Add module scores for each marker set
cat("Calculating enrichment scores for each cell type...\n")
for (cell_type in names(Renad_markers)) {
  markers <- Renad_markers[[cell_type]]  # Get marker list for the cell type
  
  #  == Add module scores for the marker set
  seurat_obj_harmony <- AddModuleScore(
    object = seurat_obj_harmony,
    features = list(markers),
    name = cell_type
  )
}

# ====================================================================== Visualize Enrichment Scores ======================================================================

#  == Aggregate scores by cluster
cat("Aggregating enrichment scores by cluster...\n")
enrichment_scores <- data.frame(
  cluster = unique(seurat_obj_harmony$seurat_clusters)
)

for (cell_type in names(Renad_markers)) {
  score_col <- paste0(cell_type, "1")  e
  scores <- aggregate(
    seurat_obj_harmony[[score_col]],
    by = list(cluster = seurat_obj_harmony$seurat_clusters),
    FUN = mean
  )
  colnames(scores)[2] <- cell_type
  enrichment_scores <- merge(enrichment_scores, scores, by = "cluster")
}

#  == Save enrichment scores 
scores_file <- file.path(enrichment_dir, "Enrichment_Scores_By_Cluster.csv")
write.csv(enrichment_scores, file = scores_file, row.names = FALSE)
cat("Enrichment scores saved to:", scores_file, "\n")

#  == Convert enrichment scores to long format
enrichment_scores_long <- enrichment_scores %>%
  pivot_longer(
    cols = -cluster,  # All columns except "cluster"
    names_to = "cell_type",  # New column for cell type names
    values_to = "enrichment_score"  # New column for enrichment scores
  )

#  == Generate dot plot 
dot_plot <- ggplot(enrichment_scores_long, aes(x = cell_type, y = factor(cluster))) +
  geom_point(aes(size = enrichment_score, color = enrichment_score)) +
  scale_color_viridis_c(option = "plasma", name = "Enrichment Score") +
  scale_size(range = c(2, 10), name = "Score Magnitude") +
  labs(
    title = "Dot Plot of Enrichment Scores by Cluster",
    x = "Cell Type",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

#  == Save dot plot to file
dot_plot_file <- file.path(enrichment_dir, "Enrichment_Scores_DotPlot.png")
ggsave(filename = dot_plot_file, plot = dot_plot, width = 15, height = 8)
cat("Dot plot saved to:", dot_plot_file, "\n")

#  == Define a new directory for saving individual dot plots
individual_dot_plots_dir <- file.path(enrichment_dir, "Individual_Dot_Plots")
if (!dir.exists(individual_dot_plots_dir)) dir.create(individual_dot_plots_dir, recursive = TRUE)

#  == Loop through each marker set and generate individual dot plots
for (cell_type in names(Renad_markers)) {
  # Enrichment score column for the current marker set
  score_col <- paste0(cell_type, "1")  
  
  #  == Check if the score column exists in the metadata
  if (!(score_col %in% colnames(seurat_obj_harmony@meta.data))) {
    cat("Warning: Enrichment score column", score_col, "not found for", cell_type, "\n")
    next
  }
  
  #  == Extract enrichment scores and cluster information
  enrichment_data <- data.frame(
    cluster = seurat_obj_harmony$seurat_clusters,
    enrichment_score = seurat_obj_harmony@meta.data[[score_col]]
  )
  
  #  == Ensure enrichment_score is numeric
  if (!is.numeric(enrichment_data$enrichment_score)) {
    cat("Warning: Enrichment scores for", cell_type, "are not numeric. Skipping.\n")
    next
  }
  
  #  == Calculate the average enrichment score for the marker set
  avg_score <- round(mean(enrichment_data$enrichment_score, na.rm = TRUE), 2)
  
  #  == Generate the dot plot
  individual_dot_plot <- ggplot(enrichment_data, aes(x = cell_type, y = factor(cluster))) +
    geom_point(aes(size = enrichment_score, color = enrichment_score)) +
    scale_color_viridis_c(option = "magma", name = "Enrichment Score") +
    scale_size(range = c(2, 10), name = "Score Magnitude") +
    labs(
      title = paste("Dot Plot for", cell_type),
      subtitle = paste("Average Enrichment Score:", avg_score),
      x = "Marker Set",
      y = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  #  == Save the individual dot plot
  individual_plot_file <- file.path(individual_dot_plots_dir, paste0(cell_type, "_DotPlot.png"))
  ggsave(filename = individual_plot_file, plot = individual_dot_plot, width = 10, height = 8)
  
  #  == Inform about the saved plot
  cat("Individual dot plot saved for", cell_type, "to:", individual_plot_file, "\n")
}


# ====================================================================== Assign Initial Labels then compare with Renad Plots ======================================================================

cat("Assigning initial broad labels to clusters...\n")
initial_annotations <- enrichment_scores %>%
  tidyr::pivot_longer(-cluster, names_to = "cell_type", values_to = "score") %>%
  group_by(cluster) %>%
  slice_max(score, n = 1) %>%  # Select the cell type with the highest score
  select(cluster, cell_type)

#  == Save initial annotations to a CSV
annotations_file <- file.path(output_dir, "Initial_Broad_Annotations.csv")
write.csv(initial_annotations, file = annotations_file, row.names = FALSE)
cat("Initial annotations saved to:", annotations_file, "\n")

#  == Add annotations to Seurat metadata
seurat_obj_harmony$Broad_Annotations <- sapply(
  seurat_obj_harmony$seurat_clusters,
  function(cluster) {
    annotation <- initial_annotations$cell_type[initial_annotations$cluster == cluster]
    if (length(annotation) > 0) return(annotation) else return("Unassigned")
  }
)

# =========================== Load and Annotate the dataset following inspection ===============================

cat("Loading the Harmony-processed Seurat object...\n")
final_harmony_file <- "D:/Masters/Processing Codes/Mouse_GSE227075_Initial_Processing_Plots/Harmony/Final_Seurat_Harmony_Clustered.rds"
seurat_obj <- readRDS(final_harmony_file)
cat("Successfully loaded the Seurat object.\n") 

#  == Define manual annotations for the clusters directly with updated labels
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

#  == Save the Seurat object
annotated_rds_path <- file.path(seurat_dir, "Southerland_Overall_Seurat_Annotated.rds")
saveRDS(seurat_obj, annotated_rds_path)
cat("Annotated Seurat object saved to:", annotated_rds_path, "\n")

#  == Generate the New UMAP
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

#  == Save UMAP plot (PDF and PNG)
umap_pdf_path <- file.path(overall_plots_dir, "UMAP_Manual_Annotation.pdf")
umap_png_path <- file.path(overall_plots_dir, "UMAP_Manual_Annotation.png")

ggsave(filename = umap_pdf_path, plot = umap_plot, width = 8, height = 6)
ggsave(filename = umap_png_path, plot = umap_plot, width = 8, height = 6, dpi = 300)

cat("UMAP plot saved as PDF and PNG in:", overall_plots_dir, "\n")
