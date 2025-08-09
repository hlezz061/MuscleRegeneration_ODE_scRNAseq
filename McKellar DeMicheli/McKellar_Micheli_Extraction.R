# =====================================================
# Load and inspect scMuscle_mm10_slim_v1-1.RData 
# to identify and extract McKellar and De Micheli data 
# =====================================================

# ─────────────────────────────────────────────────────────────────────────────
# Load libraries
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(crayon)
  library(dplyr)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# Define paths and create directory structure
# ─────────────────────────────────────────────────────────────────────────────
ref_path <- "D:/Masters/Renad_Annotation_V2/Data/scMuscle_mm10_slim_v1-1.RData"
base_dir <- "D:/Masters/Processing Codes/Cosgrove_Annotation"

log_dir      <- file.path(base_dir, "Logs")
meta_dir     <- file.path(base_dir, "Metadata")
structure_dir <- file.path(base_dir, "Structure")

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE)
dir.create(meta_dir, showWarnings = FALSE)
dir.create(structure_dir, showWarnings = FALSE)

cat(blue$bold("\n Directory structure initialized.\n"))

# ─────────────────────────────────────────────────────────────────────────────
# Load the reference atlas
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold(" Loading reference Seurat object...\n"))
load(ref_path)  
seurat_vars <- ls()
ref_obj_name <- seurat_vars[grepl("seurat", seurat_vars, ignore.case = TRUE)]

if (length(ref_obj_name) != 1) {
  stop(red$bold("Could not uniquely identify Seurat object in the .RData file."))
}
ref <- get(ref_obj_name)
rm(list = setdiff(seurat_vars, ref_obj_name))  # Clean up unused variables

cat(green$bold(paste0(" Loaded object: ", ref_obj_name, "\n")))
cat(green$bold(paste0(" Total cells: ", ncol(ref), " | Genes: ", nrow(ref), "\n")))

# ─────────────────────────────────────────────────────────────────────────────
# Save basic structure info
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Saving object structure summary...\n"))

sink(file.path(structure_dir, "object_structure.txt"))
str(ref, max.level = 2)
sink()

assays_used <- Assays(ref)
reductions_used <- Reductions(ref)
layers_present <- Layers(ref[["RNA"]])

cat(cyan$bold("Assays:"), paste(assays_used, collapse = ", "), "\n")
cat(cyan$bold("Reductions:"), paste(reductions_used, collapse = ", "), "\n")
cat(cyan$bold("RNA Layers:"), paste(layers_present, collapse = ", "), "\n")

write_lines(assays_used, file.path(structure_dir, "assays.txt"))
write_lines(reductions_used, file.path(structure_dir, "reductions.txt"))
write_lines(layers_present, file.path(structure_dir, "RNA_layers.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# Inspect and save metadata
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Inspecting metadata columns...\n"))

meta_df <- ref@meta.data
meta_cols <- colnames(meta_df)
cat(cyan$bold("Metadata columns:"), paste(meta_cols, collapse = ", "), "\n")

write_lines(meta_cols, file.path(meta_dir, "metadata_columns.txt"))
write_csv(meta_df, file.path(meta_dir, "metadata_full.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# Preview key metadata fields (if present or inferred)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Previewing key metadata fields...\n"))

#  == Manually define likely fields based on actual column names
preview_cols <- c("sample", "injury", "injury.days", "seurat_clusters", 
                  "harmony_factorIDs", "harmony_res.1.2_IDs", "bbknn_factorIDs", 
                  "DF.individual", "DF.pANN.individual")

found_cols <- preview_cols[preview_cols %in% meta_cols]

if (length(found_cols) == 0) {
  cat(red$bold(" No matching metadata fields found for preview.\n"))
} else {
  summary_file <- file.path(meta_dir, "metadata_previews.txt")
  sink(summary_file)
  
  for (col in found_cols) {
    cat(blue$bold(paste0("\n>> ", col, " (", class(meta_df[[col]]), "):\n")))
    print(table(meta_df[[col]]))
  }
  
  sink()
  cat(green$bold(paste0("Saved preview summary to: ", summary_file, "\n")))
}

#  == Save separate CSV of just these preview columns
if (length(found_cols) > 0) {
  preview_df <- meta_df[, found_cols, drop = FALSE]
  write_csv(preview_df, file.path(meta_dir, "metadata_preview_columns.csv"))
  cat(green$bold(" Saved metadata preview columns as CSV.\n"))
}

cat(green$bold("\n Inspection complete. All results saved under:\n"))
cat(blue$bold(paste0("→ ", base_dir, "\n\n")))

# ─────────────────────────────────────────────────────────────────────────────
# Inspect origin of cells (source studies and dataset labels)
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Inspecting dataset origin: `source` and `source.label`...\n"))

origin_fields <- c("sample", "source", "source.label")
origin_cols <- origin_fields[origin_fields %in% meta_cols]

if (length(origin_cols) == 0) {
  cat(red$bold("No dataset source columns found.\n"))
} else {
  #  == Save full table of sample → source mapping
  origin_df <- meta_df[, origin_cols, drop = FALSE]
  origin_summary <- origin_df %>%
    group_by(sample, source, source.label) %>%
    summarise(NumCells = n(), .groups = "drop") %>%
    arrange(desc(NumCells))
  
  #  == Save CSV
  write_csv(origin_summary, file.path(meta_dir, "sample_to_source_mapping.csv"))
  cat(green$bold(" Saved `sample` → `source`/`source.label` mapping to CSV.\n"))
  
  #  == Print a sample of the mapping
  cat(cyan$bold("\nSample of mapping:\n"))
  print(head(origin_summary, 10))
}

# ─────────────────────────────────────────────────────────────────────────────
# Filter McKellar and De Micheli cells and inspect
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Inspecting McKellar and De Micheli subsets...\n"))

#  == Define source labels of interest
mckellar_label   <- "McKellar/Walter 2020"
demicheli_label  <- "De Micheli 2020a"

#  == Subset metadata
mckellar_cells   <- rownames(meta_df[meta_df$source.label == mckellar_label, ])
demicheli_cells  <- rownames(meta_df[meta_df$source.label == demicheli_label, ])

#  == Summary of cell counts
cat(green$bold(paste0(" McKellar cells:  ", length(mckellar_cells), "\n")))
cat(green$bold(paste0(" De Micheli cells:", length(demicheli_cells), "\n\n")))

#  == Unique samples in each
cat(cyan$bold(" Unique McKellar samples:\n"))
print(unique(meta_df[mckellar_cells, "sample"]))

cat(cyan$bold("\n Unique De Micheli samples:\n"))
print(unique(meta_df[demicheli_cells, "sample"]))

#  == Unique annotations in each
cat(cyan$bold("\n️ Cell type annotations in McKellar (harmony_factorIDs):\n"))
print(table(meta_df[mckellar_cells, "harmony_factorIDs"]))

cat(cyan$bold("\n️ Cell type annotations in De Micheli (harmony_factorIDs):\n"))
print(table(meta_df[demicheli_cells, "harmony_factorIDs"]))

# ─────────────────────────────────────────────────────────────────────────────
# SUBSET + SAVE McKellar and De Micheli Seurat Objects
# ─────────────────────────────────────────────────────────────────────────────
cat(blue$bold("\n Saving Seurat objects for McKellar and De Micheli datasets...\n"))

#  == Define output directory
seurat_dir <- file.path(base_dir, "Seurat_Object")
dir.create(seurat_dir, showWarnings = FALSE)

#  == Subset full reference object
ref_mckellar  <- subset(ref, cells = mckellar_cells)
ref_demicheli <- subset(ref, cells = demicheli_cells)

#  == Confirm object dimensions
cat(green$bold(paste0("McKellar object: ", ncol(ref_mckellar), " cells\n")))
cat(green$bold(paste0("De Micheli object: ", ncol(ref_demicheli), " cells\n")))

#  == Save objects
saveRDS(ref_mckellar,  file.path(seurat_dir, "ref_McKellar.rds"))
saveRDS(ref_demicheli, file.path(seurat_dir, "ref_DeMicheli.rds"))

cat(green$bold("\n Saved Seurat objects to:\n"))
cat(cyan$bold(paste0("→ ", file.path(seurat_dir, "ref_McKellar.rds"), "\n")))
cat(cyan$bold(paste0("→ ", file.path(seurat_dir, "ref_DeMicheli.rds"), "\n")))
