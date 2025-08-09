# ================================================
# Code to Calculate the Southerland Annotations 
# True proportions using new data 
# Overall comparison with the model 
# ================================================

# ======================================================================Load any potential Libraries ====================================================================== 

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

# ====================================================================== Define File Paths ======================================================================

# == Define file paths
model_file   <- "D:/Masters/Processing Codes/Model Sim/Mouse_GSE227075_Sim 2/proportions/Trajectory_Percentages.csv"
main_file    <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/Proportions/Main Annotation/Data/Main_Cell_Type_Percentages_Collapsed_Overall.csv"
myeloid_file <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/Proportions/Myeloid Annotation/Data/Cell_Type_Percentages_Collapsed_Myeloid.csv"

# == Read CSV files
model_df   <- read.csv(model_file)
main_df    <- read.csv(main_file)
myeloid_df <- read.csv(myeloid_file)

# == Quick look at the data structures
head(model_df)    # Model output (time vs cell type percentages)
head(main_df)     # Main annotated cell types (proportion by time)
head(myeloid_df)  # Myeloid subclusters (proportion by time)

colnames(main_df)
colnames(myeloid_df)

# == Rename column for clarity before conversion
main_df <- dplyr::rename(main_df, RawProportion = percentage)
myeloid_df <- dplyr::rename(myeloid_df, RawProportion = percentage)

colnames(main_df)
colnames(myeloid_df)

head(main_df)
head(myeloid_df)

# ====================================================================== Standardize Time Labels ======================================================================

# == Convert character time labels to numeric days
main_df <- main_df %>% mutate(
  time_point = as.character(time_point),
  Days = case_when(
    str_to_lower(time_point) %in% c("sham", "d0", "0") ~ 0,
    grepl("day ?0$", str_to_lower(time_point)) ~ 0,
    TRUE ~ as.numeric(str_extract(time_point, "\\d+"))
  )
)

myeloid_df <- myeloid_df %>% mutate(
  time_point = as.character(time_point),
  Days = case_when(
    str_to_lower(time_point) %in% c("sham", "d0", "0") ~ 0,
    TRUE ~ as.numeric(str_extract(time_point, "\\d+"))
  )
)


# == Check converted timepoints and model time range
unique(main_df$Days)
unique(myeloid_df$Days)
range(model_df$time)

# ====================================================================== Prepare Model Outputs (Main Data Categories) ======================================================================

model_main <- model_df %>%
  filter(time %in% main_df$Days) %>%
  transmute(
    Days = time,
    Neutrophils        = if ("N"   %in% names(.)) N   else 0,
    Monocytes          = if ("M"   %in% names(.)) M   else 0,
    `M1 Macrophages`   = if ("M1"  %in% names(.)) M1  else 0,
    `M2 Macrophages`   = if ("M2"  %in% names(.)) M2  else 0,
    QSCs               = if ("QSC" %in% names(.)) QSC else 0,
    ASCs               = if ("ASC" %in% names(.)) ASC else 0,
    Myocytes           = if ("Mc"  %in% names(.)) Mc  else 0
  )

# ====================================================================== Prepare Main Data ======================================================================

modeled_categories <- c("Neutrophils", "Myeloid Cells", "QSCs", "ASCs", "Myocytes")
main_df$CellType <- main_df$cell_annotation

# == Print available unique cell types in main data (before filtering)
cat("\nUnique cell types in main_df before filtering:\n")
print(unique(main_df$CellType))

# == Combine ASCs and Prlf.ASCs under one label
main_df$CellType <- main_df$cell_annotation
main_df$CellType <- ifelse(main_df$CellType == "Prlf.ASCs", "ASCs", main_df$CellType)

main_df <- main_df %>%
  group_by(time_point, Days, CellType) %>%
  summarise(
    RawProportion = sum(RawProportion),
    .groups = "drop"
  )

# == Print available unique cell types in main data (after merging)
cat("\nUnique cell types in main_df after merging ASCs:\n")
print(unique(main_df$CellType))

# == Filter only modeled categories and recalculate proportions
main_filtered_raw <- main_df %>%
  filter(CellType %in% modeled_categories)

main_filtered <- main_filtered_raw %>%
  group_by(Days) %>%
  mutate(Percentage = 100 * RawProportion / sum(RawProportion)) %>%
  ungroup()

# == Print filtered result
cat("\nPreview of main_filtered (after filtering to modeled categories):\n")
print(head(main_filtered))

combined_main <- main_filtered %>%
  select(Days, CellType, DataPercent = Percentage) %>%
  inner_join(
    model_main %>%
      pivot_longer(cols = -Days, names_to = "CellType", values_to = "ModelPercent"),
    by = c("Days", "CellType")
  )

combined_main <- combined_main %>%
  filter(CellType != "Myeloid Cells")

# == Show merged results
cat("\nCombined model vs data (main categories):\n")
print(head(combined_main))

# == Sanity check: how many rows merged
cat("\nTotal merged rows in combined_main:", nrow(combined_main), "\n")

# == Save combined_main to CSV
output_path <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Data/Combined_Main_Model_Comparison.csv"
write.csv(combined_main, output_path, row.names = FALSE)

cat("\n CSV saved to:\n", output_path, "\n")

# ====================================================================== Adjust Myeloid Subcluster Percentages (Scale to Whole Dataset) ======================================================================

cat("\n==== Adjusting Myeloid Subclusters Relative to Whole Dataset ====\n")

# == Get total percentage of 'Myeloid Cells' (already adjusted based on modeled categories)
myeloid_total_df <- main_filtered %>%
  filter(CellType == "Myeloid Cells") %>%
  select(Days, MyeloidTotal = Percentage)  # <-- Already scaled properly earlier

cat("\nAdjusted Myeloid totals from main_filtered (used for scaling):\n")
print(myeloid_total_df)

# == Scale subclusters to whole dataset
myeloid_scaled <- myeloid_df %>%
  left_join(myeloid_total_df, by = "Days") %>%
  mutate(
    ScaledPercent = (RawProportion / 100) * MyeloidTotal
  )

cat("\nPreview of scaled myeloid subclusters:\n")
print(head(myeloid_scaled))

# == Keep only modeled myeloid types
myeloid_scaled$CellType <- myeloid_scaled$final_annotation
modeled_myeloid <- c("Neutrophils", "Monocytes", "M1 Macrophages", "M2 Macrophages")

myeloid_filtered <- myeloid_scaled %>%
  filter(CellType %in% modeled_myeloid)

cat("\nFiltered scaled myeloid subclusters (only modeled types):\n")
print(unique(myeloid_filtered$CellType))

# == Prepare model outputs for these myeloid types
model_myeloid <- model_df %>%
  filter(time %in% myeloid_df$Days) %>%
  transmute(
    Days = time,
    Neutrophils        = if ("N"   %in% names(.)) N   else 0,
    Monocytes          = if ("M"   %in% names(.)) M   else 0,
    `M1 Macrophages`   = if ("M1"  %in% names(.)) M1  else 0,
    `M2 Macrophages`   = if ("M2"  %in% names(.)) M2  else 0
  )

cat("\nPreview of model_myeloid (for matching days):\n")
print(head(model_myeloid))

# == Merge scaled empirical data with model predictions
combined_myeloid <- myeloid_filtered %>%
  select(Days, CellType, DataPercent = ScaledPercent) %>%
  inner_join(
    model_myeloid %>%
      pivot_longer(cols = -Days, names_to = "CellType", values_to = "ModelPercent"),
    by = c("Days", "CellType")
  )

# == Final output
cat("\nCombined model vs data (myeloid categories):\n")
print(head(combined_myeloid))
cat("\nTotal merged rows in combined_myeloid:", nrow(combined_myeloid), "\n")

output_path_myeloid <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Data/Combined_Myeloid_Model_Comparison.csv"
write.csv(combined_myeloid, output_path_myeloid, row.names = FALSE)
cat("\nCSV saved to:\n", output_path_myeloid, "\n")

# ====================================================================== Combine Neutrophils and Create Final Dataset ======================================================================

cat("\n==== Combining Neutrophils and Preparing Final Output ====\n")

# == Separate neutrophils from both dataframes
neutrophils_main <- combined_main %>% filter(CellType == "Neutrophils")
neutrophils_myeloid <- combined_myeloid %>% filter(CellType == "Neutrophils")

# == Combine neutrophils (DataPercent only; use model from myeloid since it's more specific)
neutrophils_combined <- neutrophils_main %>%
  select(Days, DataPercent_main = DataPercent) %>%
  left_join(
    neutrophils_myeloid %>% select(Days, DataPercent_myeloid = DataPercent, ModelPercent),
    by = "Days"
  ) %>%
  mutate(
    DataPercent = DataPercent_main + DataPercent_myeloid,
    CellType = "Neutrophils"
  ) %>%
  select(Days, CellType, DataPercent, ModelPercent)

cat("\nCombined Neutrophil proportions (main + myeloid):\n")
print(neutrophils_combined)

# == Remove original neutrophils from main and myeloid
combined_main_clean <- combined_main %>% filter(CellType != "Neutrophils")
combined_myeloid_clean <- combined_myeloid %>% filter(CellType != "Neutrophils")

# == Combine everything together
combined_all <- bind_rows(
  combined_main_clean,
  combined_myeloid_clean,
  neutrophils_combined
) %>%
  arrange(Days, CellType)

cat("\nFinal combined dataset (7 modeled cell types):\n")
print(combined_all)

# == Save to CSV
final_output_path <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Data/True_Model_Data_Comparison.csv"
write.csv(combined_all, final_output_path, row.names = FALSE)

cat("\nFinal merged true proportions saved to:\n", final_output_path, "\n")

# == Rescale DataPercent so each day's modeled cell types sum to 100%
combined_all <- combined_all %>%
  group_by(Days) %>%
  mutate(
    RescaledDataPercent = 100 * DataPercent / sum(DataPercent)
  ) %>%
  ungroup()

cat("\nFinal rescaled dataset (sums to 100% per day):\n")
print(head(combined_all))

# == Overwrite or save to new file
rescaled_output_path <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Data/True_Model_Data_Comparison_Rescaled.csv"
write.csv(combined_all, rescaled_output_path, row.names = FALSE)

cat("\nRescaled true proportions saved to:\n", rescaled_output_path, "\n")

# ====================================================================== Pearson Correlation: Model vs Rescaled Data ======================================================================

# == Output directory
plots_dir <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Plots"

# == Compute correlation per cell type
cor_df <- combined_all %>%
  group_by(CellType) %>%
  summarise(
    Pearson = cor(RescaledDataPercent, ModelPercent, method = "pearson"),
    .groups = "drop"
  )

# == Preview correlations
cat("\nPearson correlation (Rescaled vs Model) per CellType:\n")
print(cor_df)

# == Save as CSV
cor_output_path <- file.path(plots_dir, "Pearson_Correlations_Rescaled_vs_Model.csv")
write.csv(cor_df, cor_output_path, row.names = FALSE)
cat("\nSaved Pearson correlation CSV to:\n", cor_output_path, "\n")

# == Bar plot of Pearson correlations
cor_plot <- ggplot(cor_df, aes(x = reorder(CellType, Pearson), y = Pearson, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(Pearson, 3)), vjust = -0.6, size = 3.5) +
  labs(
    title = "Pearson Correlation: Rescaled Data vs Model",
    x = "Cell Type", y = "Pearson r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  ylim(0, 1)

# == Save plot
cor_plot_path <- file.path(plots_dir, "Pearson_Correlations_Rescaled_vs_Model_Barplot.png")
ggsave(cor_plot_path, plot = cor_plot, width = 7, height = 5, dpi = 300)

print(cor_plot)
cat("\nSaved Pearson correlation bar plot to:\n", cor_plot_path, "\n")

# ====================================================================== Pearson Correlation: By Timepoint ======================================================================

# == Compute correlation per timepoint
cor_by_day_df <- combined_all %>%
  group_by(Days) %>%
  summarise(
    Pearson = cor(RescaledDataPercent, ModelPercent, method = "pearson"),
    .groups = "drop"
  )

# == Preview results
cat("\nPearson correlation (Rescaled vs Model) per Timepoint:\n")
print(cor_by_day_df)

# == Save to CSV
cor_by_day_output_path <- file.path(plots_dir, "Pearson_Correlations_By_Timepoint.csv")
write.csv(cor_by_day_df, cor_by_day_output_path, row.names = FALSE)
cat("\nSaved Pearson correlation per timepoint to:\n", cor_by_day_output_path, "\n")

# ====================================================================== Statistical Comparison: Scaled Data ======================================================================

# == Define RMSE and MAE
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
mae  <- function(x, y) mean(abs(x - y), na.rm = TRUE)

# == Output directory
stats_dir <- file.path(plots_dir, "stat_plots2")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

# == Safe correlation
safe_cor <- function(x, y) {
  if (length(x) < 2 || length(y) < 2 || all(x == 0) || all(y == 0)) return(NA)
  cor(x, y, method = "pearson")
}

# == By Cell Type
metrics_scaled_byCT <- combined_all %>%
  group_by(CellType) %>%
  summarise(
    RMSE = rmse(RescaledDataPercent, ModelPercent),
    MAE = mae(RescaledDataPercent, ModelPercent),
    Pearson_r = safe_cor(RescaledDataPercent, ModelPercent),
    .groups = "drop"
  )

# == Overall
metrics_scaled_overall <- combined_all %>%
  summarise(
    RMSE = rmse(RescaledDataPercent, ModelPercent),
    MAE = mae(RescaledDataPercent, ModelPercent),
    Pearson_r = safe_cor(RescaledDataPercent, ModelPercent)
  ) %>%
  mutate(CellType = "Overall") %>%
  select(CellType, RMSE, MAE, Pearson_r)

# == By Timepoint
metrics_scaled_byTime <- combined_all %>%
  group_by(Days) %>%
  summarise(
    RMSE = rmse(RescaledDataPercent, ModelPercent),
    MAE = mae(RescaledDataPercent, ModelPercent),
    .groups = "drop"
  )

write.csv(bind_rows(metrics_scaled_byCT, metrics_scaled_overall),
          file.path(stats_dir, "Metrics_byCellType_Rescaled.csv"), row.names = FALSE)

write.csv(metrics_scaled_byTime,
          file.path(stats_dir, "RMSE_MAE_byTime_Rescaled.csv"), row.names = FALSE)

summary_txt <- file.path(stats_dir, "Model_vs_RescaledData_Statistical_Summary.txt")
sink(summary_txt)

cat("========== MODEL vs RESCALED DATA :: STATISTICAL SUMMARY ==========\n\n")

cat(">>> Overall Metrics:\n")
print(metrics_scaled_overall)

cat("\n>>> RMSE / MAE / Pearson r by Cell Type:\n")
print(metrics_scaled_byCT)

cat("\n>>> RMSE / MAE by Timepoint:\n")
print(metrics_scaled_byTime)

sink()
cat("Statistical summary (rescaled data) written.\n")

#  == Function to generate plots using the calculated data  == 
plot_boxplot_comparison <- function(data, title, filename) {
  long_df <- data %>%
    pivot_longer(cols = c(RescaledDataPercent, ModelPercent),
                 names_to = "Source", values_to = "Percent") %>%
    mutate(Source = recode(Source, "RescaledDataPercent" = "Data", "ModelPercent" = "Model"))
  
  plot <- ggplot(long_df, aes(x = CellType, y = Percent, fill = Source)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7, position = position_dodge(0.75)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                alpha = 0.4, size = 1.2) +
    labs(title = title, x = "Cell Type", y = "Percentage of Total Cells", fill = "Source") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(stats_dir, filename), plot, width = 9, height = 6)
  print(plot)
}

plot_boxplot_comparison(combined_all,
                        "Boxplot: Model vs Rescaled Data (All Modeled Cell Types)",
                        "Boxplot_Model_vs_RescaledData_AllCells.png")

# == Other variation of the function == 
plot_facet_boxplot_by_day <- function(data, title, filename) {
  long_df <- data %>%
    pivot_longer(cols = c(RescaledDataPercent, ModelPercent),
                 names_to = "Source", values_to = "Percent") %>%
    mutate(Source = recode(Source, "RescaledDataPercent" = "Data", "ModelPercent" = "Model"))
  
  plot <- ggplot(long_df, aes(x = CellType, y = Percent, fill = Source)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65, position = position_dodge(0.7)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                alpha = 0.4, size = 1) +
    facet_wrap(~ Days, nrow = 2, scales = "free_y") +
    labs(title = title,
         x = "Cell Type", y = "Percentage of Total Cells", fill = "Source") +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(file.path(stats_dir, filename), plot, width = 12, height = 7)
  print(plot)
}

plot_facet_boxplot_by_day(combined_all,
                          "Faceted Boxplot: Model vs Rescaled Data (By Day)",
                          "Facet_Boxplot_Rescaled_vs_Model_byDay.png")

# ====================================================================== Final Plot: Rescaled True Proportions (All 7 Cell Types) ======================================================================

# == Set output directory
final_plot_dir <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Plots"

# Define custom color palette (feel free to change if you like)
final_colors <- c(
  "ASCs" = "#e53a46",
  "QSCs" = "#7fb918",
  "Neutrophils" = "#447a9c",
  "Myocytes" = "#f3a76c",
  "Monocytes" = "#f8844a",
  "M1 Macrophages" = "#529291",
  "M2 Macrophages" = "#7fc9a5"
)

# == Plot using RescaledDataPercent
final_plot <- ggplot(combined_all, aes(x = factor(Days), y = RescaledDataPercent, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ CellType, scales = "free_y") +
  labs(
    title = "Rescaled Proportions of Modeled Cell Types Across Timepoints",
    x = "Days Post-Injury", y = "Rescaled Percentage of Cells"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = final_colors)

# == Save plot
ggsave(
  filename = file.path(final_plot_dir, "Barplot_Rescaled_True_Proportions.png"),
  plot = final_plot, width = 10, height = 6, dpi = 300
)

print(final_plot)
cat("\nPlot saved to:", file.path(final_plot_dir, "Barplot_Rescaled_True_Proportions.png"), "\n")

# ====================================================================== Plot: Model vs Rescaled Data (All Cell Types) ======================================================================

# == Output path
plots_dir <- "D:/Masters/Processing Codes/Southerland_GSE227075_Processing/True Proportions/Plots"

# Shared aesthetics
shared_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

color_scheme <- c("Data" = "orange", "Model" = "red") 

# == Plot
plot_all <- ggplot(combined_all, aes(x = Days, group = CellType)) +
  geom_line(aes(y = RescaledDataPercent, color = "Data"), linewidth = 1) +
  geom_line(aes(y = ModelPercent, color = "Model", linetype = "Model"), linewidth = 1) +
  geom_point(aes(y = RescaledDataPercent, color = "Data", shape = "Data"), size = 2) +
  geom_point(aes(y = ModelPercent, color = "Model", shape = "Model"), size = 2) +
  facet_wrap(~ CellType, scales = "free_y") +
  labs(
    title = "Model vs Rescaled Data: All Modeled Cell Types",
    x = "Days Post-Injury", y = "Percentage of Total Modeled Cells",
    color = "Source", linetype = "Source", shape = "Source"
  ) +
  scale_color_manual(values = color_scheme) +
  shared_theme

# == Save plot
ggsave(file.path(plots_dir, "Model_vs_RescaledData_AllModeledCells.png"),
       plot_all, width = 10, height = 6, dpi = 300)

print(plot_all)
cat("\nSaved plot to:", file.path(plots_dir, "Model_vs_RescaledData_AllModeledCells.png"), "\n")

