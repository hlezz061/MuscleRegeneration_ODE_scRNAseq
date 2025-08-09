##############################################################################
# Purpose: Overlay Mckellar and Micheli replicate-level proportions into
# a single grid figure with shared legend and line types
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(crayon)
})

# =======================      Configurable Paths     =======================

# Set input paths
base_mckellar <- "D:/Masters/Processing Codes/Mckellar_Cosgrove_Annotation/V4"
base_micheli  <- "D:/Masters/Processing Codes/Micheli_Cosgrove_Annotation/V3"

input_mckellar <- file.path(base_mckellar, "True Proportions", "Proportions", "Expansion_Replicates", "True_Proportions_Final_ByReplicate_Normalized.csv")
input_micheli  <- file.path(base_micheli,  "True Proportions", "Proportions", "Expansion_Replicates", "True_Proportions_Final_ByReplicate_Normalized.csv")

output_path <- "D:/Masters/Processing Codes/Paper Figures"
output_file <- file.path(output_path, "New_Replicate_Timecourse_Overlay_GRID.png")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
}

# =======================      Load and Combine Data     =======================

cat(blue$bold("\n Loading and tagging datasets...\n"))

load_and_tag <- function(path, tag) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(dataset = tag,
           injury.days = as.numeric(injury.days))
}

df_mckellar <- load_and_tag(input_mckellar, "Mckellar")
df_micheli  <- load_and_tag(input_micheli, "Micheli")

combined_df <- bind_rows(df_mckellar, df_micheli)

cell_types_raw <- unique(combined_df$true_type)
cat(green$bold(paste0("Found ", length(cell_types_raw), " unique raw cell types\n")))

# =======================  Standardize labels + set desired order  =======================

desired_order <- c("Neutrophils", "Monocytes", "M1", "M2", "QSCs", "ASCs", "Differentiating Myocytes")

combined_df <- combined_df %>%
  mutate(
    # Map internal labels to display labels
    type_display = dplyr::recode(
      true_type,
      "N"  = "Neutrophils",
      "Neutrophils" = "Neutrophils",
      "M"  = "Monocytes",
      "Monocytes"   = "Monocytes",
      "M1" = "M1",
      "M2" = "M2",
      "QSC"  = "QSCs",
      "QSCs" = "QSCs",
      "ASC"  = "ASCs",
      "ASCs" = "ASCs",
      "Mc"   = "Differentiating Myocytes",
      "Myocytes" = "Differentiating Myocytes",
      "Differentiating Myocytes" = "Differentiating Myocytes",
      .default = true_type
    ),
    type_display = factor(type_display, levels = desired_order)
  )

# Only keep panels that exist in the data, but respect the desired order
cell_types <- desired_order[desired_order %in% levels(combined_df$type_display)]

# =======================      Generate Overlay Plots     =======================

plot_list <- list()

for (ct in cell_types) {
  cat(blue$bold(paste0(" Plotting: ", ct, "\n")))
  
  df_ct <- combined_df %>%
    filter(type_display == ct, !is.na(injury.days), !is.na(final_proportion))
  
  if (nrow(df_ct) == 0) {
    cat(yellow$bold(paste0(" Skipping ", ct, ": no valid data\n")))
    next
  }
  
  means_df <- df_ct %>%
    group_by(dataset, injury.days) %>%
    summarise(mean_prop = mean(final_proportion), .groups = "drop")
  
  p <- ggplot() +
    # Replicate dots
    geom_point(data = df_ct,
               aes(x = injury.days, y = final_proportion, color = dataset),
               size = 2.7, alpha = 0.9) +
    # Mean trend lines
    geom_line(data = means_df,
              aes(x = injury.days, y = mean_prop, linetype = dataset),
              color = "black", linewidth = 1.2) +
    geom_point(data = means_df,
               aes(x = injury.days, y = mean_prop),
               color = "black", shape = 18, size = 2.7) +
    scale_color_manual(values = c("Mckellar" = "blue", "Micheli" = "red")) +
    scale_linetype_manual(values = c("Mckellar" = "solid", "Micheli" = "dashed")) +
    theme_minimal() +
    labs(
      title = ct,
      x = "Days Post Injury",
      y = "Proportion",
      color = "Dataset",
      linetype = "Dataset"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.position = if (ct == cell_types[1]) "bottom" else "none",
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.major.x = element_blank()
    )
  
  plot_list[[ct]] <- p
}

# =======================      Assemble Final Grid      =======================

cat(green$bold("\n Combining plots into grid...\n"))

combined_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

combined_plot <- combined_plot +
  plot_annotation(
    title = "Time Course of Cell Type Proportions: Mckellar vs. Micheli",
    caption = "Dots = replicates | Line = mean | Solid = Mckellar | Dashed = Micheli",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 10))
    )
  )

# =======================      Save Output      =======================

print(combined_plot)

ggsave(output_file, combined_plot, width = 14, height = 10, dpi = 300)

# Save PDF version
ggsave(file.path(output_path, "Replicate_Timecourse_Overlay_GRID.pdf"), combined_plot, width = 14, height = 10)

cat(magenta$bold(paste0("\n Final overlay grid saved to: ", output_file, "\n")))















































