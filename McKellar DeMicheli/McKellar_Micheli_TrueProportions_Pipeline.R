# ===========================================================================
# Purpose: Calculate "true" ODE model proportions for relevant cell types
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
  library(tidyr)
})

# =======================      Directories      =======================

#base_dir <- "D:/Masters/Processing Codes/Mckellar_Cosgrove_Annotation"
base_dir <- "D:/Masters/Processing Codes/Micheli_Cosgrove_Annotation"

annot_dir <- file.path(base_dir, "Annotation", "Final UMAP Annotation")
subcluster_dir <- file.path(annot_dir, "MuSC_Subclustering")

trueprop_dir <- file.path(base_dir, "True Proportions")
log_dir <- file.path(trueprop_dir, "Logs")
prop_dir <- file.path(trueprop_dir, "Proportions")
plots_dir <- file.path(trueprop_dir, "Plots")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prop_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# =======================     Load Objects      =======================

cat(blue$bold("\n Loading final annotated main object...\n"))
#main_obj <- readRDS(file.path(annot_dir, "mckellar_annotated_final.rds"))
main_obj <- readRDS(file.path(annot_dir, "micheli_annotated_final.rds"))

cat(blue$bold("\n Loading MuSC subclustered object...\n"))
#musc_obj <- readRDS(file.path(subcluster_dir, "mckellar_MuSCs_annotated_final.rds"))
musc_obj <- readRDS(file.path(subcluster_dir, "micheli_MuSCs_annotated_final.rds"))

# =======================      Raw counts + props      =======================

compute_raw_counts_and_props <- function(main_obj, prop_dir, plots_dir, log_dir) {
  # ================================
  # Create output folders
  # ================================
  starting_prop_dir <- file.path(prop_dir, "Starting")
  starting_plots_dir <- file.path(plots_dir, "Starting")
  starting_log_file <- file.path(log_dir, "Starting", "Starting_Log.txt")
  
  dir.create(starting_prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(starting_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(starting_log_file), recursive = TRUE, showWarnings = FALSE)
  
  sink(starting_log_file, split = TRUE)  # Log everything printed
  
  cat(crayon::blue$bold("\n Starting compute_raw_counts_and_props()...\n"))
  
  # ================================
  # Extract metadata & raw counts
  # ================================
  meta_main <- main_obj@meta.data
  
  raw_counts <- meta_main %>%
    group_by(sample, injury.days, final_annotation) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  cat(crayon::green$bold("\nRaw counts computed (head):\n"))
  print(head(raw_counts))
  
  # ================================
  # Raw proportions
  # ================================
  raw_props <- raw_counts %>%
    group_by(sample, injury.days) %>%
    mutate(proportion = cell_count / sum(cell_count)) %>%
    ungroup()
  
  cat(crayon::green$bold("\n Raw proportions computed (head):\n"))
  print(head(raw_props))
  
  # ================================
  # Sanity check: sums ~1 per sample/day
  # ================================
  sums_check <- raw_props %>%
    group_by(sample, injury.days) %>%
    summarise(total = sum(proportion), .groups = "drop")
  
  cat(crayon::blue$bold("\n Proportion sums per sample/day (should be ~1):\n"))
  print(sums_check)
  
  # ================================
  # Mean by day
  # ================================
  props_by_day <- raw_props %>%
    group_by(injury.days, final_annotation) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop")
  
  cat(crayon::blue$bold("\n Proportions mean by day:\n"))
  print(props_by_day)
  
  # ================================
  # Mean by replicate
  # ================================
  props_by_replicate <- raw_props %>%
    group_by(sample, final_annotation) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop")
  
  cat(crayon::blue$bold("\n Proportions mean by replicate:\n"))
  print(props_by_replicate)
  
  # ================================
  # Save CSVs
  # ================================
  write.csv(raw_props, file.path(starting_prop_dir, "Raw_Proportions_AllAnnotations.csv"), row.names = FALSE)
  write.csv(props_by_day, file.path(starting_prop_dir, "Raw_Proportions_ByDay.csv"), row.names = FALSE)
  write.csv(props_by_replicate, file.path(starting_prop_dir, "Raw_Proportions_ByReplicate.csv"), row.names = FALSE)
  write.csv(sums_check, file.path(starting_prop_dir, "Proportion_Sums_Check.csv"), row.names = FALSE)
  
  cat(crayon::green$bold("\n All CSVs saved to: "), starting_prop_dir, "\n")
  
  # ================================
  # Plots — replicate-by-day
  # ================================

  p1 <- ggplot(raw_props, aes(x = final_annotation, y = proportion, fill = final_annotation)) +
    geom_bar(stat = "identity") +
    facet_grid(sample ~ injury.days) +
    scale_x_discrete(drop = TRUE) +   # Drop unused levels!
    theme_minimal() +
    ggtitle("Raw Proportions — Replicate by Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(starting_plots_dir, "Raw_Proportions_ReplicateByDay.png"), p1, width = 14, height = 8, dpi = 300)
  
  cat(crayon::green$bold("\n Saved replicate-by-day plot (with empty bars dropped)\n"))
  
  # ================================
  # Plots — mean by day
  # ================================
  p2 <- ggplot(props_by_day, aes(x = final_annotation, y = mean_proportion, fill = final_annotation)) +
    geom_bar(stat = "identity") +
    facet_wrap(~injury.days, scales = "free_y") +
    theme_minimal() +
    ggtitle("Raw Mean Proportions — By Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(starting_plots_dir, "Raw_Proportions_ByDay.png"), p2, width = 12, height = 8, dpi = 300)
  
  cat(crayon::green$bold("\n Saved mean by day plot\n"))
  
  # ================================
  # Plots — mean by replicate
  # ================================
  p3 <- ggplot(props_by_replicate, aes(x = final_annotation, y = mean_proportion, fill = final_annotation)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample, scales = "free_y") +
    theme_minimal() +
    ggtitle("Raw Mean Proportions — By Replicate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(starting_plots_dir, "Raw_Proportions_ByReplicate.png"), p3, width = 12, height = 8, dpi = 300)
  
  cat(crayon::green$bold("\n Saved mean by replicate plot\n"))
  
  # ================================
  # Final
  # ================================
  cat(crayon::magenta$bold("\n Done: Raw counts & proportions complete!\n"))
  
  sink()  # Stop logging
  
  # Return all outputs
  return(list(
    raw_props = raw_props,
    props_by_day = props_by_day,
    props_by_replicate = props_by_replicate,
    sums_check = sums_check
  ))
}

# =======================      Run the function 1     =======================

result <- compute_raw_counts_and_props(
  main_obj = main_obj,
  prop_dir = prop_dir,
  plots_dir = plots_dir,
  log_dir = log_dir
)

# ================      Mean by day and subsequent isolation...      ================

compute_day_and_replicate_lumped_props <- function(raw_props, prop_dir, plots_dir, log_dir) {
  # ================================
  # Create output folders
  # ================================
  day_prop_dir <- file.path(prop_dir, "DayLevel")
  day_plots_dir <- file.path(plots_dir, "DayLevel")
  day_log_file <- file.path(log_dir, "DayLevel", "DayLevel_Log.txt")
  
  dir.create(day_prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(day_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(day_log_file), recursive = TRUE, showWarnings = FALSE)
  
  sink(day_log_file, split = TRUE)  # Log everything printed
  
  cat(crayon::blue$bold("\n Starting compute_day_and_replicate_lumped_props()...\n"))
  
  cat(crayon::yellow$bold("\n️ Equation: For each (injury.day, final_annotation), mean proportion = average of replicate-level props.\n"))
  
  # ================================
  # Raw mean by day
  # ================================
  true_props_day <- raw_props %>%
    group_by(injury.days, final_annotation) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop")
  
  cat(crayon::green$bold("\n Mean proportions by day (head):\n"))
  print(head(true_props_day))
  
  write.csv(true_props_day, file.path(day_prop_dir, "Raw_Proportions_MeanByDay.csv"), row.names = FALSE)
  
  # ================================
  # Raw mean by replicate (double-check)
  # ================================
  true_props_replicate <- raw_props %>%
    group_by(sample, injury.days, final_annotation) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop")
  
  cat(crayon::green$bold("\n Mean proportions by replicate (head):\n"))
  print(head(true_props_replicate))
  
  write.csv(true_props_replicate, file.path(day_prop_dir, "Raw_Proportions_MeanByReplicate.csv"), row.names = FALSE)
  
  # =======================================
  #  Equation: Map fine clusters -> lumps
  # =======================================
  cat(crayon::yellow$bold("\n️ Mapping clusters to lumps:\n"))
  cat("M1a + M1b → M1\n")
  cat("M0/M + Ly6c → Monocytes\n")
  cat("Others unchanged (M2, Neutrophils, MuSCs)\n\n")
  
  # -------------------------------
  # Day-level lumps
  # -------------------------------
  lumps_day <- true_props_day %>%
    mutate(
      true_type = case_when(
        final_annotation %in% c("M1a Macrophages", "M1b Macrophages") ~ "M1",
        final_annotation %in% c("M0/M Monocytes", "Ly6c Monocytes") ~ "Monocytes",
        final_annotation == "M2 Macrophages" ~ "M2",
        final_annotation == "Neutrophils" ~ "Neutrophils",
        final_annotation == "MuSCs" ~ "MuSCs",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(true_type)) %>%
    group_by(injury.days, true_type) %>%
    summarise(mean_proportion = sum(mean_proportion), .groups = "drop")
  
  # Normalize so lumps sum to 1 per day
  lumps_day_norm <- lumps_day %>%
    group_by(injury.days) %>%
    mutate(mean_proportion = mean_proportion / sum(mean_proportion)) %>%
    ungroup()
  
  cat(crayon::yellow$bold("\n Lumped means by day (normalized):\n"))
  print(lumps_day_norm)
  
  sums_day <- lumps_day_norm %>%
    group_by(injury.days) %>%
    summarise(total = sum(mean_proportion))
  
  cat(crayon::blue$bold("\n Sums of normalized day-level lumps (should be 1):\n"))
  print(sums_day)
  
  # Fill missing lumps with 0
  days <- unique(lumps_day_norm$injury.days)
  lumps <- c("M1", "M2", "Monocytes", "Neutrophils", "MuSCs")
  
  lumps_day_filled <- lumps_day_norm %>%
    tidyr::complete(injury.days = days, true_type = lumps, fill = list(mean_proportion = 0))
  
  write.csv(lumps_day_filled, file.path(day_prop_dir, "Immune_MuSCs_Combined_MeansByDay_Normalized.csv"), row.names = FALSE)
  
  # -------------------------------
  # Replicate-level lumps
  # -------------------------------
  lumps_replicate <- raw_props %>%
    mutate(
      true_type = case_when(
        final_annotation %in% c("M1a Macrophages", "M1b Macrophages") ~ "M1",
        final_annotation %in% c("M0/M Monocytes", "Ly6c Monocytes") ~ "Monocytes",
        final_annotation == "M2 Macrophages" ~ "M2",
        final_annotation == "Neutrophils" ~ "Neutrophils",
        final_annotation == "MuSCs" ~ "MuSCs",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(true_type)) %>%
    group_by(sample, injury.days, true_type) %>%
    summarise(final_count = sum(proportion), .groups = "drop") %>%
    group_by(sample, injury.days) %>%
    mutate(final_proportion = final_count / sum(final_count)) %>%
    ungroup()
  
  sums_replicate <- lumps_replicate %>%
    group_by(sample, injury.days) %>%
    summarise(total = sum(final_proportion))
  
  cat(crayon::yellow$bold("\n Lumped replicate-level props (normalized, head):\n"))
  print(head(lumps_replicate))
  
  cat(crayon::blue$bold("\n Sums of replicate-level lumps (should be 1 per replicate/day):\n"))
  print(sums_replicate)
  
  write.csv(lumps_replicate, file.path(day_prop_dir, "Immune_MuSCs_Combined_MeansByReplicate_Normalized.csv"), row.names = FALSE)
  
  # -------------------------------
  # Plots - lumps by day
  # -------------------------------
  p_day <- ggplot(lumps_day_filled, aes(x = true_type, y = mean_proportion, fill = true_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~injury.days, scales = "free_y") +
    theme_minimal() +
    ggtitle("Lumped Normalized Means — By Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(day_plots_dir, "Immune_MuSCs_Lumped_MeanByDay.png"), p_day, width = 12, height = 8, dpi = 300)
  
  cat(crayon::green$bold("\n Saved lumped day-level plot\n"))
  
  # -------------------------------
  # Plots - lumps by replicate
  # -------------------------------
  p_replicate <- ggplot(lumps_replicate, aes(x = true_type, y = final_proportion, fill = true_type)) +
    geom_bar(stat = "identity") +
    facet_grid(sample ~ injury.days) +
    scale_x_discrete(drop = TRUE) +
    theme_minimal() +
    ggtitle("Lumped Normalized Proportions — Replicate by Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(day_plots_dir, "Immune_MuSCs_Lumped_ByReplicate.png"), p_replicate, width = 14, height = 8, dpi = 300)
  
  cat(crayon::green$bold("\n Saved lumped replicate-level plot\n"))
  
  # -------------------------------
  # Final wrap
  # -------------------------------
  cat(crayon::magenta$bold("\n Done: Day-level & replicate-level lumps computed + saved!\n"))
  
  sink()
  
  return(list(
    lumps_day_filled = lumps_day_filled,
    lumps_replicate = lumps_replicate
  ))
}

# =======================      Run the function 2     =======================

results_day <- compute_day_and_replicate_lumped_props(
  raw_props = result$raw_props,   # pull it from the previous function’s output
  prop_dir = prop_dir,
  plots_dir = plots_dir,
  log_dir = log_dir
)

# ============================= Extract subcluster fractions + Final Day True proportions =============================

expand_musc_lump_daylevel <- function(lumps_day_filled, musc_obj, prop_dir, plots_dir, log_dir) {
  # ================================
  # Setup output folders
  # ================================
  exp_prop_dir <- file.path(prop_dir, "Expansion")
  exp_plots_dir <- file.path(plots_dir, "Expansion")
  exp_log_file <- file.path(log_dir, "Expansion", "Expansion_Log.txt")
  
  dir.create(exp_prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exp_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(exp_log_file), recursive = TRUE, showWarnings = FALSE)
  
  sink(exp_log_file, split = TRUE)
  
  cat(crayon::blue$bold("\n Starting expand_musc_lump_daylevel()...\n"))
  
  # ================================
  # Step 1: Extract MuSC lumps
  # ================================
  musc_means <- lumps_day_filled %>%
    filter(true_type == "MuSCs") %>%
    rename(muSCs_mean = mean_proportion)
  
  cat(crayon::yellow$bold("\n️ Extracted MuSC lump means by day:\n"))
  print(musc_means)
  
  write.csv(musc_means, file.path(exp_prop_dir, "MuSC_Lump_Means_ByDay.csv"), row.names = FALSE)
  
  cat("\n Equation: For each day, MuSC lump = mean proportion BEFORE splitting.\n")
  
  # ================================
  # Step 2: Compute subcluster fractions
  # ================================
  meta_musc <- musc_obj@meta.data
  
  musc_total <- meta_musc %>%
    group_by(sample, injury.days) %>%
    summarise(total_count = n(), .groups = "drop")
  
  musc_subcounts <- meta_musc %>%
    group_by(sample, injury.days, muscs_annotation) %>%
    summarise(sub_count = n(), .groups = "drop")
  
  musc_subprops <- musc_subcounts %>%
    left_join(musc_total, by = c("sample", "injury.days")) %>%
    mutate(sub_fraction = sub_count / total_count)
  
  cat(crayon::yellow$bold("\nℹ️ Subcluster fractions per (sample, day):\n"))
  print(head(musc_subprops))
  
  # fractions sum to ~1 for each replicate/day
  sub_sum_check <- musc_subprops %>%
    group_by(sample, injury.days) %>%
    summarise(total_fraction = sum(sub_fraction))
  
  cat(crayon::blue$bold("\n Sanity check: fractions sum to ~1 per replicate/day:\n"))
  print(sub_sum_check)
  
  # Compute mean subcluster fraction by day
  musc_subprops_day <- musc_subprops %>%
    group_by(injury.days, muscs_annotation) %>%
    summarise(mean_sub_fraction = mean(sub_fraction), .groups = "drop")
  
  # Fill in missing subclusters with explicit zero
  musc_expected <- c("QSCs", "ASCs", "Differentiating Myocytes")
  days <- unique(musc_subprops_day$injury.days)
  
  musc_subprops_day <- musc_subprops_day %>%
    tidyr::complete(
      injury.days = days,
      muscs_annotation = musc_expected,
      fill = list(mean_sub_fraction = 0)
    )
  
  cat(crayon::yellow$bold("\n MuSC subcluster fractions (mean by day, filled):\n"))
  print(musc_subprops_day)
  
  write.csv(musc_subprops_day, file.path(exp_prop_dir, "MuSC_Subcluster_Fractions_ByDay.csv"), row.names = FALSE)
  
  cat("\n Equation: mean_sub_fraction = avg(sub_count / total MuSC) across replicates — missing combinations filled with 0.\n")
  
  # ================================
  # Expand lumps to subtypes
  # ================================
  musc_expanded <- musc_subprops_day %>%
    left_join(musc_means, by = "injury.days") %>%
    mutate(
      final_mean = muSCs_mean * mean_sub_fraction,
      true_type = muscs_annotation
    ) %>%
    select(injury.days, true_type, final_mean)
  
  cat(crayon::yellow$bold("\n Expanded MuSC subtypes:\n"))
  print(musc_expanded)
  
  write.csv(musc_expanded, file.path(exp_prop_dir, "MuSC_Expanded_ByDay.csv"), row.names = FALSE)
  
  cat("\n Equation: final_mean_subtype = (MuSC lump) × (mean_sub_fraction).\n")

  # ================================
  # Combine with immune lumps
  # ================================
  
  immune_only <- lumps_day_filled %>%
    filter(true_type != "MuSCs") %>%
    rename(final_mean = mean_proportion)
  
  # Combine and order: day first, then biologically logical order
  true_props_final <- bind_rows(immune_only, musc_expanded) %>%
    arrange(
      injury.days,
      factor(true_type, levels = c(
        "M1", "M2", "Monocytes", "Neutrophils",
        "QSCs", "ASCs", "Differentiating Myocytes"
      ))
    )
  
  cat(crayon::magenta$bold("\n Combined final true means (immune + expanded MuSC subtypes, ordered):\n"))
  print(true_props_final)
  
  # Normalize by day
  true_props_final_norm <- true_props_final %>%
    group_by(injury.days) %>%
    mutate(final_mean_norm = final_mean / sum(final_mean)) %>%
    ungroup()
  
  cat(crayon::blue$bold("\n FINAL sums by day (should be 1):\n"))
  print(true_props_final_norm %>%
          group_by(injury.days) %>%
          summarise(total = sum(final_mean_norm)))
  
  write.csv(true_props_final_norm, file.path(exp_prop_dir, "True_Proportions_Final_ByDay_Normalized.csv"), row.names = FALSE)
  
  cat(crayon::green$bold("\n Saved final normalized ODE-ready means by day.\n"))
  
  # ================================
  # Plots - Subcluster fractions
  # ================================
  p_subs <- ggplot(musc_subprops_day, aes(x = muscs_annotation, y = mean_sub_fraction, fill = muscs_annotation)) +
    geom_bar(stat = "identity") +
    facet_wrap(~injury.days) +
    theme_minimal() +
    ggtitle("MuSC Subcluster Fractions by Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(exp_plots_dir, "MuSC_Subcluster_Fractions_ByDay.png"), p_subs, width = 12, height = 8, dpi = 300)
  cat(crayon::green$bold("\n Saved subcluster fractions plot.\n"))
  
  # ================================
  # Plots - Expanded MuSC subtypes
  # ================================
  p_expanded <- ggplot(musc_expanded, aes(x = true_type, y = final_mean, fill = true_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~injury.days) +
    theme_minimal() +
    ggtitle("Expanded MuSC Subtypes by Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(exp_plots_dir, "MuSC_Expanded_ByDay.png"), p_expanded, width = 12, height = 8, dpi = 300)
  cat(crayon::green$bold("\n Saved expanded MuSC subtypes plot.\n"))
  
  # ================================
  # Plots - Final ODE-ready combined
  # ================================
  p_final <- ggplot(true_props_final_norm, aes(x = true_type, y = final_mean_norm, fill = true_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~injury.days) +
    theme_minimal() +
    ggtitle("Final ODE-Ready True Proportions by Day") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(exp_plots_dir, "True_Proportions_Final_ByDay_Normalized.png"), p_final, width = 12, height = 8, dpi = 300)
  cat(crayon::green$bold("\n Saved final combined ODE-ready plot.\n"))
  
  # ================================
  # Wrap up
  # ================================
  cat(crayon::magenta$bold("\n Done: MuSC lumps expanded & final day-level ODE-ready proportions saved!\n"))
  
  sink()
  
  return(list(
    musc_means = musc_means,
    musc_subprops_day = musc_subprops_day,
    musc_expanded = musc_expanded,
    true_props_final_norm = true_props_final_norm
  ))
}

# =======================      Run the function 3     =======================

result_expansion <- expand_musc_lump_daylevel(
  lumps_day_filled = results_day$lumps_day_filled,
  musc_obj = musc_obj,
  prop_dir = prop_dir,
  plots_dir = plots_dir,
  log_dir = log_dir
)

# =======================      Run the function 4     =======================

build_ode_input_df <- function(true_props_final_norm) {
  cat(crayon::blue$bold("\n Building ODE input dataframe...\n"))
  
  # ------------------------
  # Pick only relevant columns
  # ------------------------
  df <- true_props_final_norm %>%
    select(injury.days, true_type, final_mean_norm)
  
  # ------------------------
  # Pivot to wide: rows = days, cols = types
  # ------------------------
  df_wide <- df %>%
    tidyr::pivot_wider(
      names_from = true_type,
      values_from = final_mean_norm,
      values_fill = list(final_mean_norm = 0)  # fill missing with 0 first
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(. == 0, NA_real_, .)))
  
  # ------------------------
  # Add placeholder columns for dead myonuclei & dead neutrophils
  # ------------------------
  df_wide <- df_wide %>%
    mutate(
      Md = NA_real_,
      Nd = NA_real_
    )
  
  # ------------------------------------------------
  # Reorder columns to match Renad's ODE structure
  # order: time, Md, N, Nd, M, M1, M2, QSC, ASC, Mc
  # ------------------------------------------------
  df_wide <- df_wide %>%
    rename(
      time = injury.days,
      N = Neutrophils,
      M = Monocytes,
      QSC = QSCs,
      ASC = ASCs,
      Mc = `Differentiating Myocytes`
    ) %>%
    select(time, Md, N, Nd, M, M1, M2, QSC, ASC, Mc)
  
  cat(crayon::green$bold("\n ODE input dataframe built:\n"))
  print(df_wide)
  
  return(df_wide)
}

# usage
#McKellarPercent <- build_ode_input_df(result_expansion$true_props_final_norm)

DeMicheliPercent <- build_ode_input_df(result_expansion$true_props_final_norm)
                                      
# ======================= Replicate Function 5 ============================

expand_musc_lump_replicatelevel <- function(lumps_replicate, musc_obj, prop_dir, plots_dir, log_dir) {
  # ================================
  # Setup output folders
  # ================================
  exp_prop_dir <- file.path(prop_dir, "Expansion_Replicates")
  exp_plots_dir <- file.path(plots_dir, "Expansion_Replicates")
  exp_log_file <- file.path(log_dir, "Expansion_Replicates", "Replicates_Log.txt")
  
  dir.create(exp_prop_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exp_plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(exp_log_file), recursive = TRUE, showWarnings = FALSE)
  
  sink(exp_log_file, split = TRUE)
  
  cat(crayon::blue$bold("\n Starting expand_musc_lump_replicatelevel()...\n"))
  
  # ================================
  # Step 1: Extract MuSC lumps per replicate × day
  # ================================
  musc_means <- lumps_replicate %>%
    filter(true_type == "MuSCs") %>%
    rename(muSCs_count = final_proportion)
  
  cat(crayon::yellow$bold("\n Extracted MuSC lumps by replicate × day:\n"))
  print(head(musc_means))
  
  write.csv(musc_means, file.path(exp_prop_dir, "MuSC_Lump_Means_ByReplicate.csv"), row.names = FALSE)
  
  # =========================================================
  # Step 2: Compute subcluster fractions per replicate × day
  # =========================================================
  meta_musc <- musc_obj@meta.data
  
  musc_total <- meta_musc %>%
    group_by(sample, injury.days) %>%
    summarise(total_count = n(), .groups = "drop")
  
  musc_subcounts <- meta_musc %>%
    group_by(sample, injury.days, muscs_annotation) %>%
    summarise(sub_count = n(), .groups = "drop")
  
  musc_subprops <- musc_subcounts %>%
    left_join(musc_total, by = c("sample", "injury.days")) %>%
    mutate(sub_fraction = sub_count / total_count)
  
  cat(crayon::yellow$bold("\nℹ️ Subcluster fractions per replicate × day:\n"))
  print(head(musc_subprops))
  
  sub_sum_check <- musc_subprops %>%
    group_by(sample, injury.days) %>%
    summarise(total_fraction = sum(sub_fraction))
  
  cat(crayon::blue$bold("\n Sanity check: fractions sum to ~1 per replicate × day:\n"))
  print(sub_sum_check)
  
  musc_expected <- c("QSCs", "ASCs", "Differentiating Myocytes")
  samples <- unique(musc_subprops$sample)
  days <- unique(musc_subprops$injury.days)
  
  musc_subprops_filled <- musc_subprops %>%
    tidyr::complete(
      sample = samples,
      injury.days = days,
      muscs_annotation = musc_expected,
      fill = list(sub_fraction = 0)
    )
  
  cat(crayon::yellow$bold("\n Subcluster fractions filled for missing subtypes:\n"))
  print(head(musc_subprops_filled))
  
  write.csv(musc_subprops_filled, file.path(exp_prop_dir, "MuSC_Subcluster_Fractions_ByReplicate.csv"), row.names = FALSE)
  
  # ==========================================================
  # Step 3: Expand lumps to MuSC subtypes per replicate × day
  # ==========================================================
  musc_expanded <- musc_subprops_filled %>%
    left_join(musc_means, by = c("sample", "injury.days")) %>%
    mutate(
      final_count = muSCs_count * sub_fraction,
      true_type = muscs_annotation
    ) %>%
    select(sample, injury.days, true_type, final_count)
  
  cat(crayon::yellow$bold("\n Expanded MuSC subtypes per replicate × day:\n"))
  print(head(musc_expanded))
  
  write.csv(musc_expanded, file.path(exp_prop_dir, "MuSC_Expanded_ByReplicate.csv"), row.names = FALSE)
  
  # ==============================================================
  # Step 4: Combine with immune lumps, normalize per replicate
  # ==============================================================
  immune_only <- lumps_replicate %>%
    filter(true_type != "MuSCs") %>%
    mutate(final_count = final_proportion)
  
  true_props_final <- bind_rows(immune_only, musc_expanded) %>%
    arrange(
      sample, injury.days,
      factor(true_type, levels = c(
        "M1", "M2", "Monocytes", "Neutrophils",
        "QSCs", "ASCs", "Differentiating Myocytes"
      ))
    ) %>%
    group_by(sample, injury.days) %>%
    mutate(final_proportion = final_count / sum(final_count)) %>%
    ungroup()
  
  cat(crayon::magenta$bold("\n Final true replicate-level proportions (normalized):\n"))
  print(head(true_props_final))
  
  sums_check <- true_props_final %>%
    group_by(sample, injury.days) %>%
    summarise(total = sum(final_proportion))
  
  cat(crayon::blue$bold("\n Sums per replicate × day (should be 1):\n"))
  print(sums_check)
  
  write.csv(true_props_final, file.path(exp_prop_dir, "True_Proportions_Final_ByReplicate_Normalized.csv"), row.names = FALSE)
  
  # ========================================================================
  # Step 5: Write ONE CSV + ONE plot per replicate (across ALL timepoints!)
  # ========================================================================
  replicate_list <- unique(true_props_final$sample)
  
  for (rep in replicate_list) {
    rep_dir <- file.path(exp_prop_dir, rep)
    rep_plot_dir <- file.path(exp_plots_dir, rep)
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(rep_plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    rep_df <- true_props_final %>%
      filter(sample == rep)
    
    write.csv(rep_df, file.path(rep_dir, paste0(rep, "_True_Proportions.csv")), row.names = FALSE)
    
    p <- ggplot(rep_df, aes(x = factor(injury.days), y = final_proportion, fill = true_type)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      ggtitle(paste0("True Proportions across Timepoints — ", rep)) +
      labs(x = "Injury Day", y = "Proportion", fill = "Cell Type") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    
    ggsave(file.path(rep_plot_dir, paste0(rep, "_True_Proportions.png")), p, width = 12, height = 8, dpi = 300)
    
    cat(crayon::green$bold(paste0("\n Saved CSV + timecourse plot for replicate ", rep, "\n")))
  }
  
  # ================================
  # Wrap up
  # ================================
  cat(crayon::magenta$bold("\n  Done: Replicate-level expansion complete - ODE-ready files across timepoints saved!\n"))
  
  sink()
  
  return(list(
    musc_means = musc_means,
    musc_subprops_filled = musc_subprops_filled,
    musc_expanded = musc_expanded,
    true_props_final = true_props_final
  ))
}

result_expansion_replicates <- expand_musc_lump_replicatelevel(
  lumps_replicate = results_day$lumps_replicate,
  musc_obj = musc_obj,
  prop_dir = prop_dir,
  plots_dir = plots_dir,
  log_dir = log_dir
)





