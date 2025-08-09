# MuscleRegeneration_ODE_scRNAseq
R scripts and data processing pipelines for refining and testing an scRNA-seq-informed ODE model of skeletal muscle regeneration under toxin-induced (McKellar and Micheli) and ischemic (Southerland) injury conditions.

The scripts below were developed and tested in R (version 4.4.2) and require the Seurat, Harmony, and deSolve packages, among others.

Repository Structure and Script Descriptions
McKellar DeMicheli

McKellar_Micheli_Extraction.R: Extraction of relevant datasets (data, metadata, and proportions) from McKellar and Micheli datasets.

McKellar_New_Annotation.R: Initial annotation script for McKellar dataset.

Micheli_New_Annotation.R: Initial annotation script for Micheli dataset.

McKellar_FinalAnnotation_Subclustering.R: Final subclustering and annotation of McKellar dataset and its MuSC populations (Figure 1).

Micheli_FinalAnnotation_Subclustering.R: Final subclustering and annotation of Micheli dataset and its MuSC populations (Figure 2).

McKellar_Micheli_TrueProportions_Pipeline.R: Full pipeline for calculating relevant cell type proportions for ODE model input.

McKellar_Micheli_plot_replicate_proportions.R: Generates replicate-level cell type proportion plots for McKellar and Micheli datasets (Figure 3).

McKellar_Micheli_Regen_Model_Own_Data.R: Runs the ODE regeneration model using processed McKellar/Micheli data, along with sensitivity analysis (Figures 4 and 5).

Southerland
Southerland_GSE227075_Initial_Processing_and_Overall_Annotation.R: Initial QC, processing, and global annotation for Southerland dataset (Figure 6A).

Southerland_GSE227075_Myeloid_Subclustering_and_Annotation.R: Subclustering and annotation of myeloid populations in Southerland dataset (Figure 6B).

ode_fit_export_percents.R: Exports percentage-based ODE model fits for comparison with empirical data.

Southerland_GSE227075_All_Proportions.R: Calculates cell type proportions for all populations in Southerland dataset.

Compare_Model_vs_Southerland_Data.R: Calculates the relevant Southerland cell type proportions and compares ODE model outputs to Southerland ischemic injury dataset proportions (Figure 7).


