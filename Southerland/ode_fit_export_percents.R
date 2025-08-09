# ==================================================
# Code that runs the ODE model  
# And calculates the percentages that will be  
# Used in the comparison with the Southerland data
# ==================================================

# ====================================================================== Libraries ======================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(gplots)
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(RColorBrewer)
})

# ====================================================================== Save settings  ======================================================================

save_dir <- "D:/Masters/Processing Codes/Model Sim/Mouse_GSE227075_Sim 2"
dir.create(save_dir, showWarnings = FALSE)
dir.create(file.path(save_dir, "results"), showWarnings = FALSE)
dir.create(file.path(save_dir, "parameters"), showWarnings = FALSE)
dir.create(file.path(save_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(save_dir, "proportions"), showWarnings = FALSE)  


save_plot <- function(plot_obj, filename) {
  ggsave(filename = file.path(save_dir, "plots", filename),
         plot = plot_obj, width = 6, height = 4)
}

#  ====================================================================== Input data (own annotated data) ======================================================================

# == The data we are fitting. First, in percentages of all cells,
# == and then calibrated to cell counts per mm^3.

# == Start with the raw proportions (my own proportions)
McKellarPercent <- data.frame(
  time = c(0, 1, 2, 3.5, 5, 7),# Time in days post injury (DPI)
  Md = c(NA, NA, NA, NA, NA, NA), # Dead myonuclei not observed in data
  N = c(0.0778, 0.175, 0.0195, 0.00663, 0.00443, 0.00755), # Neutrophils
  Nd = c(NA, NA, NA, NA, NA, NA), # Dead neutrophils not observed in data, initial condition assumed zero
  M = c(0.194, 0.586, 0.185, 0.0516, 0.0348, 0.0878), # Monocytes
  M1 = c(0, 0.169, 0.719, 0.496, 0.0724, 0.0491), # M1 Macrophages
  M2 = c(0.251, 0.0666, 0.0751, 0.379, 0.695, 0.560), # M2 Macrophages
  QSC = c(0.477, 0.00150, 0.000917, 0.00414, 0.144, 0.228), # Quiescent Satellite Cells
  ASC = c(0, 0.00307, 0.000678, 0.0407, 0.0287, 0.0295), # Activated Satellite Cells
  Mc = c(0, 0, 0.000282, 0.0223, 0.0205, 0.0375) # Myocytes
)

# == Convert to percentage & round (except 'time')
McKellarPercent <- McKellarPercent %>%
  mutate(across(-time, ~ round(.x * 100, 2)))

# == McKellarPercent now holds the percentage values rounded to 2 decimal points
print(McKellarPercent)

# == Now for DeMicheli 

DeMicheliPercent <- data.frame(
  time = c(0, 2, 5, 7),
  Md = c(NA, NA, NA, NA),
  N = c(0.0487, 0.0425, 0.00149, 0.00323),
  Nd = c(NA, NA, NA, NA),
  M = c(0.0379, 0.262, 0.0430, 0.0220),
  M1 = c(0.0390, 0.543, 0.319, 0.108),
  M2 = c(0.404, 0.148, 0.497, 0.593),
  QSC = c(0.460, 0.00242, 0.0454, 0.151),
  ASC = c(0, 0, 0.0540, 0.0809),
  Mc = c(0.0103, 0.00108, 0.0404, 0.0418)
)

DeMicheliPercent <- DeMicheliPercent %>%
  mutate(across(-time, ~ round(.x * 100, 2)))

print(DeMicheliPercent)

# == This function calibrates a given cells percent dataframe, like either of the two above
# == to units of cells per cubic millimeter. It does so based on specified counts of
# == representative cell types taken from the literature that should be achieved at different
# == days post-injury.
CalibratePercents <- function(PercentDF) {
  CalibDF <- data.frame(PercentDF)
  NTimes <- nrow(PercentDF)
  NVars <- ncol(PercentDF)
  Times = PercentDF[,1]
  for (row in 1:NTimes) {
    # Figure out the scale factor for this row:
    RowTime <- PercentDF[row,1]
    RowScale <- 0
    if (RowTime==0) {
      RowScale <- 2500 / PercentDF[row,'QSC'] # QSC at Day 0 should be 2500
    }
    if (RowTime==1) {
      RowScale <- 4500 / PercentDF[row,'N'] # Neutrophils at Day 1 should be 4500
    }
    if (RowTime==2) {
      RowScale <- 8000 / PercentDF[row,'M1'] # M1 Macrophages at Day 2 should be 8000
    }
    if (RowTime==3.5) {
      RowScale <- 900 / PercentDF[row,'Mc'] # Myocytes at Day 3.5 should be 900
    }
    if (RowTime==5) {
      RowScale <- 7000 / PercentDF[row,'M2'] # M2 Macrophages at Day 5 should be 7000
    }
    if (RowTime==7) {
      RowScale <- 3000 / PercentDF[row,'M2'] # M2 Macrophages at Day 7 should be 3000
    }
    
    # == Then apply the scale to the row
    CalibDF[row,] = RowScale * CalibDF[row,]
    
  }
  
  # == Put the times -- because we scaled them above, duh
  CalibDF[,1] = Times
  return(CalibDF)
}


McKellarCellsPerMM3 <- CalibratePercents(McKellarPercent)
DeMicheliCellsPerMM3 <- CalibratePercents(DeMicheliPercent)

write.csv(McKellarCellsPerMM3, file.path(save_dir, "results", "McKellar_calibrated.csv"), row.names = FALSE)
write.csv(DeMicheliCellsPerMM3, file.path(save_dir, "results", "DeMicheli_calibrated.csv"), row.names = FALSE)

# ====================================================================== ODE Model ======================================================================

# == Initial states of cell type variables
InitialState <- c(Md = 30000, # Dead mynuclei
                  N = 0, # Neutrophils
                  Nd = 0, # Dead neutrophils
                  M = 0, # Monocytes
                  M1 = 0, # M1 Macrophages
                  M2 = 0, # M2 Macrophages
                  QSC = 2500, # Quiescent satellite cells
                  ASC = 0, # Activated satellite cells
                  Mc = 0 # Myocytes
)

# == Initial guesses at parameter values for differential equation model
# == (arrived at by a combination of reasoning of time-scales for
# == infiltration and exfiltration processes, as well as some amount of
# == "by hand" tweaking based on qualitative expectations and comparison
# == to McKellar data, but not DeMicheli data)
ParamsInit <- c(cNMd=0.0002, # Md removal by N
                cM1Md=0.00003, # Md removal by M1
                cNin=10, # N infiltration prompted by Md
                cNout=10, # N exfiltration
                cM1Nd=0.0001, # Nd removal by M1
                cMin=4, # M inflitration prompted by N
                cMM1=0.00005, # M differentiation to M1, prompted by Nd or Md
                cMout=1.2, # M exfiltration
                cM1out=1e-10, # M1 exflitration # Do we need this?
                cM1M2=1000, # M1 differentiation to M2
                cM2inhib=1000, # inhibition of M1 differentiation to M2
                cM2out=0.3, # M2 exfiltration
                cQSCN=0.0002, # QSC activation to ASC prompted by N
                cQSCMd=0.00002, # QSC activation to ASC prompted by Md
                cASCM2=0.0001, # ASC deactivation to QSC prompted by M2
                cASCpro=0.000005, # ASC proliferation prompted by M1
                cASCdiff=0.00005, # ASC differentiation to Mc prompted by M2
                cASCout=1e-10, # ASC exfiltration # Do we need this?
                cMcout=0.2 # Mc exit from system (really fusion into myotubes or something)
)

write.csv(data.frame(Parameter = names(ParamsInit), Value = ParamsInit),
          file.path(save_dir, "parameters", "Params_Initial.csv"), row.names = FALSE)

# == Observation times -- these are the McKellar days post injury, which
# == conveniently include the DeMicheli days post injury as a subset
ObservationTimes <- c(0, 1, 2, 3.5, 5, 7)


# == Function capturing regeneration dynamics  == 
MuscleRegenDerivs <- function(time, state, params) {
  with(as.list(c(state, params)), { 
    
    dMd <- -cNMd * N * Md -cM1Md * M1 * Md  
    dN <- cNin * Md - cNout * N - cNMd * N * Md 
    dNd <- cNMd * N * Md - cM1Nd * M1 * Nd
    dM <- cMin * N - cMM1 * M * (Nd + Md) - cMout * M
    dM1 <- cMM1 * M * (Nd + Md) - cM1out * M1 - (cM1M2 / (cM2inhib + Nd + Md)) * M1
    dM2 <- M1 * cM1M2 / (cM2inhib + Nd + Md) - cM2out * M2
    dQSC <- -cQSCN * QSC * N - cQSCMd * QSC * Md + cASCM2 * M2 * ASC
    dASC <- cQSCN * QSC * N + cQSCMd * QSC * Md - cASCM2 * ASC * M2 + cASCpro * M1 * ASC - cASCdiff * M2 * ASC - cASCout * ASC
    dMc <- cASCdiff * ASC * M2 - cMcout * Mc
    
    return(list(c(dMd, dN, dNd, dM, dM1, dM2, dQSC, dASC, dMc)))
  })
}


# ====================================================================== Function to perform one simulation trajectory for a given set of parameters ======================================================================

#  ==  Simulation Function  == 
MuscleSimulation <- function(Params,OutputTimes=c()) {
  
  # Figure out time points
  if (length(OutputTimes)==0) {
    OutputTimes = seq(from = 0, to = 7, by = 1/24)
  }
  
  # Call ODE simulator to solve dynamics
  Trajectory <- ode(y = InitialState, times = OutputTimes, func = MuscleRegenDerivs, parms = Params, method = 'lsoda')
  
  # Return solution
  return(Trajectory)
}

# == Function to plot the outcome of a simulation, optionally showing also McKellar and/or DeMicheli data
ShowSimulation <- function(TrajData,ShowMcK=FALSE,ShowDeM=FALSE) {
  
  # == Rearrange data into long lists with time and variable name as keys, and cells/mm3 values in Traj, McK or DeM columns
  TrajDF <- as.data.frame(TrajData)
  write.csv(TrajDF, file.path(save_dir, "results", "Trajectory_Init.csv"), row.names = FALSE)
  
  TrajLong <- TrajDF %>% pivot_longer(cols = -time, names_to = "variable", values_to = "Traj")
  McKLong <- McKellarCellsPerMM3 %>% pivot_longer(cols = -time, names_to = "variable", values_to = "McK")
  DeMLong <- DeMicheliCellsPerMM3 %>% pivot_longer(cols = -time, names_to = "variable", values_to = "DeM")
  
  # == Then join those, essentially filling in NA for McK and DeM wherever they don't have a value
  Joined_Temp <- left_join(TrajLong, McKLong, by = c("time", "variable")) 
  Joined <- left_join(Joined_Temp, DeMLong, by = c("time", "variable"))
  
  # == Finally, separate them again! Create plots for each variable
  plots <- lapply(unique(Joined$variable), function(var) {
    
    # == Extract trajectory and observed just for that variable
    JoinedVar = filter(Joined,variable==var)
    
    # == Plot it, depending on which observed data was requested
    if (!ShowMcK & !ShowDeM) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    } else if (ShowMcK & !ShowDeM) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + geom_point(aes(y=McK), color="blue") + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    } else if (!ShowMcK & ShowDeM) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + geom_point(aes(y=DeM), color="red") + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    } else {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + 
        geom_point(aes(y=McK), color="blue", alpha=0.5) + 
        geom_point(aes(y=DeM), color="red", alpha=0.5) + 
        labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw() + theme(text = element_text(size=7))
    }
    return(ggp)
  })
  
  # == Combine all the plots into one grid
  plot_grid <- do.call(grid.arrange, c(plots, ncol = 3, nrow = 3))
  ggsave(filename = file.path(save_dir, "plots", "Sim_Default_vs_Data.png"), plot = plot_grid, width = 10, height = 8)
  
  # == Print + save the combined plot grid
  print(plot_grid)
}


#  == How does it look with default parameters? == 
Traj <- MuscleSimulation(ParamsInit)
ShowSimulation(Traj,TRUE,FALSE)

# ====================================================================== Model fitting-related functions ======================================================================

# == Function to minimize - sum squared error of residuals of each
# == state variable, normalized by the maximum observed value of that state variable.
# == By default, the Params are evaluated with respect to McKellar dataset.
Error <- function(Params,CompareTo="McKellar") {
  
  # Simulate trajectory with the parameters Params
  if (CompareTo=="McKellar") {
    Traj <- MuscleSimulation(Params,OutputTimes=c(0,1,2,3.5,5,7))
    TrajDF <- as.data.frame(Traj)
    Empirical <- McKellarCellsPerMM3
  } else {
    if (CompareTo=="DeMicheli") {
      Traj <- MuscleSimulation(Params,OutputTimes=c(0,2,5,7))
      TrajDF <- as.data.frame(Traj)
      Empirical <- DeMicheliCellsPerMM3
    } else {
      print("Unknown dataset")
      return(+inf)
    }
  }
  
  # == Calculate error for each variable
  errorN <- sum(sqrt(((Empirical$N - TrajDF$N) / max(Empirical$N, na.rm = TRUE))^2))
  errorM <- sum(sqrt(((Empirical$M - TrajDF$M) / max(Empirical$M, na.rm = TRUE))^2))
  errorM1 <- sum(sqrt(((Empirical$M1 - TrajDF$M1) / max(Empirical$M1, na.rm = TRUE))^2))
  errorM2 <- sum(sqrt(((Empirical$M2 - TrajDF$M2) / max(Empirical$M2, na.rm = TRUE))^2))
  errorQSC <- sum(sqrt(((Empirical$QSC - TrajDF$QSC) / max(Empirical$QSC, na.rm = TRUE))^2))
  errorASC <- sum(sqrt(((Empirical$ASC - TrajDF$ASC) / max(Empirical$ASC, na.rm = TRUE))^2))
  errorMc <- sum(sqrt(((Empirical$Mc - TrajDF$Mc) / max(Empirical$Mc, na.rm = TRUE))^2))
  
  total_error <- errorN + errorM + errorM1 + errorM2 + errorQSC + errorASC + errorMc
  return(total_error)
}

# == Similar to the above function Error with parameter "McKellar" for dataset  == 
# == However, this function assumes that the logarithms of the parameters are input 
# == It takes the exponential of that input to compute the regular-scale parameters,
# == which it then feeds into Error. This is useful because it ensures regular-scale
# == parameters are always positive.
Error_LogParams <- function(LogParams,CompareTo="McKellar") {
  Params <- exp(LogParams)
  return(Error(Params,CompareTo))
}


# == Error of initial parameters == 
McK_ErrorInit <- Error(ParamsInit,"McKellar")

DeM_ErrorInit <- Error(ParamsInit,"DeMicheli")

writeLines(c(
  paste("Error on McKellar data with initial parameters:", McK_ErrorInit),
  paste("Error on DeMicheli data with initial parameters:", DeM_ErrorInit)
), file.path(save_dir, "results", "Error_Initial.txt"))

# == Parameter optimization by Nelder-Mead method  == 
y = optim(log(ParamsInit),Error_LogParams,method="Nelder-Mead")
ParamsOpt <- exp(y$par)
TrajOpt <- MuscleSimulation(ParamsOpt)
ShowSimulation(TrajOpt,TRUE,TRUE)

# == Error after optimization == 
opt_result <- optim(log(ParamsInit), Error_LogParams, method="Nelder-Mead")
ParamsOpt <- exp(opt_result$par)
write.csv(data.frame(Parameter = names(ParamsOpt), Value = ParamsOpt),
          file.path(save_dir, "parameters", "Params_Optimized.csv"), row.names = FALSE)


TrajOpt <- MuscleSimulation(ParamsOpt)
write.csv(as.data.frame(TrajOpt), file.path(save_dir, "results", "Trajectory_Optimized.csv"), row.names = FALSE)
ShowSimulation(TrajOpt,TRUE,TRUE) 

McK_ErrorOpt <- Error(ParamsOpt,"McKellar")
DeM_ErrorOpt <- Error(ParamsOpt,"DeMicheli")
writeLines(c(
  paste("Error on McKellar data with optimized parameters:", McK_ErrorOpt),
  paste("Error on DeMicheli data with optimized parameters:", DeM_ErrorOpt)
), file.path(save_dir, "results", "Error_Optimized.txt"))

# ====================================================================== Convert to Percentages ======================================================================

# == Output to be used for the Southerland data comparison == 

# == Load the optimized trajectory
Trajectory_Optimized <- read.csv(file.path(save_dir, "results", "Trajectory_Optimized.csv"))

# == Select relevant columns (excluding Md and Nd)
cell_columns <- c("N", "M", "M1", "M2", "QSC", "ASC", "Mc")
TrajCells <- Trajectory_Optimized[, cell_columns]

# == Compute row-wise total
RowTotal <- rowSums(TrajCells, na.rm = TRUE)

# == Convert to percentages (instead of raw proportions)
Percentages <- sweep(TrajCells, 1, RowTotal, "/") * 100
Percentages <- cbind(time = Trajectory_Optimized$time, Percentages)

# == Save as CSV
write.csv(Percentages, file.path(save_dir, "proportions", "Trajectory_Percentages.csv"), row.names = FALSE)

# == Melt 
Percentages_melted <- melt(Percentages, id.vars = "time", variable.name = "CellType", value.name = "Percentage")

# ====================================================================== Plots Using Percentages ======================================================================


# == Plot each cell type's proportion
prop_plot <- ggplot(Percentages_melted, aes(x = time, y = Percentage, color = CellType)) +
  geom_line(linewidth = 1) +  # Updated argument name
  labs(title = "Simulated Cell Type Proportions Over Time",
       x = "Time (DPI)", y = "Proportion of Total Cells") +
  theme_bw()

# == Save plot
ggsave(file.path(save_dir, "plots", "Trajectory_Percentage_Over_Time.png"), plot = prop_plot, width = 8, height = 5)



# == Line Plot: Cell Type Percentages Over Time
percent_line_plot <- ggplot(Percentages_melted, aes(x = time, y = Percentage, color = CellType)) +
  geom_line(linewidth = 1) +
  labs(title = "Simulated Cell Type Percentages Over Time",
       x = "Time (DPI)", y = "Percentage of Total Cells") +
  theme_bw()

ggsave(file.path(save_dir, "plots", "Trajectory_Percentages_Over_Time.png"),
       plot = percent_line_plot, width = 8, height = 5)
print(percent_line_plot)

# == Stacked Area Plot
percent_area_plot <- ggplot(Percentages_melted, aes(x = time, y = Percentage, fill = CellType)) +
  geom_area(alpha = 0.8) +
  labs(title = "Stacked Percentages of Simulated Cell Types",
       x = "Time (DPI)", y = "Percentage of Total Cells") +
  theme_bw()

ggsave(file.path(save_dir, "plots", "Sim_Stacked_Percentages.png"),
       plot = percent_area_plot, width = 8, height = 5)
print(percent_area_plot)

# == Faceted Bar Plot: Percentages Over Time
percent_bar_plot <- ggplot(Percentages_melted, aes(x = time, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~CellType, scales = "free_y") +
  labs(title = "Simulated Cell Type Percentages by Type",
       x = "Time (DPI)", y = "Percentage of Total Cells") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  scale_x_continuous(breaks = Percentages$time[seq(1, nrow(Percentages), by = 24)])

ggsave(file.path(save_dir, "plots", "Faceted_Percent_CellTypes.png"),
       plot = percent_bar_plot, width = 10, height = 6)
print(percent_bar_plot)

# == Faceted Bar Plot: Cell Type Percentages at Selected Time Points
key_times <- c(0, 1, 2, 3.5, 5, 7)

Percentages_melted_filtered <- Percentages_melted %>%
  filter(time %in% key_times)

time_facet_percent_plot_filtered <- ggplot(Percentages_melted_filtered, aes(x = CellType, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~time, scales = "free_y", ncol = 3) +
  labs(title = "Cell Type Percentages at Selected Time Points",
       x = "Cell Type", y = "Percentage of Total Cells") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(save_dir, "plots", "Faceted_Percentages_By_KeyTimepoints.png"),
       plot = time_facet_percent_plot_filtered, width = 10, height = 6)
print(time_facet_percent_plot_filtered)



