# =================================================
# Re-runing the same Al-Ghazawi code with the new 
# Re-annotated data included for simulation 
# =================================================


#################### Libraries ####################

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

# ===================================== Own Annotated Data (2025-07-20) ===================================================

# Start with the raw proportions (my own proportions)
McKellarPercent <- data.frame(
  time = c(0, 1, 2, 3.5, 5, 7),
  Md = c(NA, NA, NA, NA, NA, NA),
  N = c(0.0778, 0.175, 0.0195, 0.00663, 0.00443, 0.00755),
  Nd = c(NA, NA, NA, NA, NA, NA),
  M = c(0.194, 0.586, 0.185, 0.0516, 0.0348, 0.0878),
  M1 = c(0, 0.169, 0.719, 0.496, 0.0724, 0.0491),
  M2 = c(0.251, 0.0666, 0.0751, 0.379, 0.695, 0.560),
  QSC = c(0.477, 0.00150, 0.000917, 0.00414, 0.144, 0.228),
  ASC = c(0, 0.00307, 0.000678, 0.0407, 0.0287, 0.0295),
  Mc = c(0, 0, 0.000282, 0.0223, 0.0205, 0.0375)
)

# Convert to percentage & round (except 'time')
McKellarPercent <- McKellarPercent %>%
  mutate(across(-time, ~ round(.x * 100, 2)))

# McKellarPercent now holds the percentage values rounded to 2 decimal points
print(McKellarPercent)

# Now for DeMicheli 

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

# =========================================== Converting function to pre-defined scales ==========================================

# This function calibrates a given cells percent dataframe, like either of the two above
# to units of cells per cubic millimeter. It does so based on specified counts of
# representative cell types taken from the literature that should be achieved at different
# days post-injury.
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
    
    # Then apply the scale to the row
    CalibDF[row,] = RowScale * CalibDF[row,]
    
  }
  
  # Put but the times -- because we scaled them above, duh!
  CalibDF[,1] = Times
  return(CalibDF)
}

# Converted McKellar 
McKellarCellsPerMM3 <- CalibratePercents(McKellarPercent)
print(McKellarCellsPerMM3)

# COnverted Micheli 
DeMicheliCellsPerMM3 <- CalibratePercents(DeMicheliPercent)
print(DeMicheliCellsPerMM3)

# ============================================= ODE Model ==================================================

# Initial states of cell type variables
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

# Initial guesses at parameter values for differential equation model
# (arrived at by a combination of reasoning of time-scales for
# infiltration and exfiltration processes, as well as some amount of
# "by hand" tweaking based on qualitative expectations and comparison
# to McKellar data, but not DeMicheli data)
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

# ================================= ODE model dynamics ===============================================

# Observation times -- these are the McKellar days post injury, which
# conveniently include the DeMicheli days post injury as a subset
ObservationTimes <- c(0, 1, 2, 3.5, 5, 7)

# Function capturing regeneration dynamics
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

# ============================================== Simulation Function ======================================

# Function to perform one simulation trajectory for a given set of parameters
# A second parameter OutputTimes is optional, and allows user to specify
# which time points should be reported. By default, the function returns
# time points from 0 to 7 days in increments of 1 hour
MuscleSimulation <- function(Params,OutputTimes=c()) {
  
  # Figure out time points
  if (length(OutputTimes)==0) {
    OutputTimes = seq(from = 0, to = 7, by = 1/24)
  }
  
  # Call ODE simulator to solve dynamics
  Trajectory <- ode(y = InitialState,
                    times = OutputTimes,
                    func = MuscleRegenDerivs,
                    parms = Params,
                    method = 'lsoda')
  
  # Return solution
  return(Trajectory)
}

# ======================================= Plotting Function =========================================================

# Function to plot the outcome of a simulation, optionally showing also 
# McKellar and/or DeMicheli data
ShowSimulation <- function(TrajData,ShowMcK=FALSE,ShowDeM=FALSE) {
  # Rearrange data into long lists with time and variable name as keys, and cells/mm3 values in Traj, McK or DeM columns
  TrajDF <- as.data.frame(TrajData)
  TrajLong <- TrajDF %>% 
    pivot_longer(cols = -time, names_to = "variable", values_to = "Traj")
  McKLong <- McKellarCellsPerMM3 %>% 
    pivot_longer(cols = -time, names_to = "variable", values_to = "McK")
  DeMLong <- DeMicheliCellsPerMM3 %>% 
    pivot_longer(cols = -time, names_to = "variable", values_to = "DeM")
  
  # Then join those, essentially filling in NA for McK and DeM wherever they don't have a value
  Joined_Temp <- left_join(TrajLong, McKLong, by = c("time", "variable")) 
  Joined <- left_join(Joined_Temp, DeMLong, by = c("time", "variable"))
  
  # Finally, separate them again!
  
  # Create plots for each variable
  plots <- lapply(unique(Joined$variable), function(var) {
    
    # Extract trajectory and observed just for that variable
    JoinedVar = filter(Joined,variable==var)
    
    # Plot it, depending on which observed data was requested
    if (ShowMcK==FALSE & ShowDeM==FALSE) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    }
    if (ShowMcK==TRUE & ShowDeM==FALSE) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + geom_point(aes(y=McK), color="blue") + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    }
    if (ShowMcK==FALSE & ShowDeM==TRUE) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + geom_point(aes(y=DeM), color="red") + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw()
    }
    if (ShowMcK==TRUE & ShowDeM==TRUE) {
      ggp <- ggplot(JoinedVar, aes(x=time)) + geom_line(aes(y=Traj)) + geom_point(aes(y=McK), color="blue", alpha=0.5) + geom_point(aes(y=DeM), color="red", alpha=0.5) + labs(x = "time (DPI)", y = paste(var, "(cells/mm3)")) + theme_bw() + theme(text = element_text(size=7))
    }
    return(ggp)
  })
  
  # Combine all the plots into one grid
  plot_grid <- do.call(grid.arrange, c(plots, ncol = 3, nrow = 3))
  
  # Save the combined plot grid to the specified location
  ggsave("D:/Masters/Processing Codes/Masters Paper Figure/ODE/NoUse_simulation_results.png", 
         plot = plot_grid, 
         width = 12, height = 12, dpi = 300)
  
  # Print the combined plot grid
  print(plot_grid)
}

# How does it look with default parameters?
Traj <- MuscleSimulation(ParamsInit)
ShowSimulation(Traj,TRUE,FALSE)

########################################## Model fitting-related functions ##########################################

# Objective function to minimize - sum squared error of residuals of each
# state variable, normalized by the maximum observed value of that state variable.
# By default, the Params are evaluated with respect to McKellar dataset.
# Enter
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
  
  # Calculate error for each variable
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

# Similar to the above function Error with parameter "McKellar" for dataset.
# However, this function assumes that the logarithms of the parameters are input.
# It takes the exponential of that input to compute the regular-scale parameters,
# which it then feeds into Error. This is useful because it ensures regular-scale
# parameters are always positive.
Error_LogParams <- function(LogParams,CompareTo="McKellar") {
  Params <- exp(LogParams)
  return(Error(Params,CompareTo))
}

# Error of initial parameters
McK_ErrorInit <- Error(ParamsInit,"McKellar")
print(paste("Error on McKellar data with initial parameters:",McK_ErrorInit))
DeM_ErrorInit <- Error(ParamsInit,"DeMicheli")
print(paste("Error on DeMicheli data with initial parameters:",DeM_ErrorInit))

# Parameter optimization by Nelder-Mead method
y = optim(log(ParamsInit),Error_LogParams,method="Nelder-Mead")
ParamsOpt <- exp(y$par)
TrajOpt <- MuscleSimulation(ParamsOpt)
ShowSimulation(TrajOpt,TRUE,TRUE)

# Error after optimization?
McK_ErrorOpt <- Error(ParamsOpt,"McKellar")
print(paste("Error on McKellar data with optimized parameters:",McK_ErrorOpt))
DeM_ErrorOpt <- Error(ParamsOpt,"DeMicheli")
print(paste("Error on DeMicheli data with optimized parameters:",DeM_ErrorOpt))

# ================================= Initial and Optimized Parameter table for the paper ======================

# Print and save initial and optimized parameters
cat("Initial Parameters (ParamsInit):\n")
print(ParamsInit)

cat("\nOptimized Parameters (ParamsOpt):\n")
print(ParamsOpt)

# Save both to CSV for documentation
ParamTable <- data.frame(
  Parameter = names(ParamsInit),
  Initial = round(ParamsInit, 6),
  Optimized = round(ParamsOpt, 6)
)

# Set output path
output_dir <- "D:/Masters/Processing Codes/Masters Paper Figure/"
output_file <- file.path(output_dir, "ODE_Model_Parameter_Comparison.csv")

# Write CSV
write.csv(ParamTable, output_file, row.names = FALSE)

cat(paste("\nParameter table saved to:", output_file, "\n"))

# Save both to CSV for documentation — with better small value preservation
ParamTable <- data.frame(
  Parameter = names(ParamsInit),
  Initial = format(ParamsInit, scientific = TRUE, digits = 6),
  Optimized = format(ParamsOpt, scientific = TRUE, digits = 6)
)

# Output location
output_dir <- "D:/Masters/Processing Codes/Masters Paper Figure/"
output_file <- file.path(output_dir, "Accurate_ODE_Model_Parameter_Comparison.csv")

# Save to CSV
write.csv(ParamTable, output_file, row.names = FALSE)

cat(paste("\nParameter table saved to:", output_file, "\n"))


############################
### Sensitivity analysis ###
############################

# First, we compute the standard trajectory at the optimized parameters
# at 1 hour time interval
TrajOpt = as.data.frame(MuscleSimulation(ParamsOpt))

# Now, for each parameter, we simulate what would happen when that parameter
# is increased by 1%
TrajsPerturb <- c()
for (ParamIndex in 1:length(ParamsOpt)) {
  ParamsPerturb <- ParamsOpt
  ParamsPerturb[ParamIndex] <- ParamsPerturb[ParamIndex]*1.01
  TrajPerturb <- as.data.frame(MuscleSimulation(ParamsPerturb))
  TrajsPerturb <- c(TrajsPerturb,TrajPerturb)
}

# A function for plotting Sensitivity results after grouping by variable, 
# which are computed below
OurBluRed <- c('#6666FF','#7777FF','#8888FF','#9999FF','#AAAAFF','#BBBBFF','#CCCCFF','#DDDDFF','#EEEEFF','#FFFFFF','#FFEEEE','#FFDDDD','#FFCCCC','#FFBBBB','#FFAAAA','#FF9999','#FF8888','#FF7777','#FF6666')
ShowSensitivity <- function(PerturbMat,VarName='') {
  MaxAbs <- max(abs(PerturbMat))
  ToShow <- 0.5*PerturbMat/MaxAbs + 0.5 # Scaling between zero and one, with zero in the center
  NRows <- nrow(ToShow)
  NCols <- ncol(ToShow)
  for (i in 1:NRows) {
    ToShow[i,1:NCols] <- ToShow[i,NCols:1]
  }
  image(ToShow,col=OurBluRed,zlim=c(0,1),xaxt="n",yaxt="n")
  title(VarName)
  angleAxis(side=4, at=((0:(length(ParamsInit)-1))-0.0)/(length(ParamsInit)-1), labels=rev(names(ParamsInit)), srt=0, cex=0.75, xpd=TRUE)
  axis(side=1, at=(0:7)/7, labels=c('0','1','2','3','4','5','6','7'), cex=0.8, lwd.ticks=0, padj=-1.6)
}

# Now, for each state variable, compute a matrix of the proportionate differences
VarNames = names(InitialState)

# Save the grid of heatmaps as a PNG file
png(filename = "D:/Masters/Processing Codes/Masters Paper Figure/ODE/sensitivity_analysis_grid.png", 
    width = 2400, height = 2400, res = 300)  # Adjust the size as needed

# Set up a 3x3 grid of plots for the heatmaps
par(mfrow=c(3,3),mai=c(0.15,0.1,0.4,0.4))

# Loop through the variables and generate the sensitivity heatmaps
for (StateIndex in 1:length(InitialState)) {
  VarMax <- max(TrajOpt[[StateIndex+1]])
  PerturbMat <- matrix(data=0,nrow=length(TrajsPerturb[[1]]),ncol=length(ParamsOpt))
  for (ParamIndex in 1:length(ParamsOpt)) {
    PerturbMat[,ParamIndex] <- (TrajsPerturb[[10*(ParamIndex-1)+StateIndex+1]]-TrajOpt[[StateIndex+1]])/VarMax
  }
  PerturbMaxAbs <- max(abs(PerturbMat))
  PerturbMat <- PerturbMat / PerturbMaxAbs
  ShowSensitivity(PerturbMat,VarNames[StateIndex])
} 

# Close the PNG device to save the image
dev.off()
