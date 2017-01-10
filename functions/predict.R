#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Participant 1 regression
# Description: Functions for Experiment 3 analysis - P1      
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 8/02/2016                                                                            
# Modifications: 11/11/2016
#-----------------------------------------------------------------------------------------------------
curve.fitting <- function(i, p.agg, ap.agg, p.glm, ap.glm){
  #--------------------------------------------
  # Plot psychometric curves and extract threshold and slope values
  #--------------------------------------------
  opar <- par(mfrow = c(1,2), lwd = 3)
  out1 <- fitPsychfunctions(i, p.agg, ap.agg, p.glm, ap.glm)
  out2 <- fitBootstrap(i, p.agg, ap.agg, p.glm, ap.glm)
  par(opar)
  #--------------------------------------------
  # Return fitted threshold and slope values
  #--------------------------------------------
  return.list <- list(out1, out2)
  return(return.list)
}