#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Participant 1 regression
# Description: Functions for Experiment 3 analysis - P1      
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 8/02/2016                                                                            
# Modifications:
#-----------------------------------------------------------------------------------------------------
curve.fitting.1 <- function(i, p.agg, ap.agg, p.glm, ap.glm){
  #--------------------------------------------
  # Plot psychometric curves and extract threshold and slope values
  #--------------------------------------------
  opar <- par(mfrow = c(1,2), lwd = 3)
  out1 <- periodicity_psy_curves(i, p.agg, ap.agg, p.glm, ap.glm)
  out2 <- bootstrap.periodicity.testLevel(i, p.agg, ap.agg, p.glm, ap.glm)
  par(opar)
  #--------------------------------------------
  # Return fitted threshold and slope values
  #--------------------------------------------
  return.list <- list(out1, out2)
  return(return.list)
}