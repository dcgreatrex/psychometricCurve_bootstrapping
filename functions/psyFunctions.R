#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Periodicity vs testLevel psychometric curves
# Description: Takes an aggregated table as input for both periodicity conditions, plots the best
# fitting psychometric curves, computes and returns threshold and slope values
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 8/02/2016                                                                            
# Modifications: 11/11/2016
#-----------------------------------------------------------------------------------------------------
fitPsychfunctions <- function(i, p.agg, ap.agg, p.glm, ap.glm){
  
  #--------------------------------------------
  # Periodic - predict responses using model
  #--------------------------------------------
  p.mns <- with(p.agg, tapply(Left/(Left+Right), testLevel, mean))
  testLevel.p <- seq(min(p.agg$testLevel), max(p.agg$testLevel), len = 1000)
  p.pred <- data.frame(testLevel = testLevel.p)
  p.pred$y <- predict(p.glm, p.pred, type="response")
  
  #--------------------------------------------
  # Aperiodic - predict responses using model
  #--------------------------------------------
  ap.mns <- with(ap.agg, tapply(Left/(Left+Right), testLevel, mean))
  testLevel.ap <- seq(min(ap.agg$testLevel), max(ap.agg$testLevel), len = 1000)
  ap.pred <- data.frame(testLevel = testLevel.ap)
  ap.pred$y <- predict(ap.glm, ap.pred, type="response")
  
  #--------------------------------------------
  # Plot psychometric curves
  #--------------------------------------------
  plot(as.numeric(names(p.mns)), p.mns, col = 'black', xaxt="n", yaxt="n", ylab = NA, xlab = NA, xaxlabels = "n", yaxlabels = "n",
       main = paste("Psychometric curves \nParticipant", i), ylim=c(0, 1))
  axis(1,cex.axis=1.25)
  axis(2,cex.axis=1.25)
  mtext(1, text = "\nDisplacement of test tone away from average mid point
        Interaural intensity difference dB", 4, cex = 1.15)
  mtext(2, text = "% 'left' responses", 3, cex = 1.25)
  abline(v = 0, lty = 2, col = "grey")
  # Periodici standard errors
  confid.p <- predict(p.glm, newdata = list(testLevel = testLevel.p), type = "response", se.fit = TRUE)
  polygon(c(testLevel.p, rev(testLevel.p)), 
          with(confid.p, c(fit + 1 * se.fit, rev(fit - 1 * se.fit))), 
          col = alpha("grey", 0.5), border = "white", lty = 1)
  # Aperiodici standard errors
  confid.ap <- predict(ap.glm, newdata = list(testLevel = testLevel.ap), type = "response", se.fit = TRUE)
  polygon(c(testLevel.ap, rev(testLevel.ap)), 
          with(confid.ap, c(fit + 1 * se.fit, rev(fit - 1 * se.fit))), 
          col = alpha("green3", 0.5), border = "white", lty = 1)
  # add lines
  lines(testLevel.p, p.pred$y, col = "black", lwd = 8)
  lines(testLevel.ap, ap.pred$y, col = "green3", lwd = 8)
  # add points
  points(as.numeric(names(p.mns)), p.mns, col = 'black', pch = 21, bg = "black", cex = 2)
  points(as.numeric(names(ap.mns)), ap.mns, col = 'black', pch = 21, bg = "green3", cex = 2)
  # add legend
  legend(min(testLevel.p),1.02, c('Periodic','Aperiodic'),pch = 19,
         col=c('black','green3'), bty = "n", cex=1.25)
  
  
  #--------------------------------------------
  # Exptract threshold and slope values
  #--------------------------------------------
  # Gold, J., Nodal, F., Peters, F., King, A., & Bajo, V. (2015). Auditory gap-in noise detection in ferrets and humans. Behavioural Neuroscience. 129(4), 473-490.
  # http://www.dlinares.org/psychopract.html
  
  # threshold and slope function
  thresholdslope <- function(model, thresh){
    mean <- -coef(model)[[1]] / coef(model)[[2]] 
    sd <- abs(1 / coef(model)[[2]])
    # threshold
    threshold <- qnorm(thresh, mean, abs(sd))
    # slope
    a <- qnorm(p = (thresh - 0.001), mean, sd)
    b <- qnorm(p = (thresh + 0.001), mean, sd)
    slope <- ((thresh - 0.001) - (thresh + 0.001))/(b-a)
    return(c(threshold, slope))
  }
  
  threshold.out = c(); slope.out = c()
  for (i in 1:2){
    if (i == 1){
      model = p.glm
      mf <- threshold_slope(p.pred$y, testLevel.p)
      mf <- c(mf$x_th, mf$slope)
    }else{
      model = ap.glm
      mf <- threshold_slope(ap.pred$y, testLevel.ap)
      mf <- c(mf$x_th, mf$slope)
    }
    dg.out <- thresholdslope(model, 0.5)
    # glm fit threshold_slope
    glm.out <- c((-coef(model)[[1]] / coef(model)[[2]]), coef(model)[[2]])
    # arrange as array for output
    threshold.out <- rbind(threshold.out, c(mf[1], dg.out[1], glm.out[1]))
    slope.out <- rbind(slope.out, c(mf[2], dg.out[2], glm.out[2]))
  }
  
  #--------------------------------------------
  # rename threshold and slope tables
  #--------------------------------------------
  colnames(threshold.out) <- c("th.modelfree", "th.dg.estimate", "th.glm.estimate")
  rownames(threshold.out) <- c("periodic", "aperiodic")
  colnames(slope.out) <- c("sl.modelfree", "sl.dg.estimate", "sl.glm.estimate")
  rownames(threshold.out) <- c("periodic", "aperiodic")
  thresh.sl.out <- cbind(threshold.out,slope.out)
  
  #--------------------------------------------
  # return threshold and slope values as a list
  #--------------------------------------------
  return(thresh.sl.out)
}