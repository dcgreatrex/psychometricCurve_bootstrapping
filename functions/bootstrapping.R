#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Periodicity vs testLevel Bootstrapping tables and curves
# Description: Takes an aggregated table as input for both periodicity conditions, plots the best
# fitting psychometric curves, computes and returns threshold and slope values
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 8/02/2016                                                                            
# Modifications: 11/11/2016
#-----------------------------------------------------------------------------------------------------
fitBootstrap <- function(i, p.agg, ap.agg, p.glm, ap.glm){
  #--------------------------------------------
  # Run boostrapping procedures to estimate model fit
  #--------------------------------------------
  library(boot)
  
  psyfun.stat <- function(d){
    nn <- rowSums(d[,1:2])
    t.glm <- glm(cbind(d[,4], nn - d[,4]) ~ d[, 3], binomial("probit"))
    as.vector(coef(t.glm))
  }
  psyfun.gen <- function(d, mle){
    nn <- with(d, Right + Left)
    d$resp <- rbinom(nrow(d), nn, mle)
    d
  }
  
  # periodic bootstrapping
  p.glm.d <- cbind(p.agg, resp = p.agg$Left)
  p.glm.d[,1] <- NULL
  p.glm.d <- p.glm.d[c(3,2,1,4)]
  p.boot <- boot(p.glm.d, statistic = psyfun.stat,
                 R = 10000, sim = "parametric",
                 ran.gen = psyfun.gen,
                 mle = fitted(p.glm))
  p.boot
  # plot periodic fits
  #plot(p.boot, index = 1)
  #plot(p.boot, index = 2)
  
  # periodic bootstrapping
  ap.glm.d <- cbind(ap.agg, resp = ap.agg$Left)
  ap.glm.d[,1] <- NULL
  ap.glm.d <- ap.glm.d[c(3,2,1,4)]
  ap.boot <- boot(ap.glm.d, statistic = psyfun.stat,
                  R = 10000, sim = "parametric",
                  ran.gen = psyfun.gen,
                  mle = fitted(ap.glm))
  ap.boot
  # plot aperiodic fits
  #plot(ap.boot, index = 1)
  #plot(ap.boot, index = 2)
  
  # extract mean threshold and slope values from each bootstrap
  # Periodic
  boot.model <- p.boot
  bootstrap.out1 <- cbind(summary(boot.model)[1,c(2:4)], mean(boot.model$t[,1]), summary(boot.model)[2,c(2:4)], mean(boot.model$t[,2]))
  # calculate 95% confidence intervals
  boot.ci1 <- boot.ci(boot.model, type =c("norm"), index = 1)
  boot.ci2 <- boot.ci(boot.model, type =c("norm"), index = 2)
  bootstrap.out1 <- cbind(bootstrap.out1, boot.ci1$normal[,2], boot.ci1$normal[,3], boot.ci2$normal[,2], boot.ci2$normal[,3])
  colnames(bootstrap.out1) <- c("B0.original", "B0.bootBias", "B0.bootSE", "B0.bootMEAN", 
                                "B1.original", "B1.bootBias",  "B1.bootSE", "B1.bootMEAN",
                                "thresh.95%c.int.1", "thresh.95%c.int.2", "slope.95%c.int.1", "slope.95%c.int.2")
  bootstrap.out1
  
  # Aperiodic
  boot.model <- ap.boot
  bootstrap.out2 <- cbind(summary(boot.model)[1,c(2:4)], mean(boot.model$t[,1]), summary(boot.model)[2,c(2:4)], mean(boot.model$t[,2]))
  boot.ci1 <- boot.ci(boot.model, type =c("norm"), index = 1)
  boot.ci2 <- boot.ci(boot.model, type =c("norm"), index = 2)
  bootstrap.out2 <- cbind(bootstrap.out2, boot.ci1$normal[,2], boot.ci1$normal[,3], boot.ci2$normal[,2], boot.ci2$normal[,3])
  colnames(bootstrap.out2) <- c("B0.original", "B0.bootBias", "B0.bootSE", "B0.bootMEAN", 
                                "B1.original", "B1.bootBias",  "B1.bootSE", "B1.bootMEAN",
                                "thresh.95%c.int.1", "thresh.95%c.int.2", "slope.95%c.int.1", "slope.95%c.int.2")
  bootstrap.out2
  
  # final bootstrap table
  bootstrap.out <- rbind(bootstrap.out1, bootstrap.out2)
  rm(bootstrap.out1, bootstrap.out2)
  rownames(bootstrap.out) <- c("periodic", "aperiodic")
  bootstrap.out
  
  # bootstrap threshold slope table
  p.thresh <- -bootstrap.out[1,4]/bootstrap.out[1,8]
  p.slope <- bootstrap.out[1,8]
  ap.thresh <- -bootstrap.out[2,4]/bootstrap.out[2,8]
  ap.slope <- bootstrap.out[2,8]
  boot.t.sl.out <- rbind(c(p.thresh, p.slope), c(ap.thresh, ap.slope))
  colnames(boot.t.sl.out) <- c("boot.thresh", "boot.slope")
  rownames(boot.t.sl.out) <- c("periodic", "aperiodic")
  boot.t.sl.out
  
  #--------------------------------------------
  # Plot Bootstrapped psychometric curves
  #--------------------------------------------
  # periodic bootstrapping
  p.boot <- boot(p.glm.d, statistic = psyfun.stat,
                 R = 2000, sim = "parametric",
                 ran.gen = psyfun.gen,
                 mle = fitted(p.glm))
  # aperiodic bootstrapping
  ap.boot <- boot(ap.glm.d, statistic = psyfun.stat,
                  R = 2000, sim = "parametric",
                  ran.gen = psyfun.gen,
                  mle = fitted(ap.glm))
  # plot
  cxx <- seq(min(p.agg$testLevel), max(p.agg$testLevel), len = 1000)
  cc.bt <- p.boot$t
  with(p.agg, plot(testLevel, Left/(Left+Right), type = "n", xaxt="n", yaxt="n", ylab = NA, xlab = NA, xaxlabels = "n", yaxlabels = "n", 
                   main = paste("Bootstrap (n=2000) psychometric curves\nParticipant", i), ylim=c(0, 1)))
  axis(1,cex.axis=1.25)
  axis(2,cex.axis=1.25)
  mtext(1, text = "\nDisplacement of test tone away from average mid point
        Interaural intensity difference dB", 4, cex = 1.15)
  mtext(2, text = "% 'left' responses", 3, cex = 1.25)
  abline(v = 0, lty = 2, col = "grey")
  invisible(apply(cc.bt, 1, function(x){
    lines(cxx, pnorm(cxx,-(x[1]/x[2]), 1/x[2]), 
          col = rgb(0,0,0,0.01))
  }))
  invisible(apply(ap.boot$t, 1, function(x){
    lines(cxx, pnorm(cxx,-(x[1]/x[2]), 1/x[2]), 
          col = rgb(0,0.8,0,0.01))
  }))
  points(I(Left/(Right+Left)) ~ testLevel, p.agg, pch = 21, bg = "black", cex = 2, col = "grey")
  points(I(Left/(Right+Left)) ~ testLevel, ap.agg, pch = 21, bg = "green3", cex = 2, col = "black")
  # add legend
  legend(min(cxx),1.02, c('Periodic','Aperiodic'),pch = 19,
         col=c('black','green3'), bty = "n", cex=1)
  
  #--------------------------------------------
  # Return fitted threshold and slope values
  #--------------------------------------------
  return.list <- list(bootstrap.out, boot.t.sl.out)
  return(return.list)
}