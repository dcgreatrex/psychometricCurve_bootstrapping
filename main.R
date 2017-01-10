#-----------------------------------------------------------------------------------------------------
# Author: David Greatrex, PhD Candidate, University of Cambridge.  
# Date: 01/03/2016, 11/11/2016
# Name: Fitting psychometric functions using a bootstrapping procedure
#-----------------------------------------------------------------------------------------------------
library(lattice)
library(leaps)
library(scales)
library(boot)
library(psyphy)
library(modelfree)
library(base)
library(ggplot2)

#--------------------------------------------
# User input variables - define the path to the location of the downloaded folder
#--------------------------------------------
mainWd = "/Users/..../psychometricCurve_bootstrapping"

#--------------------------------------------
# Set working directories and fixed variables
#--------------------------------------------
# working directories and file names
functionsWd = paste(mainWd, "functions", sep = "/")
plotsWd = paste(mainWd, "plots", sep = "/")
dataWd <- paste(mainWd, "data/participantData", sep = "/")
tablesWd <- paste(mainWd, "data/tables", sep = "/")
# number of participant in the experiment  
noPs <- 24

#--------------------------------------------
#- load and list functions
#--------------------------------------------
# 1) main Expt 3 function
source(paste(functionsWd,"prepare.r",sep="/")) # load function script
ListSegmentationFunctions()
# 2) participant specific data aggregation and model fit
source(paste(functionsWd,"glmfit.r",sep="/")) # load function script
# 3) peridoicity vs testLevel aggregation tables
source(paste(functionsWd,"predict.r",sep="/")) # load function script
# 3.1) periodicity vs testLevel curve fitting function (called from within curve_fit_1)
source(paste(functionsWd,"psyFunctions.r",sep="/")) # load function script
# 3.2) periodicity vs testLevel bootstrapping function (called from within curve_fit_1)
source(paste(functionsWd,"bootstrapping.r",sep="/")) # load function script

#--------------------------------------------
#- load and clean data
#--------------------------------------------
datMaster <- loadData(dataWd, noPs) 
dat <- cleanData(datMaster)
dat <- removeLongTrials(dat)
dat <- renameData(dat)
dat <- add.sdtflags(dat)
dat$answerLeft <- ifelse(dat$answer == -1, 1, 0)
dat$answerRight <- ifelse(dat$answer == 1, 1, 0)

#--------------------------------------------
# 1 - Fit psychometric functions to both experimental conditions for each participant and extract threshold and slope values
#--------------------------------------------
S1.outputFolder <- "PsychometricFits"
dir.create(paste(plotsWd, S1.outputFolder, sep = "/"))
out.stage1 <- data.frame()
for (i in 1:noPs){
  # set plot features
  plotTitle <- paste0("P", i, "_periodicity_curves.png")
  png(paste(plotsWd, S1.outputFolder, plotTitle, sep = "/"), width = 800, height = 400)
  # fit P i with glm models and remove outliers
  p <- participant(dat, i)
  # plot P i with psychometric curves and bootstrap
  p.out <- curve.fitting(p[1], data.frame(p[2]), data.frame(p[3]), p[[4]], p[[5]])
  dev.off()
  rm(p)
  # process output so that it can be appended to output data.frame
  out.tmp <- cbind(p.out[[1]], p.out[[2]][[1]], p.out[[2]][[2]])
  rm(p.out)
  # define descriptor columns for output dataset
  desc <- matrix(c(i, i, "periodic", "aperiodic"), nrow = 2, ncol = 2)
  colnames(desc) <- c("subNo", "periodicity")
  # append P i outout to the out.stage1 data.frame
  out.stage1 <- rbind(out.stage1,  data.frame(cbind(desc, out.tmp)))
  rm(out.tmp, desc)
}
# save the out.stage1 table as a csv file for further reference
rownames(out.stage1) <- NULL
write.csv(out.stage1, file = paste(tablesWd, "stage1Table.csv", sep = "/"))

#--------------------------------------------
# 2 - Plot group slope and threshold values
#--------------------------------------------
# open file to export plots
plotTitle <- "threhsold and slope values.png"
png(paste(plotsWd, plotTitle, sep = "/"), width = 400, height = 800)
opar <- par(mfrow = c(2,1), lwd = 3)

# remove outlying participants from the out.stage1 file - subjects 15 and 18 to be removed
out.stage1 <- subset(out.stage1, (out.stage1$subNo != 15) & (out.stage1$subNo != 18))

# threshold
p.thresh <- out.stage1$boot.thresh[out.stage1$periodicity=="periodic"]
ap.thresh <- out.stage1$boot.thresh[out.stage1$periodicity=="aperiodic"]
thresh <- t(rbind(p.thresh, ap.thresh))
# plot threshold values for each for the periodicity conditions
plot(thresh[,1], thresh[,2], xlim = c(-2.3,2.3), ylim = c(-2.8,2.8), pch = 19, xlab = NA, ylab = NA, cex.main = 2,
     main = "Threshold",
     axes=FALSE, frame.plot=TRUE)
axis(1,cex.axis=1.25)
axis(2,cex.axis=1.25)
mtext(1, text = "Periodic", 3, cex = 1.5, font = 2)
mtext(2, text = "Aperiodic", 3, cex = 1.5, font = 2)
abline(0, 1, lty = 2)
# add average threshold to plot
points(mean(thresh[,1]), mean(thresh[,2]), col = "dodgerblue", pch = 19, cex = 2)
points(mean(thresh[,1]), mean(thresh[,2]), col = "dodgerblue", pch = 3, cex = 4)
rm(p.thresh, ap.thresh)

# slope
p.slope <- out.stage1$boot.slope[out.stage1$periodicity=="periodic"]
ap.slope <- out.stage1$boot.slope[out.stage1$periodicity=="aperiodic"]
slope <- t(rbind(p.slope, ap.slope))
# plot threshold values for each for the periodicity conditions
plot(slope[,1], slope[,2], xlim = c(0.1,0.6), ylim = c(0.1,0.6), pch = 19, xlab = NA, ylab = NA, cex.main = 2,
     main = "Slope",
     axes=FALSE, frame.plot=TRUE)
axis(1,cex.axis=1.25)
axis(2,cex.axis=1.25)
mtext(1, text = "Periodic", 3, cex = 1.5, font = 2)
mtext(2, text = "Aperiodic", 3, cex = 1.5, font = 2)
abline(0, 1, lty = 2)
# add average slope to plot
points(mean(slope[,1]), mean(slope[,2]), col = "dodgerblue", pch = 19, cex = 2)
points(mean(slope[,1]), mean(slope[,2]), col = "dodgerblue", pch = 3, cex = 4)
rm(p.slope, ap.slope)

# close write to file
par(opar)
dev.off()

#--------------------------------------------
# run one-way t-tests to confirm that B1 (slope) is positive across both the periodic and aperiodic conditions
#--------------------------------------------
# periodic slope
mean(slope[,1]) # 0.2700701
sd(slope[,1]) #  0.0840744
t.test(slope[,1],mu= 0)
# OUTPUT
# One Sample t-test
# 
# data:  slope[, 1]
# t = 15.067, df = 21, p-value = 9.874e-13
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.2327936 0.3073466
# sample estimates:
#   mean of x 
# 0.2700701 

# aperiodic slope
mean(slope[,2]) # 0.2836022
sd(slope[,2]) #  0.09881552
t.test(slope[,2],mu= 0)
# OUTPUT
# One Sample t-test
# 
# data:  slope[, 2]
# t = 13.462, df = 21, p-value = 8.503e-12
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.2397899 0.3274146
# sample estimates:
#   mean of x 
# 0.2836022 

#--------------------------------------------
# run t.tests to compare threshold and slope values across the group
#--------------------------------------------
# Threshold
mean(thresh[,1]); mean(thresh[,2]) # -0.06566494 # 0.130707
sd(thresh[,1]); sd(thresh[,2]) # 1.442921 # 1.672011
t.test(thresh[,1], thresh[,2])
# OUTPUT
# Welch Two Sample t-test
# 
# data:  thresh[, 1] and thresh[, 2]
# t = -0.41705, df = 41.12, p-value = 0.6788
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.1472128  0.7544688
# sample estimates:
#   mean of x   mean of y 
# -0.06566494  0.13070704 

# Slope
t.test(slope[,1], slope[,2])
# OUTPUT
# Welch Two Sample t-test
# 
# data:  slope[, 1] and slope[, 2]
# t = -0.48921, df = 40.95, p-value = 0.6273
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.06939697  0.04233269
# sample estimates:
#   mean of x mean of y 
# 0.2700701 0.2836022  

#--------------------------------------------
rm(thresh, slope)
#--------------------------------------------