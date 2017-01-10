#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Participant specific model fit and data generation
# Description: Functions for Experiment 3 analysis     
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 8/02/2016                                                                            
# Modifications: 11/11/2016
#-----------------------------------------------------------------------------------------------------
participant <- function(dat, i){
  #--------------------------------------------
  # Load participant data
  #--------------------------------------------
  #tmp <- get.cleansubset(dat, i)
  tmp <- subset(dat, dat$subNo == i);
  #--------------------------------------------
  # Periodicity data subset
  #--------------------------------------------
  pdat <- subset(tmp, tmp$periodicity == "Periodic")
  apdat <- subset(tmp, tmp$periodicity == "Aperiodic")
  
  #--------------------------------------------
  # Aggregated periodicity tables
  #--------------------------------------------
  pdat.agg <- data.frame(ddply(pdat,~periodicity + testLevel,summarise,Right=sum(answerRight),Left=sum(answerLeft)))
  apdat.agg <- data.frame(ddply(apdat,~periodicity + testLevel,summarise,Right=sum(answerRight),Left=sum(answerLeft)))
  
  #--------------------------------------------
  # Periodic - Select the best fitting GLM sigmoid function
  #--------------------------------------------
  pdat.agg <- data.frame(ddply(pdat,~periodicity * testLevel,summarise,Right=sum(answerRight),Left=sum(answerLeft)))
  
  #--------------------------------------------
  # Periodic - refit model
  #--------------------------------------------
  v <- list()
  v[["logit"]] <- glm(cbind(Left, Right) ~ testLevel, family = binomial, data = pdat.agg)
  v[["probit"]] <- glm(cbind(Left, Right) ~ testLevel, family = binomial(probit), data = pdat.agg)
  v[["cauchit"]] <- glm(cbind(Left, Right) ~ testLevel,  family = binomial(cauchit), data = pdat.agg)
  v[["weibull"]] <- glm(cbind(Left, Right) ~ testLevel,  family = binomial(cloglog), data = pdat.agg)
  aicScores <- sapply(v, AIC)
  print(aicScores)
  p.glm <- v[[which(aicScores == min(aicScores), arr.ind = TRUE)]]
  print(summary(p.glm))
  
  #--------------------------------------------
  # Aperiodic - remove outliers
  #--------------------------------------------
  apdat.agg <- data.frame(ddply(apdat,~periodicity * testLevel,summarise,Right=sum(answerRight),Left=sum(answerLeft)))
  
  #--------------------------------------------
  # Aperiodic - refit model
  #--------------------------------------------
  v <- list()
  v[["logit"]] <- glm(cbind(Left, Right) ~ testLevel, family = binomial, data = apdat.agg)
  v[["probit"]] <- glm(cbind(Left, Right) ~ testLevel, family = binomial(probit), data = apdat.agg)
  v[["cauchit"]] <- glm(cbind(Left, Right) ~ testLevel,  family = binomial(cauchit), data = apdat.agg)
  v[["weibull"]] <- glm(cbind(Left, Right) ~ testLevel,  family = binomial(cloglog), data = apdat.agg)
  aicScores <- sapply(v, AIC)
  print(aicScores)
  ap.glm <- v[[which(aicScores == min(aicScores), arr.ind = TRUE)]]
  print(summary(ap.glm))
  
  #--------------------------------------------
  # Return aggregated data and fitted glm models
  #--------------------------------------------
  return.list <- list(i, pdat.agg, apdat.agg, p.glm, ap.glm)
  return(return.list)
}