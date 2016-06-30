#-----------------------------------------------------------------------------------------------------
# Title: PhD experimental analysis - Experiment 3
# Description: Functions for Experiment 3 analysis      
# Author: David Greatrex, PhD Candidate, University of Cambridge. 
# Author responsibility: PhD candidate
# Date: 13/01/2016                                                                            
# Modifications:
#-----------------------------------------------------------------------------------------------------

#-------------------------
#- load depenent libraries
#-------------------------
library("data.table")
library("plyr")
library("ggplot2")
library("lme4")
library("lattice")
library("sqldf")
library("reshape")
library("car")
library("grid")
library("ggthemes") 
library("reshape2")
library("sqldf")
options(sqldf.driver = "SQLite")

#=================================================================================
# Print a list of all functions within this file to screen
#=================================================================================
ListSegmentationFunctions <- function(){
  print("Title: PhD Experiment 2 Analytical Functions",quote=FALSE)
  print("Author: David Greatrex, PhD Candidate, University of Cambridge.",quote=FALSE)
  print("This script contains the following functions:",quote=FALSE)
  print("loadData() : loads a master data file containing appended participant data",quote=FALSE)
  print("cleanData() : removes NAs, and trials determined to be invalid",quote=FALSE)
  print("renameData() : recodes undesirable data encoding in original data",quote=FALSE)
  print("add.sdtflags() : add signal detection theory variables to data set and return flagged data",quote=FALSE)
  print("get.cleansubset() : get subset of data for participant i (no rts > ±3sd or < 0.250 )",quote=FALSE)
  print("get.rtdata() : Clean and tag RT data for analysis",quote=FALSE)
}

loadData <- function(wd, noPs){
  #----------------------------
  # Load each participants data, append to a master data file and return the master
  #----------------------------
  for (i in 1:noPs){
    # Load participant data
    file = paste(toString(i), "/", toString(i), ".txt", sep="")  
    datTmp = read.table(paste(wd,file, sep='/'), header=TRUE)
    # append to a master data file
    if (i == 1){
      dat <- datTmp  
    }else{
      dat <- rbind(dat,datTmp)
    }
  }
  rm(datTmp)
  return(dat)
}

cleanData <- function(dat){
  #----------------------------
  # remove all NAs
  #----------------------------
  datClean <- subset(dat, dat$answer != 0)
  return(datClean)
}

removeLongTrials <- function(dat){
  dat2 <- data.frame()
  for (i in 1:max(dat$subNo)){
    tmp <- subset(dat, dat$subNo == i)
    sd <- sd(tmp$rt)
    mu <- mean(tmp$rt)
    tmp2 <- subset(tmp, tmp$rt < (mu + (3 * sd)) )
    dat2 <- rbind(dat2, tmp2)
  }
  return(dat2)
}


renameData <- function(dat){
  #----------------------------
  # rename variables
  #----------------------------
  dat$periodicity[dat$periodicity==1] <- "Periodic"
  dat$periodicity[dat$periodicity==0] <- "Aperiodic"
  #return renamed data
  return(dat)
}

add.sdtflags <- function(data){
  #----------------------------
  # Description: add.sdtflags - Add signal detection theory variables to data set and return flagged data
  #----------------------------
  # hit
  data$hit <- ifelse(data$testLevel < 0 & data$answer < 0, 1, 0)
  # miss
  data$miss <- ifelse(data$testLevel < 0 & data$answer > 0, 1, 0)
  # false alarm
  data$fa <- ifelse(data$testLevel > 0 & data$answer > 0, 1, 0)
  # correct rejection
  data$cr<- ifelse(data$testLevel > 0 & data$answer < 0, 1, 0)
  # return
  return(data)
}

get.cleansubset <- function(dat, i){
  #----------------------------
  # Description: get subset of data for participant i (no rts > ±3sd or < 0.250 )
  #----------------------------
  tmp <- subset(dat, dat$subNo == i)
  
  #----------------------------
  # Clean participant data
  #----------------------------
  # remove any RTs that are greater than + 3 sd or less than 250ms.
  rt.sd <- sd(tmp$rt)
  tmp$rtoutlier = 0
  tmp$rtoutlier[tmp$rt > (mean(tmp$rt)+(3*rt.sd))] = 1
  tmp<-tmp[(!tmp$rtoutlier == 1),]
  tmp$rtoutlier <- NULL
  #tmp<-tmp[!(tmp$rt < 0.250),]
  
  #----------------------------
  # Return clean subset for participant i
  #----------------------------
  return(tmp)
}

get.rtdata <- function(data, noPs, pException){
  
  #----------------------------
  # Description: Clean and tag RT data for analysis
  #----------------------------
  return.data <- data.frame()
  #----------------------------
  # Set up the participant look accounting for removed participant numbers
  #----------------------------
  plookup <- c(1:noPs)
  plookup <- plookup[! plookup %in% pException]
  for (i in plookup){
    print(i)
    # make participant subset
    tmp <- subset(dat, dat$subNo == i)
    # normalise testLevel
    val.array <- sort(unique(tmp$testLevel))
    testLevel.NAME <- c("tl-3","tl-2","tl-1","tl0","tl1","tl2","tl3")
    testLevel.NUM <- c(-3,-2,-1, 0, 1, 2, 3)
    for (j in 1:7){
      tmp$testLevelNAME[tmp$testLevel == val.array[j]] <- testLevel.NAME[j]
      tmp$testLevelNUM[tmp$testLevel == val.array[j]] <- testLevel.NUM[j]
    }
    
    # create IOI sd variable
    for (j in 1:length(tmp$index)){
      IOI.var <- c(tmp$t1[j], tmp$t2[j], tmp$t3[j], tmp$t4[j], tmp$t5[j]) 
      tmp$IOIvar[j] <- sd(IOI.var)
    }
    qtiles <- quantile(tmp$IOIvar[tmp$periodicity == "Aperiodic"], probs = seq(0, 1, 0.5))
    for (j in 1:length(tmp$index)){
      if(tmp$periodicity[j] == "Periodic"){
        tmp$IOI.var.lab[j] <- "zero"
      } else if(tmp$IOIvar[j] <= qtiles[2]){
        tmp$IOI.var.lab[j] <- "low"
      } else {
        tmp$IOI.var.lab[j] <- "high"
      }
    }
    
    # remove any trial in which the RT was greater than 3 sd from the rt mean
    m <- mean(tmp$rt)
    sd <- sd(tmp$rt)
    tmp <- tmp[tmp$rt <= (m + (3*sd)),]
    
    # append tmp file onto previous tmp to recreate master dat
    return.data <- rbind(return.data, tmp)
    rm(tmp, m, sd)
  }
  #----------------------------
  # Return clean and tagged rt data for analysis
  #----------------------------
  return(return.data)
}