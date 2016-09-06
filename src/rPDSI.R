## This function is a wrapper for the PDSI C++ program distributed by the University of Nebraska-Lincoln.
## It takes as input an output directory, monthly temperature and precipitation data, soil available water
## content, and the latitude of the land in question.

## Author: R. Kyle Bocinsky
## Date: 11/26/2012

## The C++ source code can be downloaded from
## http://greenleaf.unl.edu/downloads/scPDSI.zip
## The C++ source must be compiled for the system on which this function will be run, and the compiled program
## must be placed in the same directory as the other R scripts. This script formats the input files appropriately.

## This script only computes the monthly PDSI values, both the original ones, and the self-calibrated values.

rPDSI <- function(output.dir, monthly_T, monthly_P, mon_T_normal, awc, lat, scPDSI.path){
  ## Create the output directory if it does not already exist
  suppressWarnings(dir.create(output.dir))
  
  ## Create simple parameter file
  parameter <- data.frame(awc,lat)
  
  ## Write tables to the output directory
  write.table(as.matrix(monthly_T),file=paste0(output.dir,"monthly_T"),row.names=FALSE, col.names=FALSE, quote=F)
  write.table(as.matrix(monthly_P),file=paste0(output.dir,"monthly_P"),row.names=FALSE, col.names=FALSE, quote=F)
  write.table(as.matrix(mon_T_normal),file=paste0(output.dir,"mon_T_normal"),row.names=FALSE, col.names=FALSE, quote=F)
  write.table(parameter,file=paste0(output.dir,"parameter"),row.names=FALSE, col.names=FALSE)
  
  system(paste0(scPDSI.path," -s -x potentials -i",output.dir," -o",output.dir))
  
  pdsi <- read.table(paste0(output.dir,"monthly/original/PDSI.tbl"))
  pdsi.vect <- c(t(as.matrix(pdsi[,2:13])))
  pdsi.vect[pdsi.vect == -99] <- NA
  
  p_pt <- read.table(paste0(output.dir,"monthly/original/potentials"), header=T)
  p_pt.vect <- c(t(as.matrix(p_pt[,8])))
  p_pt.vect[p_pt.vect == -99] <- NA
  
  pdsi.out <- data.frame(pdsi.vect,p_pt.vect)
  names(pdsi.out) <- c("PDSI","P_PE")

  return(pdsi.out)
}