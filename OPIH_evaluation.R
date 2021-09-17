#  Based on OPIH_forecast.R by Pete Lawson, 5 January 2012
#  9/17/21
#  Erik Suring
#  Input file is tab delimited.


require(stats); require(graphics)
OPIHData <- read.table("PUB2021.txt", header = TRUE, sep = "\t", strip.white = TRUE, comment.char = "#")

nYears <- length(OPIHData[,1])
earliestYear <- OPIHData[1,1]
latestYear <- OPIHData[nYears,1]

OPIHData$AdltAll <- c(OPIHData$AdltSRS[1:17],OPIHData$AdltMSM[18:nYears])   # row 18  <-  1986. Use SRS data prior, MSM data post
OPIHData$AdltAllno83 <- OPIHData$AdltAll
OPIHData$AdltAllno83[15] <-  NA   #  1983 is removed from model input data
OPIHData[,"AdltAll"] <- OPIHData$AdltAll
OPIHData[,"AdltAllno83"] <- OPIHData$AdltAllno83

# Lag Jack and Smolt data to return year
OPIHData[,"lagJackCR"] <- c(NA, OPIHData$JackCR[1:nYears-1])
OPIHData[,"lagJackOC"] <- c(NA, OPIHData$JackOC[1:nYears-1])
OPIHData[,"lagJackOPI"] <- OPIHData$lagJackCR + OPIHData$lagJackOC

OPIHData[,"lagSmD"] <- c(NA, OPIHData$SmD[1:nYears-1])
OPIHData[,"lagSmCR"] <- c(NA, OPIHData$SmCR[1:nYears-1])
OPIHData[,"lagSmAdj"] <- OPIHData$lagJackCR*(OPIHData$lagSmD/OPIHData$lagSmCR)   #  Calculate Smolt Adjustment

# Create dataframe for moving window results
if(exists("OPIHModel.MA.results")) rm(OPIHModel.MA.results)

# Run regression for a moving window of 30 years.
for(i in 2:(nYears-30)){  # Starts at row two because first year (1969) is empty as predictor data are lagged
  OPIHData.subset <- OPIHData[i:(i+29),]
  OPIHModel.subset <- lm(AdltAllno83 ~ lagJackOPI + lagSmAdj, data=OPIHData.subset)
  OPIHModel.subset.summary <- summary(OPIHModel.subset)
  
  # Model evaluation results
  RMSE <- sqrt(sum(OPIHModel.subset$residuals^2) / OPIHModel.subset$df)
  Adj.R2 <- OPIHModel.subset.summary$adj.r.squared
  Fstat <- OPIHModel.subset.summary$fstatistic[1]
  JackOPI.tvalue <- OPIHModel.subset.summary$coefficients[2,3]
  SmAdj.tvalue <- OPIHModel.subset.summary$coefficients[3,3]
  Subset.year <- OPIHData.subset[1,1]
  
  # Write results for each subset to a data.frame
  if(exists("OPIHModel.MA.results")){
    YearRow <- data.frame(Subset.year, RMSE, Adj.R2, Fstat, JackOPI.tvalue, SmAdj.tvalue)
    OPIHModel.MA.results <- rbind(OPIHModel.MA.results, YearRow)
  } else {
    OPIHModel.MA.results <- data.frame(Subset.year, RMSE, Adj.R2, Fstat, JackOPI.tvalue, SmAdj.tvalue)
  }

} # for(i in 2:(nYears-30)){

write.table(OPIHModel.MA.results, file = "OPIH_MA_Results.csv", sep = ", ", row.names=FALSE)
