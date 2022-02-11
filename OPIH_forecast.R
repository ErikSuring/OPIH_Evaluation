#  OPIH_forecast.R
#  Pete Lawson
#  5 January 2012
#  19 January 2012 
#  October 2021 - Erik Suring - Code cleaned, added output figures
#  January 2022 - Erik Suring - Added prediction output tables


require(stats); require(graphics); require(ggplot2)
TimeStamp <- format(Sys.time(), "%a %b %d %X %Y")
results <- "OPIH_results_2022.txt"
if(exists(results)) file.remove(results)

OPIHData <- read.table("PUB2022.txt", header = TRUE, sep = "\t", strip.white = TRUE, comment.char = "#")

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

# Output OPIHData as a .csv file
write.table(OPIHData, file = "OPIHData.csv", sep = ", ", row.names=FALSE)

# Adopted OPIH regression
adopted.OPIH <- lm(AdltAllno83 ~ lagJackOPI + lagSmAdj, data=OPIHData)
pred <- round(predict.lm(adopted.OPIH, OPIHData, interval="prediction"), 1)
output <- cbind(OPIHData$Year, OPIHData$lagJackOPI, OPIHData$lagSmAdj, OPIHData$AdltAll, pred)
colnames(output) <- c("Year", "OPIH_Jacks", "Smolt_Adjustment", "Observed_Abundance", "Predicted_Abundance", "Lower_CI", "Upper_CI")

# Write text file with regression summary and predicted values
sink(file = results, append = FALSE, type = ("output"), split = TRUE)
cat("\n","Run date and time:", TimeStamp,"\n")
print(summary(adopted.OPIH))
print(output)
sink()
write.table(output, file = "OPIH_Forecast.csv", sep = ", ", row.names=FALSE)

#Prediction error graph
prediction.error <- as.data.frame(output)
prediction.error <- subset(prediction.error, Year!=1983)
prediction.error$pred.error <- prediction.error$Observed_Abundance-prediction.error$Predicted_Abundance
prediction.error$pred.proportion <- prediction.error$pred.error/prediction.error$Observed_Abundance*100

coeff <- 0.5
ggplot(data = prediction.error, aes(x=Year, y=pred.error)) +  
  geom_col(aes(y=pred.proportion/coeff), fill = "orange") +
  geom_point(color="black", size=2) +
  ylab("Observed minus forecast (1000s adult abundance)") +
  xlab("Return year") +
  scale_y_continuous(
    breaks = seq(-700,700,100),
    sec.axis = sec_axis(~.*coeff , name = "Error Proportion", , breaks =seq(-200,200,50)) #seq(-220,220,40)
  ) +
  theme_bw()

ggsave("error.png", scale = 1.5, width = 6, height = 3, units = "in")
  
# Observed and predicted abundance graph  
ggplot(data = prediction.error, aes(x=Year, y=Predicted_Abundance)) +  
  geom_col(aes(y=Observed_Abundance), fill = "orange") +
  geom_line(color="grey", size=1.5) +
  geom_point(color="black", size=2) +
  ylab("Coho Salmon abundance (1000s of fish)") +
  xlab("Return year") +
  scale_y_continuous(
    breaks = seq(0,4000,500)) +
  theme_bw()

ggsave("forecast_series.png", scale = 1.5, width = 6, height = 3, units = "in")
