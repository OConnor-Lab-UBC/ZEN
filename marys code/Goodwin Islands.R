### Goodwin Islands EMS 
#########################

library(plyr)
library(reshape2)

data <- read.csv("../species data/GoodwinIslands.csv")
head(data)

#sum to get total abundance across sizes
data$total <- rowSums((data[,c("X8", "X5.6", "X4", "X2.8", "X2", "X1.4", "X1", "X0.71", "X0.5")]))

idata <- data[(data$inshore.offshore == "Inshore"),]
head(idata)
idata1 <- idata[,c(2,4,7,10,23)]

## sum across replicates
idata2 <- ddply(idata1, .(month, year, species.name.revised), summarise, sum(total))
  
idata3 <- melt(idata2, id = 1:2, measure = 3)
