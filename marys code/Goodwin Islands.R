### Goodwin Islands EMS 
### Data from Emmett Duffy's Goodwin Islands eelgrass epiphyte surveys
### Code writte by Mary O'Connor
### May 19 2016
###
### Purpose is to do an Elements of Metacommunity Structure analysis on 
#########################

library(plyr)
library(reshape2)
library(stringr)

data <- read.csv("../species data/GoodwinIslands.csv")
head(data)

#sum to get total abundance across sizes
data$total <- rowSums((data[,c("X8", "X5.6", "X4", "X2.8", "X2", "X1.4", "X1", "X0.71", "X0.5")]))

#remove unnamed taxa
data1 <- data[-which(data$species.name.revised == ""),]

#fix spelling of a taxon
#data[(data$species.name.revised == 'Astyris\xe6lunata')]
#data$species.name.revised[[Astyris\xe6lunata]] <- "Astyrisea lunata"

### maybe split the species name and then paste with .
#data1[,24:26] <- (str_split_fixed(data1$species.name.revised, pattern = ' ', n = 3))
#data1$species <- paste(data1[,(24:26)], sep = '.') this was crashing R; come bact to it later

idata <- data[(data$inshore.offshore == "Inshore"),]
head(idata)
idata1 <- idata[,c(2,4,7,10,23)]

## sum across replicates
idata2 <- ddply(idata1, .(month, year, species.name.revised), summarise, sum(total))
idata2$time.ID <- paste(idata2$month, idata2$year, sep = '.')
 
idata3 <- melt(idata2, id = c(1,2,3,5), measure.vars = "..1")

idata3 <- idata3[,-(1:2)]

idata4 <- dcast(idata3, time.ID ~ species.name.revised, mean)

idata4 <- as.data.frame(idata4)

#remove NaN and Nas. 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
idata4[is.nan.data.frame(idata4)] <- 0
idata4[is.na(idata4)] <- 0

#identify and remove species with no observations, and samples with no species
rowSums(idata4[,-1]) -> idata4$totals
idata5 <- idata4[which(idata4$totals != "0"),] 
#idata5 <- idata5[-is.na(idata5$totals),] 
rownames(idata5) <- idata5[,1]
idata5 <- idata5[, -c(1, ncol(idata5))]
idata5[idata5 > 0] <- 1

## something is weird with the data; start here and remove species for which there are no observations
cols.to.delete <- which(colSums(idata5) == '0')
idata5 <- idata5[, -(cols.to.delete+1)]

# For EMS analysis
library(metacom)
Metacommunity(idata5[,1:40], verbose = TRUE, allowEmpty = TRUE) -> meta 

## assuming this works, make a plot
a <- as.data.frame(meta[1])

pdf('GoodwinIslands.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="GoodwinIslands",  border="black", scales = list(cex = c(0.5, 0.5), x = list(rot = c(90))))
dev.off()

