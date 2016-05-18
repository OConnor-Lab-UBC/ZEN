######################################################
#'
#' @title  EMS analysis for ZEN data
#'
#' @author Mary O'Connor
#' @author oconnor@zoology.ubc.ca
#'
#' @log
#' 
#' 5/17/2016 drafted code
#'
######################################################## 


## load libraries

library(metacom)
library(plyr)
library(dplyr)
library(broom)
library(reshape2)
library(Matrix)
library(lattice)

## load data
subsite <- read.csv("../species data/ZEN_2014_subsitesummary_2016_05_17.csv")
head(subsite)

plot <- read.csv("../species data/ZEN_2014_plot_2016_05_17.csv")
head(plot)

ZENverts <- read.csv("../Emmett files/ZEN_2014_mesograzers_2016-05-17.csv")
head(ZENverts)
names(ZENverts)

### identify sites by their coast and create a merged file
site.geog <- read.csv("../site data/geography.csv")
site.geog$Site <- paste(site.geog$Site.code, site.geog$locale, sep = '.')
lat.long <- read.csv("../site data/ZEN_2014_Lat&Long.csv")
site.info <- merge(site.geog, lat.long, by = 'Site')

plot.info <- plot[,1:5] 
plot.geog <- merge(plot.info, site.info, by = 'Site')

write.csv(site.info, 'siteinfo.csv')
write.csv(plot.geog, 'plotgeog.csv')

ZENverts.sites <- merge(site.info, ZENverts, by = 'Site')


## group sites by ocean basin
sitesP <- ZENverts.sites[(ZENverts.sites$Ocean=='Pacific'),]
sitesA <- ZENverts.sites[(ZENverts.sites$Ocean=='Atlantic'),]
sitesM <- ZENverts.sites[(ZENverts.sites$Ocean=='Med'),]

sitesO <- sitesA

## remove columns, collapse plots by site
sites <- aggregate(sitesO[,-c(1:16, 214:215)], list(sitesO$Site), sum)
dim(sites)
#str(sites)

#Trim the species with 0 occurrences in this dataset
cols.to.delete <- which(colSums(sites[,2:198]) == 0)
sites1 <- sites[!sites$Group.1%in%c('CR.B', 'CR.A'),-(cols.to.delete+1)]  # 
rowSums(sites1[,-1])
rownames(sites1) <- sites1[,1]
sites2 <- sites1[,-1]
sites2[sites2 > 0] <- 1
Metacommunity(sites2, verbose = TRUE, allowEmpty = TRUE) -> meta # could use allowEmpty == TRUE to allow Empty rows and cols in null matrix



## assuming this works, make a plot
a <- as.data.frame(meta[1])

pdf('Pacific.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Pacific",  border="black", scales = list(cex = c(0.5, 0.5), rot = c(90, 90)))
dev.off()

pdf('Atlantic.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Atlantic",  border="black", scales = list(cex = c(0.5, 0.5), rot = c(90, 90)))
dev.off()
