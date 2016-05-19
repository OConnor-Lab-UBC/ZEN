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

sitesO <- ZENverts.sites

## remove columns, collapse plots by site
sites <- aggregate(sitesO[,-c(1:17, 214:216)], list(sitesO$Site), sum)
dim(sites)
#str(sites)

# Trim the species with 0 occurrences in this dataset
cols.to.delete <- which(colSums(sites[,2:197]) == 0)
sites1 <- sites[!sites$Group.1%in%c('CR.B', 'CR.A'),-(cols.to.delete+1)]  # 
rowSums(sites1[,-1])
rownames(sites1) <- sites1[,1]
sites2 <- sites1[,-1]

# for a bray curtis similarity matrix
div <- sites2
div.mat <- vegdist(div)
c

write.csv(as.data.frame(div.mat), 'ZEN2014epifaunadissim.csv')

# For EMS analysis
sites2[sites2 > 0] <- 1
Metacommunity(sites2, verbose = TRUE, allowEmpty = TRUE) -> meta 



## assuming this works, make a plot
a <- as.data.frame(meta[1])

pdf('Pacific.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Pacific",  border="black", scales = list(cex = c(0.5, 0.5), x = list(rot = c(90))))
dev.off()

pdf('Atlantic.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Atlantic",  border="black", scales = list(cex = c(0.5, 0.5), x = list(rot = c(90))))
dev.off()

pdf('Global.pdf', width = 7, height = 9)
levelplot(as.matrix(a), col.regions=c(0,1), region = TRUE, colorkey=FALSE, ylab = '', xlab = '', main="Global",  border="black", scales = list(cex = c(0.5, 0.5), x = list(rot = c(90))))
dev.off()



library(BiodiversityR)
sites <- as.data.frame(unique(levels(as.factor((sitesP$Site)))))
colnames(sites) <- 'site'

BB.A <- rankabundance(sites2, y = sites, factor = "site", level = 'BB.A')
BB.A <- rankabunplot(BB.A)

BB.B <- rankabundance(sites2, y = sites, factor = "site", level = 'BB.B')
BB.B <- rankabunplot(BB.B)

BC.A <- rankabundance(sites2, y = sites, factor = "site", level = 'BC.A')
BC.A <- rankabunplot(BC.A)

BC.B <- rankabundance(sites2, y = sites, factor = "site", level = 'BC.B')
BC.B <- rankabunplot(BC.B)
