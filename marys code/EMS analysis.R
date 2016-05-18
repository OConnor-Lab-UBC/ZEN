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

## remove columns
sites <- aggregate(ZENverts[,-c(1:3, 201:202)], list(ZENverts$Site), sum)
dim(sites)
str(sites)

#Trim the species with 0 occurrances in this dataset
cols.to.delete <- which(colSums(sites[,2:198]) == 0)
sites1 <- sites[-(sites1$Group.1=='CR.B')&(sites1$Group.1=='CR.A'),-(cols.to.delete+1)]

Metacommunity(sites1) -> meta
meta

