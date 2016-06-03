###################################################################################
#                                                                                ##
# ZEN 2014 data analysis (JED)                                                   ##
# Data are current as of 18 May 2016                                             ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2016-06-03                                                        ##
#                                                                                ##
###################################################################################

# TO DO:

# Model outputs have changed probably because data set has been updated. Need to 
# begin with latest data set and rerun everything and check model outputs. 

# Make sure we are using the latest data file on OSF
# Do standard model suite for: shoot density, allelic richness
# integrate thes models into SEM?
# 2) Integrate depth/intertidal variable
# 3) Integrate LOI data where available

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA                                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# CREATE MATRICES: MESOGRAZER RICHNESS SPECIES X PLOT AND X SITE                  #
# SUMMARIZE MESOGRAZER RICHNESS BY PLOT AND SITE                                  #
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
#                                                                                 #
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
# LOG TRANSFORMS                                                                  #
# OBTAIN SITE MEANS                                                               #
# EXPLORE DISTRIBUTIONS OF VARIABLES (SITE LEVEL)                                 #
# EXPLORE DATA COMPLETENESS                                                       #
# IMPUTE MISSING DATA                                                             #
#                                                                                 #
# EXPLORATORY PLOTS: PAIRS PANELS                                                 #
# BIVARIATE PLOTS                                                                 #
#                                                                                 #
# IDENTIFY BEST EXPLANATORY VARIABLES USING RANDOM FORESTS - SITE MEANS           #
# PERFORM VIF ANALYSIS: ZOSTERA ABOVE-GROUND MASS (NO OTHER ZOSTERA PREDICTORS)   #
# PERFORM VIF ANALYSIS: PERIPHYTON BIOMASS                                        #
# PERFORM VIF ANALYSIS: CRUSTACEAN ABUNDANCE                                      #
# HIERARCHICAL MODELING APPROACH: RATIONALE                                       #
# BUILD HIERARCHICAL MODELS: ZOSTERA ABOVE-GROUND BIOMASS                         #
#                                                                                 #
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: RANDOM FORESTS                         #
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: PLOT LEVEL (HIERARCHICAL)              #
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: SITE LEVEL (SIMPLE)                    #
#                                                                                 #
# MODELS OF EELGRASS SHOOT DENSITY: PLOT LEVEL (HIERARCHICAL)                     #
# MODELS OF EELGRASS SHOOT DENSITY: SITE LEVEL (SIMPLE)                           #
#                                                                                 #
# MODELS OF EELGRASS PRODUCTIVITY: RANDOM FORESTS                                 #
# MODELS OF EELGRASS PRODUCTIVITY: PLOT LEVEL (HIERARCHICAL)                      #
# MODELS OF EELGRASS PRODUCTIVITY: SITE LEVEL (SIMPLE)                            #
#                                                                                 #
# MODELS OF PERIPHYTON BIOMASS: RANDOM FORESTS                                    #
# MODELS OF PERIPHYTON BIOMASS: PLOT LEVEL (HIERARCHICAL)                         #
# MODELS OF PERIPHYTON BIOMASS: SITE LEVEL (SIMPLE)                               #
#                                                                                 #
# MODELS OF CRUSTACEAN BIOMASS: RANDOM FORESTS                                    #
# MODELS OF CRUSTACEAN BIOMASS: PLOT LEVEL (HIERARCHICAL)                         #
# MODELS OF CRUSTACEAN BIOMASS: SITE LEVEL (SIMPLE)                               #
#                                                                                 #
# MODELS OF GAMMARID BIOMASS: RANDOM FORESTS                                      #
# MODELS OF GAMMARID BIOMASS: PLOT LEVEL (HIERARCHICAL)                           #
# MODELS OF GAMMARID BIOMASS: SITE LEVEL (SIMPLE)                                 #
#                                                                                 #
# MODELS OF GASTROPOD BIOMASS: RANDOM FORESTS                                     #
# MODELS OF GASTROPOD BIOMASS: PLOT LEVEL (HIERARCHICAL)                          #
# MODELS OF GASTROPOD BIOMASS: SITE LEVEL (SIMPLE)                                #
#                                                                                 #
# SEM: CRUSTACEAN BIOMASS                                                         #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA                                                                        #
###################################################################################

# Source data: ZEN_2014_Site&PlotData_2016_05_17_Released.xlsx                    
                                                                                
# This script analyzes data from the ZEN 2014 field sampling program in 2014. 
# The script started at ZEN workshop, Dvais, CA, May 2016.

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(plyr)         
library(reshape)
library(reshape2)
library(plotrix)
library(ggplot2)
library(psych)
library(stringr)
library(nlme)
library(piecewiseSEM)
library(randomForest)
library(car)
library(MuMIn)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################


# Read in summary SITE data: 
SubsiteSummary <- read.csv(file = "ZEN_2014_Site&PlotData_2016_05_17_Released_subsite_summary.csv", header = TRUE)
names(SubsiteSummary)
str(SubsiteSummary)
# NOTE: Variables of interest from this file are primarily site-level metadata, plus 
# genetic data (which are site-level metrics). The many site-level means are recalucklated below.  

# Add variable to site-level data set designating whether the site was glaciated (1) or not (0) 
# during Pleistocene. Note: this is a quick and dirty estimate from a low-res map of max glacial extent. 
levels(SubsiteSummary$Site.Code)
SubsiteSummary$glaciation <- c(0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0)

# Read in site-level data on fetch and human population density compiled by Johan Eklof, 2016-05-18:
fetch <- read.csv(file = "ZEN 2014 Fetch estimates per subsite.csv", header = TRUE)
# Estimation of population density based on SEDAC predictions of pop density in 2015: 
# http://sedac.ciesin.columbia.edu/data/set/gpw-v3-population-density-future-estimates
popn <- read.csv(file = "ZEN 2014 Pop dens estimates per subsite.csv", header = TRUE)

# Integrate fetch and population density data into subsite summary dataframe
SubsiteSummary$mean.fetch <- fetch$Mean.Fetch[match(SubsiteSummary$Site, fetch$Site)]
SubsiteSummary$pop.density.2015 <- popn$PopDens2[match(SubsiteSummary$Site, popn$Site)]


# Read in plot data: 
ZEN_2014_plot_data <- read.csv(file = "ZEN_2014_Site&PlotData_2016_05_17_Released.csv", header = TRUE)
names(ZEN_2014_plot_data)

# Remove redundant variables
ZEN_2014_plot_data <- ZEN_2014_plot_data[,!colnames(ZEN_2014_plot_data) %in% c("X")]
str(ZEN_2014_plot_data)

# Rename misspelled or confusing variables
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Sheath.Width.cm."] <- "Zostera.sheath.width"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Shealth.Length.cm."] <- "Zostera.sheath.length"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Longest.Leaft.Length.cm."] <- "Zostera.longest.leaf.length"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Amphipod"] <- "amphipod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Caprellid"] <- "caprellid.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Gastropod"] <- "gastropod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Isopod"] <- "isopod.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Lettuce"] <- "lettuce.survival.24hr"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Pdn.Squid"] <- "squid.survival.24hr"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Above.Zmarina.g"] <- "Zostera.aboveground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Above.Total.macrophytes.core.g"] <- "total.macrophytes.aboveground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Mean.Below.Zmarina.g"] <- "Zostera.belowground.mean.mass"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Shoots.Zmarina.per.m2"] <- "Zostera.shoots.per.m2.core"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Epibiota.Periphyton"] <- "periphyton.mass.per.g.zostera"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Mesograzers"] <- "mesograzer.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Molluscs"] <- "mollusc.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Crustaceans"] <- "crustacean.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Abund.Amphipods"] <- "amphipod.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Gammarids"] <- "gammarid.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Caprellids"] <- "caprellid.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Gastropods"] <- "gastropod.abund.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Abund.Isopods"] <- "isopod.abund.per.g.plant"

names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Mesograzers"] <- "mesograzer.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Molluscs"] <- "mollusc.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Crustaceans"] <- "crustacean.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Total.Biomass.Amphipods"] <- "amphipod.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Gammarids"] <- "gammarid.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Caprellids"] <- "caprellid.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Gastropods"] <- "gastropod.mass.per.g.plant"
names(ZEN_2014_plot_data)[names(ZEN_2014_plot_data)=="Std.Biomass.Isopods"] <- "isopod.mass.per.g.plant"

names(ZEN_2014_plot_data)

# Read in plot-level (mesograzer) data to calculate site-level grazer richness.
epifauna <- read.csv(file = "ZENnverts.v8.csv", header = TRUE)

names(epifauna)

epifauna_clean <- epifauna[,!colnames(epifauna) %in% c("Data.Entered.by", "Measured.size",
  "Measured.Size","Mesured.Size","Mesured.Size.1","Mesured.Size.2","Mesured.Size.3",
  "Mesured.Size.4","Mesured.Size.5","Mesured.Size.6","Mesured.Size.7","Mesured.Size.8",
  "individuals.SizeCount","DW..g.","Vegetative.or.Flower","Shoot.ID", "Total.Abundance.to.1.mm",
  "Notes")]

names(epifauna_clean)

# Remove Sampling.Time = 2
epifauna_clean <- droplevels(subset(epifauna_clean, Sampling.Time == "1"))


###################################################################################
# CREATE MATRICES: MESOGRAZER RICHNESS SPECIES X PLOT AND X SITE                  #
###################################################################################

# PLOT level:

# Create a dataframe with biomass of MESOGRAZER SPECIES as columns
# Melt the dataframe so that each PLOT is a row:

names(epifauna_clean)

# Remove any non-grazers
grazers <- droplevels(subset(epifauna_clean, Type == "Mesograzer"))
levels(grazers$Type) # good

mesograzers_plot <- melt(grazers, id.var = c("Site", "Plot.ID", "unique.ID","Species"), measure.var = c("Total.Abundance"))
# Cast the molten data frame so that each SPECIES is a column across the top:
mesograzers_plot <- dcast(mesograzers_plot, Site + Plot.ID + unique.ID ~ Species, sum)
str(mesograzers_plot)
names(mesograzers_plot)

# Export data frame as csv file: grazer species x plot ID matrix:
write.csv(mesograzers_plot, "ZEN_2014_mesograzers_PLOT-LEVEL_2016-05-17.csv", row.names = F)


# SITE level:

# Create a dataframe with biomass of MESOGRAZER SPECIES as columns
# Melt the dataframe so that each SITE is a row:

mesograzers_site <- melt(grazers, id.var = c("Site","Species"), measure.var = c("Total.Abundance"))
# Cast the molten data frame so that each SPECIES is a column across the top:
mesograzers_site <- dcast(mesograzers_site, Site ~ Species, sum)
str(mesograzers_site)
names(mesograzers_site)

# Export data frame as csv file: grazer species x plot ID matrix:
write.csv(mesograzers_site, "ZEN_2014_mesograzers_SITE-LEVEL_2016-05-17.csv", row.names = F)


###################################################################################
# SUMMARIZE MESOGRAZER RICHNESS BY PLOT AND SITE                                  #
###################################################################################

# Create variable: PLOT-level grazer richness
mesograzers_plot$richness.plot = 
  rowSums(mesograzers_plot[, !colnames(mesograzers_plot) %in% c("Site", "Plot.ID", "unique.ID")] > 0, na.rm = T)

# Create variable: SITE-level grazer richness 
richness.site <- ddply(mesograzers_plot, c("Site"), function(x) {
  # Sum abundances 
  sum.x = colSums(x[, !colnames(x) %in% c("Site", "Plot.ID", "unique.ID")], na.rm = T)
  # Get richness
  sum(sum.x > 0, na.rm = T)
})

# change name of newly created variable 
colnames(richness.site)[2] <- c("richness.site")

# Create a dataframe with the richness values by site
temp <- mesograzers_plot[c("Site", "Plot.ID", "unique.ID", "richness.plot")]
grazer.richness <- merge(temp, richness.site, by = c("Site"))
names(grazer.richness)


###################################################################################
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
###################################################################################

names(grazer.richness)
colnames(grazer.richness)[4:5] <- c("grazer.richness.plot", "grazer.richness.site") 

# Bind the data frames into a master data frame
ZEN_2014_master_data = ZEN_2014_plot_data
ZEN_2014_master_data$grazer.richness.plot <- grazer.richness$grazer.richness.plot[match(ZEN_2014_master_data$Site, grazer.richness$Site)]
ZEN_2014_master_data$grazer.richness.site <- grazer.richness$grazer.richness.site[match(ZEN_2014_master_data$Site, grazer.richness$Site)]

names(ZEN_2014_master_data)

# Select those environmental variables we want to include:
names(SubsiteSummary)

envdata <- SubsiteSummary[c("Site", "Subsite", "Subsite.Name", "Ocean", "Coast", "Basin", "Latitude", 
                            "Longitude", "Date.Collected", "Temperature.C", "Salinity.ppt",
                            "Day.length.hours", "Depth.Categorical", "Depth.m", "Perc.Silt", "Perc.Sand", "Perc.Gravel", 
                            "glaciation", "mean.fetch", "pop.density.2015" ,"GenotypicRichness", "AllelicRichness", "Inbreeding")]

ZEN_2014_master_data <- merge(ZEN_2014_master_data, envdata, by = c("Site"))
names(ZEN_2014_master_data)
nrow(ZEN_2014_master_data) # 1000



###################################################################################
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
###################################################################################

# Examine frequency distribution of sites by environmental factor
# par(mfrow = c(1,1))
par(mfrow = c(2,4))
hist(ZEN_2014_master_data$Latitude, col = "cyan", main = "Surveys by latitude")    
hist(ZEN_2014_master_data$Longitude, col = "cyan", main = "Surveys by longitude")    
hist(ZEN_2014_master_data$Temperature.C, col = "cyan", main = "Surveys by temperature")    
hist(ZEN_2014_master_data$Salinity.ppt, col = "cyan", main = "Surveys by salinity")    

hist(ZEN_2014_master_data$Depth.m, col = "cyan", main = "Surveys by depth")    
hist(ZEN_2014_master_data$Day.length.hours, col = "cyan", main = "Surveys by day length")    
hist(ZEN_2014_master_data$mean.fetch, col = "cyan", main = "Surveys by mean fetch")    
hist(ZEN_2014_master_data$pop.density.2015, col = "cyan", main = "Surveys by population density")    

hist(ZEN_2014_master_data$Zostera.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$total.macrophytes.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$periphyton.mass.per.g.zostera, col = "cyan")    
hist(ZEN_2014_master_data$crustacean.abund.per.g.plant, col = "cyan")    

hist(ZEN_2014_master_data$gastropod.abund.per.g.plant, col = "cyan")    
hist(ZEN_2014_master_data$grazer.richness.plot, col = "cyan")    
hist(ZEN_2014_master_data$amphipod.survival, col = "cyan")    

hist(ZEN_2014_master_data$mean.fetch, col = "cyan")    
hist(ZEN_2014_master_data$pop.density.2015, col = "cyan")    


###################################################################################
# LOG TRANSFORMS                                                                  #
###################################################################################

names(ZEN_2014_master_data)

# NOTE: For many variables I add a constant roughly equal to the smallest value recorded 

ZEN_2014_master_data$log10.mean.fetch <- log10(ZEN_2014_master_data$mean.fetch) 

ZEN_2014_master_data$log10.Zostera.AG.mass <- log10(ZEN_2014_master_data$Zostera.aboveground.mean.mass + 1) 
ZEN_2014_master_data$log10.macrophytes.total.AG.mass.core <- log10(ZEN_2014_master_data$total.macrophytes.aboveground.mean.mass + 1) 
ZEN_2014_master_data$log10.Zostera.BG.mass <- log10(ZEN_2014_master_data$Zostera.belowground.mean.mass + 1) 
ZEN_2014_master_data$log10.Zostera.shoots.core <- log10(ZEN_2014_master_data$Zostera.shoots.per.m2.core + 1) 
ZEN_2014_master_data$log10.Zostera.sheath.width <- log10(ZEN_2014_master_data$Zostera.sheath.width) 
ZEN_2014_master_data$log10.Zostera.sheath.length <- log10(ZEN_2014_master_data$Zostera.sheath.length) 
ZEN_2014_master_data$log10.Zostera.longest.leaf.length <- log10(ZEN_2014_master_data$Zostera.longest.leaf.length) 
ZEN_2014_master_data$log10.periphyton.mass.per.g.zostera <- log10(ZEN_2014_master_data$periphyton.mass.per.g.zostera + 0.0001) 

ZEN_2014_master_data$log10.mesograzer.abund.per.g.plant <- log10(ZEN_2014_master_data$mesograzer.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.crustacean.abund.per.g.plant <- log10(ZEN_2014_master_data$crustacean.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.amphipod.abund.per.g.plant <- log10(ZEN_2014_master_data$amphipod.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gammarid.abund.per.g.plant <- log10(ZEN_2014_master_data$gammarid.abund.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gastropod.abund.per.g.plant <- log10(ZEN_2014_master_data$gastropod.abund.per.g.plant + 0.1) 

ZEN_2014_master_data$log10.mesograzer.mass.per.g.plant <- log10(ZEN_2014_master_data$mesograzer.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.crustacean.mass.per.g.plant <- log10(ZEN_2014_master_data$crustacean.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.amphipod.mass.per.g.plant <- log10(ZEN_2014_master_data$amphipod.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gammarid.mass.per.g.plant <- log10(ZEN_2014_master_data$gammarid.mass.per.g.plant + 0.1) 
ZEN_2014_master_data$log10.gastropod.mass.per.g.plant <- log10(ZEN_2014_master_data$gastropod.mass.per.g.plant + 0.1) 

ZEN_2014_master_data$log10.grazer.richness.plot <- log10(ZEN_2014_master_data$grazer.richness.plot + 1) 
ZEN_2014_master_data$log10.grazer.richness.site <- log10(ZEN_2014_master_data$grazer.richness.site + 1) 

ZEN_2014_master_data$log10.Leaf.PercN <- log10(ZEN_2014_master_data$Leaf.PercN) 
ZEN_2014_master_data$log10.GenotypicRichness <- log10(ZEN_2014_master_data$GenotypicRichness) 
ZEN_2014_master_data$log10.AllelicRichness <- log10(ZEN_2014_master_data$AllelicRichness) 

# Obtain values for mean body mass per individual
ZEN_2014_master_data$crustacean.mean.body.mass <- ZEN_2014_master_data$crustacean.mass.per.g.plant / ZEN_2014_master_data$crustacean.abund.per.g.plant 
ZEN_2014_master_data$log10.crustacean.mean.body.mass <- log10(ZEN_2014_master_data$crustacean.mean.body.mass)

ZEN_2014_master_data$gastropod.mean.body.mass <- ZEN_2014_master_data$gastropod.mass.per.g.plant / ZEN_2014_master_data$gastropod.abund.per.g.plant 
ZEN_2014_master_data$log10.gastropod.mean.body.mass <- log10(ZEN_2014_master_data$gastropod.mean.body.mass)

ZEN_2014_master_data$mesograzer.mean.body.mass <- ZEN_2014_master_data$mesograzer.mass.per.g.plant / ZEN_2014_master_data$mesograzer.abund.per.g.plant 
ZEN_2014_master_data$log10.mesograzer.mean.body.mass <- log10(ZEN_2014_master_data$mesograzer.mean.body.mass)

ZEN_2014_master_data$basin <- ZEN_2014_master_data$Coast
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='East.Atlantic'] <- 'East'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='East.Pacific'] <- 'East'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='West.Atlantic'] <- 'West'
levels(ZEN_2014_master_data$basin)[levels(ZEN_2014_master_data$basin)=='West.Pacific'] <- 'West'


# Change values of NaN to NA:
ZEN_2014_master_data[ZEN_2014_master_data == "NaN"] = NA 
ZEN_2014_master_data[ZEN_2014_master_data == "Inf"] = NA 
names(ZEN_2014_master_data)

# Export the data frame
write.csv(ZEN_2014_master_data, "ZEN_2014_master_data_2016-05-19.csv", row.names = F)


###################################################################################
# OBTAIN SITE MEANS                                                               #
###################################################################################

# Obtain mean values per site
ZEN_2014_site_means <- ddply(ZEN_2014_master_data, c("Site"), summarize, 
                             Zostera.AG.mass.site = mean(Zostera.aboveground.mean.mass, na.rm = T), 
                             # macroalgae.total.AG.mass.core.site = mean(macroalgae.aboveground.mean.mass, na.rm = T),                      
                             # seagrass.othersp.AG.mass.site = mean(seagrass.aboveground.mean.mass, na.rm = T),                      
                             macrophytes.total.AG.mass.core.site = mean(total.macrophytes.aboveground.mean.mass, na.rm = T),                      
                             Zostera.BG.mass.site = mean(Zostera.belowground.mean.mass, na.rm = T),                      
                             # macrophytes.total.BG.mass.core.site = mean(total.macrophytes.belowground.mean.mass, na.rm = T),                      
                             Zostera.shoots.core.site = mean(Zostera.shoots.per.m2.core, na.rm = T),                      
                             Zostera.sheath.width.site = mean(Zostera.sheath.width, na.rm = T),                      
                             Zostera.sheath.length.site = mean(Zostera.sheath.length, na.rm = T),                      
                             Zostera.longest.leaf.length.site = mean(Zostera.longest.leaf.length, na.rm = T),                      
                             # pct.cover.macroalgae.site = mean(pct.cover.macroalgae, na.rm = T),                      
                             # pct.cover.seagrass.site = mean(pct.cover.seagrass, na.rm = T),                      
                             periphyton.mass.per.g.zostera.site = mean(periphyton.mass.per.g.zostera, na.rm = T),                      
                             
                             mesograzer.abund.per.g.plant.site = mean(mesograzer.abund.per.g.plant, na.rm = T),                      
                             crustacean.abund.per.g.plant.site = mean(crustacean.abund.per.g.plant, na.rm = T), 
                             amphipod.abund.per.g.plant.site = mean(amphipod.abund.per.g.plant, na.rm = T), 
                             gammarid.abund.per.g.plant.site = mean(gammarid.abund.per.g.plant, na.rm = T), 
                             caprellid.abund.per.g.plant.site = mean(caprellid.abund.per.g.plant, na.rm = T), 
                             isopod.abund.per.g.plant.site = mean(isopod.abund.per.g.plant, na.rm = T), 
                             gastropod.abund.per.g.plant.site = mean(gastropod.abund.per.g.plant, na.rm = T), 
                             # polychaete.abund.per.g.plant.site = mean(polychaete.abund.per.g.plant, na.rm = T), 
                             
                             mesograzer.mass.per.g.plant.site = mean(mesograzer.mass.per.g.plant, na.rm = T),                      
                             crustacean.mass.per.g.plant.site = mean(crustacean.mass.per.g.plant, na.rm = T), 
                             amphipod.mass.per.g.plant.site = mean(amphipod.mass.per.g.plant, na.rm = T), 
                             gammarid.mass.per.g.plant.site = mean(gammarid.mass.per.g.plant, na.rm = T), 
                             caprellid.mass.per.g.plant.site = mean(caprellid.mass.per.g.plant, na.rm = T), 
                             isopod.mass.per.g.plant.site = mean(isopod.mass.per.g.plant, na.rm = T), 
                             gastropod.mass.per.g.plant.site = mean(gastropod.mass.per.g.plant, na.rm = T), 
                             # polychaete.mass.per.g.plant.site = mean(polychaete.mass.per.g.plant, na.rm = T), 
                             
                             crustacean.mean.body.mass.site = mean(crustacean.mean.body.mass, na.rm = T), 
                             log10.crustacean.mean.body.mass.site = mean(log10.crustacean.mean.body.mass, na.rm = T), 
                             gastropod.mean.body.mass.site = mean(gastropod.mean.body.mass, na.rm = T), 
                             log10.gastropod.mean.body.mass.site = mean(log10.gastropod.mean.body.mass, na.rm = T), 
                             mesograzer.mean.body.mass.site = mean(mesograzer.mean.body.mass, na.rm = T), 
                             log10.mesograzer.mean.body.mass.site = mean(log10.mesograzer.mean.body.mass, na.rm = T), 
                             
                             grazer.richness.plot.site = mean(grazer.richness.plot, na.rm = T), 
                             amphipod.survival.24hr.site = mean(amphipod.survival.24hr, na.rm = T), 
                             gastropod.survival.24hr.site = mean(gastropod.survival.24hr, na.rm = T), 
                             squid.survival.24hr.site = mean(squid.survival.24hr, na.rm = T), 
                             Leaf.PercN.site = mean(Leaf.PercN, na.rm = T), 
                             
                             log10.Zostera.AG.mass.site = mean(log10.Zostera.AG.mass, na.rm = T), 
                             log10.macrophytes.total.AG.mass.core.site = mean(log10.macrophytes.total.AG.mass.core, na.rm = T), 
                             log10.Zostera.BG.mass.site = mean(log10.Zostera.BG.mass, na.rm = T), 
                             log10.Zostera.shoots.core.site = mean(log10.Zostera.shoots.core, na.rm = T), 
                             log10.Zostera.sheath.width.site = mean(log10.Zostera.sheath.width, na.rm = T), 
                             log10.Zostera.sheath.length.site = mean(log10.Zostera.sheath.length, na.rm = T), 
                             log10.Zostera.longest.leaf.length.cm.site = mean(log10.Zostera.longest.leaf.length, na.rm = T), 
                             # log10.pct.cover.macroalgae.site = mean(log10.pct.cover.macroalgae, na.rm = T), 
                             # log10.pct.cover.seagrass.site = mean(log10.pct.cover.seagrass, na.rm = T), 
                             log10.periphyton.mass.per.g.zostera.site = mean(log10.periphyton.mass.per.g.zostera, na.rm = T), 
                             
                             log10.mesograzer.abund.per.g.plant.site = mean(log10.mesograzer.abund.per.g.plant, na.rm = T), 
                             log10.crustacean.abund.per.g.plant.site = mean(log10.crustacean.abund.per.g.plant, na.rm = T), 
                             log10.amphipod.abund.per.g.plant.site = mean(log10.amphipod.abund.per.g.plant, na.rm = T), 
                             log10.gammarid.abund.per.g.plant.site = mean(log10.gammarid.abund.per.g.plant, na.rm = T), 
                             log10.gastropod.abund.per.g.plant.site = mean(log10.gastropod.abund.per.g.plant, na.rm = T), 
                             
                             log10.mesograzer.mass.per.g.plant.site = mean(log10.mesograzer.mass.per.g.plant, na.rm = T), 
                             log10.crustacean.mass.per.g.plant.site = mean(log10.crustacean.mass.per.g.plant, na.rm = T), 
                             log10.amphipod.mass.per.g.plant.site = mean(log10.amphipod.mass.per.g.plant, na.rm = T), 
                             log10.gammarid.mass.per.g.plant.site = mean(log10.gammarid.mass.per.g.plant, na.rm = T), 
                             log10.gastropod.mass.per.g.plant.site = mean(log10.gastropod.mass.per.g.plant, na.rm = T), 
                             
                             log10.grazer.richness.plot.site = mean(log10.grazer.richness.plot, na.rm = T), 
                             log10.Leaf.PercN.site = mean(log10.Leaf.PercN, na.rm = T) 
)

# Change values of NaN to NA:
ZEN_2014_site_means[ZEN_2014_site_means == "NaN"] = NA 

# Add site-level environmental (and other) variables back in
ZEN_2014_site_means <- merge(envdata, ZEN_2014_site_means, by = c("Site"))
ZEN_2014_site_means$log10.mean.fetch <- ZEN_2014_master_data$log10.mean.fetch[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$pop.density.2015 <- ZEN_2014_master_data$pop.density.2015[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.GenotypicRichness <- ZEN_2014_master_data$log10.GenotypicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.AllelicRichness <- ZEN_2014_master_data$log10.AllelicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$grazer.richness.site <- ZEN_2014_master_data$grazer.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.grazer.richness.site <- ZEN_2014_master_data$log10.grazer.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]

names(ZEN_2014_site_means)

# Export the data frame
write.csv(ZEN_2014_site_means, "ZEN_2014_site_means_2016-05-19.csv", row.names = F)

# Create separate data sets for Atlantic and Pacific sites (subsite summary)
ZEN_2014_site_means.Atl <- droplevels(subset(ZEN_2014_site_means, Ocean == "Atlantic"))
ZEN_2014_site_means.Pac <- droplevels(subset(ZEN_2014_site_means, Ocean == "Pacific"))

# Create separate data sets for Atlantic and Pacific sites (subsite summary)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))


###################################################################################
# EXPLORE DISTRIBUTIONS OF VARIABLES (SITE LEVEL)                                 #
###################################################################################

names(ZEN_2014_site_means)

# Examine frequency distribution of sites by environmental factor
# par(mfrow = c(1,1))
par(mfrow = c(2,4))
hist(ZEN_2014_site_means$Latitude, col = "cyan", main = "Surveys by latitude")    
hist(ZEN_2014_site_means$Longitude, col = "cyan", main = "Surveys by longitude")    
hist(ZEN_2014_site_means$Temperature.C, col = "cyan", main = "Surveys by temperature")    
hist(ZEN_2014_site_means$Salinity.ppt, col = "cyan", main = "Surveys by salinity")    

# hist(ZEN_2014_site_means$Depth.m, col = "cyan")    

hist(ZEN_2014_site_means$log10.Zostera.AG.mass.site, col = "cyan")    
hist(ZEN_2014_site_means$log10.macrophytes.total.AG.mass.core.site, col = "cyan")    
hist(ZEN_2014_site_means$log10.Zostera.shoots.core.site, col = "cyan")    
hist(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site, col = "cyan")    

hist(ZEN_2014_site_means$log10.crustacean.mass.per.g.plant.site, col = "cyan")    
hist(ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site, col = "cyan")    
hist(ZEN_2014_site_means$log10.grazer.richness.plot.site, col = "cyan")    
hist(ZEN_2014_site_means$amphipod.survival.24hr.site, col = "cyan")    


###################################################################################
# EXPLORE DATA COMPLETENESS                                                       #
###################################################################################

# NOTE: AIC comparisons among models are invalid unless exactly the same number of plots
# are used in each comparison, because the DF influences calculation of the AIC score. 
# This means that we need data on all plots and might need to impute missing data for 
# valid AIC model comparisons. 

# How many observations are missing for each variable?
sum(is.na(ZEN_2014_master_data$log10.Zostera.AG.mass)) # 24
sum(is.na(ZEN_2014_master_data$log10.Zostera.shoots.core)) # 15
sum(is.na(ZEN_2014_master_data$Zostera.longest.leaf.length)) # 0
# sum(is.na(ZEN_2014_master_data$pct.cover.macroalgae)) # 20
# sum(is.na(ZEN_2014_master_data$pct.cover.seagrass)) # 20
sum(is.na(ZEN_2014_master_data$grazer.richness.site)) # 0
sum(is.na(ZEN_2014_master_data$Leaf.PercN)) # 14
sum(is.na(ZEN_2014_master_data$Temperature.C)) # 0
sum(is.na(ZEN_2014_master_data$Salinity.ppt)) # 0
sum(is.na(ZEN_2014_master_data$Day.length.hours)) # 0
sum(is.na(ZEN_2014_master_data$Perc.Sand)) # 420
sum(is.na(ZEN_2014_master_data$GenotypicRichness)) # 0
sum(is.na(ZEN_2014_master_data$AllelicRichness)) # 0 
sum(is.na(ZEN_2014_master_data$Inbreeding)) # 0
sum(is.na(ZEN_2014_master_data$log10.periphyton.mass.per.g.zostera)) # 36
sum(is.na(ZEN_2014_master_data$log10.mesograzer.abund.per.g.plant)) # 22
sum(is.na(ZEN_2014_master_data$log10.crustacean.abund.per.g.plant)) # 22
sum(is.na(ZEN_2014_master_data$log10.gastropod.abund.per.g.plant)) # 22
sum(is.na(ZEN_2014_master_data$amphipod.survival.24hr)) # 191

# How many observations are available for each variable?
sum(!is.na(ZEN_2014_master_data$log10.Zostera.AG.mass)) # 976
sum(!is.na(ZEN_2014_master_data$log10.Zostera.shoots.core)) # 985
sum(!is.na(ZEN_2014_master_data$Zostera.longest.leaf.length)) # 1000
# sum(!is.na(ZEN_2014_master_data$pct.cover.macroalgae)) # 968
# sum(!is.na(ZEN_2014_master_data$pct.cover.seagrass)) # 968
sum(!is.na(ZEN_2014_master_data$grazer.richness.site)) # 1000
sum(!is.na(ZEN_2014_master_data$amphipod.survival.24hr)) # 809
sum(!is.na(ZEN_2014_master_data$squid.survival.24hr)) # 976
sum(!is.na(ZEN_2014_master_data$Leaf.PercN)) # 986
sum(!is.na(ZEN_2014_master_data$Temperature.C)) # 1000
sum(!is.na(ZEN_2014_master_data$Salinity.ppt)) # 1000
sum(!is.na(ZEN_2014_master_data$Day.length.hours)) # 1000
sum(!is.na(ZEN_2014_master_data$Perc.Sand)) # 580
sum(!is.na(ZEN_2014_master_data$GenotypicRichness)) # 1000
sum(!is.na(ZEN_2014_master_data$AllelicRichness)) # 1000
sum(!is.na(ZEN_2014_master_data$Inbreeding)) # 1000
sum(!is.na(ZEN_2014_master_data$log10.periphyton.mass.per.g.zostera)) # 964
sum(!is.na(ZEN_2014_master_data$log10.mesograzer.abund.per.g.plant)) # 978
sum(!is.na(ZEN_2014_master_data$log10.crustacean.abund.per.g.plant)) # 978
sum(!is.na(ZEN_2014_master_data$log10.gastropod.abund.per.g.plant)) # 978
sum(!is.na(ZEN_2014_master_data$amphipod.survival.24hr)) # 809

# NOTE: One big reason analyses below are dropping so many plots is because nearly
# half are missing data for sediment grain size. So let's omit grain size until this is fixed.

# Look at percentage of values missing for each variable
# First create function to calculate % of missing values infor each variable in a data frame… 
pMiss <- function(x){sum(is.na(x))/length(x)*100}
# Now apply it to the data frame: 
apply(ZEN_2014_master_data,2,pMiss)

# Some relevant results:
#   Most variables have fewer < 4% missing -  good. 
#   Exception is predation data since many sites did not have one or another tupe of grazer
#   sediment grain size (e.g., Perc.Sand = 42.510121)


###################################################################################
# IMPUTE MISSING DATA                                                             #
###################################################################################

# First the rationale: There is no getting around imputation of missing values if we
# want to do model comparisons. Which we do. This is because models compared with AIC 
# MUST be based on the same data set and have same DF. If values are missing such that one 
# model has different DF than another with which it is compared, the AIC scores used 
# to compare them will be invalid. The alternative to imputation would be to remove all 
# rows (plots) that do not contain a value for EVERY variable in the data set. This is
# far the worse alternative since it will end up throwing away perhaps 20% of the entire 
# data set. 

# First, realize that we cannot include the following predictors because biased by 
# missing from entire site: periphyton (missing from SW.A). Also BC.A is missing 12 
# of 20 samples for Zostera biomass and shoot density. I am going ahead to impute the 
# missing Zostera data for BC.A but we may want to exclude this site from any analysis 
# that uses Zostera mass or shoot density as either a predictor or  response. 

# Note that data set for RF imputation model cannot include any rows that lack data for one of 
# the predictor variables.

# We first impute missing values for Zostera variables (AG biomass, shoot density, %N), then
# proceed to the grazer variables that depend on these grass variables

y <- ZEN_2014_master_data


# LEAF % NITROGEN
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.Leaf.PercN)) # 14
leafN.rf = randomForest(log10.Leaf.PercN ~ Ocean + Coast + Latitude + Longitude + 
                          Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness 
                        # + log10.Zostera.AG.mass + log10.Zostera.shoots.core 
                        + log10.Zostera.sheath.width + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length 
                        + log10.grazer.richness.site + log10.AllelicRichness,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Leaf.PercN), 
  "log10.Leaf.PercN"] = predict(leafN.rf, y[is.na(y$log10.Leaf.PercN), ])
sum(is.na(y$log10.Leaf.PercN)) # 0


# ZOSTERA SHOOT DENSITY
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.Zostera.shoots.core)) # 15
shootdensity.rf = randomForest(log10.Zostera.shoots.core ~ Ocean + Coast + Latitude + Longitude + 
                                 Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness 
                               # + log10.Zostera.AG.mass 
                               + log10.Leaf.PercN 
                               + log10.Zostera.sheath.width + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length 
                               + log10.grazer.richness.site + log10.AllelicRichness,
                               na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.shoots.core), 
  "log10.Zostera.shoots.core"] = predict(shootdensity.rf, y[is.na(y$log10.Zostera.shoots.core), ])
sum(is.na(y$log10.Zostera.shoots.core)) # 0


# ZOSTERA ABOVE-GROUND BIOMASS
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.Zostera.AG.mass)) # 24
ZAG.rf = randomForest(log10.Zostera.AG.mass ~ Ocean + Coast + Latitude + Longitude + 
                        Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness 
                      + log10.Zostera.shoots.core + log10.Leaf.PercN 
                      + log10.Zostera.sheath.width + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length 
                      + log10.grazer.richness.site + log10.AllelicRichness,
                      na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.AG.mass), 
  "log10.Zostera.AG.mass"] = predict(ZAG.rf, y[is.na(y$log10.Zostera.AG.mass), ])
sum(is.na(y$log10.Zostera.AG.mass)) # 0


# CRUSTACEAN mesograzers
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.crustacean.mass.per.g.plant)) # 23
crust.rf = randomForest(log10.crustacean.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude + 
                          Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness + log10.Zostera.AG.mass 
                        + log10.Zostera.shoots.core + log10.Zostera.sheath.width + log10.Zostera.sheath.length + 
                          log10.Zostera.longest.leaf.length + log10.grazer.richness.site + log10.AllelicRichness,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.crustacean.mass.per.g.plant), 
  "log10.crustacean.mass.per.g.plant"] = predict(crust.rf, y[is.na(y$log10.crustacean.mass.per.g.plant), ])
sum(is.na(y$log10.crustacean.mass.per.g.plant)) # 0 


# GAMMARID mesograzers
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.gammarid.mass.per.g.plant)) # 25
gamm.rf = randomForest(log10.gammarid.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude + 
                         Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness + 
                         log10.Zostera.AG.mass + log10.Zostera.shoots.core + log10.Zostera.sheath.width + log10.Zostera.sheath.length + 
                         log10.Zostera.longest.leaf.length + log10.grazer.richness.site + log10.AllelicRichness,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gammarid.mass.per.g.plant), 
  "log10.gammarid.mass.per.g.plant"] = predict(gamm.rf, y[is.na(y$log10.gammarid.mass.per.g.plant), ])
sum(is.na(y$log10.gammarid.mass.per.g.plant)) # 0


# GASTROPOD mesograzers
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.gastropod.mass.per.g.plant)) # 26
gast.rf = randomForest(log10.gastropod.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude + 
                         Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness + 
                         log10.Zostera.AG.mass + log10.Zostera.shoots.core + log10.Zostera.sheath.width + log10.Zostera.sheath.length + 
                         log10.Zostera.longest.leaf.length + log10.grazer.richness.site + log10.AllelicRichness,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gastropod.mass.per.g.plant), 
  "log10.gastropod.mass.per.g.plant"] = predict(gast.rf, y[is.na(y$log10.gastropod.mass.per.g.plant), ])
sum(is.na(y$log10.gastropod.mass.per.g.plant)) # 0


# total MESOGRAZERS
# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data$log10.mesograzer.mass.per.g.plant)) # 24
meso.rf = randomForest(log10.mesograzer.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude + 
                         Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness + 
                         log10.Zostera.AG.mass + log10.Zostera.shoots.core + log10.Zostera.sheath.width + log10.Zostera.sheath.length + 
                         log10.Zostera.longest.leaf.length + log10.grazer.richness.site + log10.AllelicRichness,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.mass.per.g.plant), 
  "log10.mesograzer.mass.per.g.plant"] = predict(meso.rf, y[is.na(y$log10.mesograzer.mass.per.g.plant), ])
sum(is.na(y$log10.mesograzer.mass.per.g.plant)) # 0


# PERIPHYTON
# NOTE: Here we need a reduced data set of 49 sites (i.e., excluding SW.A) for predictive model and imputation,
# because ALL plots from SW.A. had no periphyton values so it s not valid to impute values for that site.
ZEN_2014_master_data_49 <- droplevels(subset(y, Site != "SW.A"))
x <- ZEN_2014_master_data_49

# Use random forest to impute missing values. First build predictive model:
sum(is.na(ZEN_2014_master_data_49$log10.periphyton.mass.per.g.zostera)) # 16
peri.rf = randomForest(log10.periphyton.mass.per.g.zostera ~ Ocean + Coast + Latitude + Longitude + 
                         Temperature.C + Salinity.ppt + Day.length.hours + Depth.Categorical + AllelicRichness + 
                         log10.Zostera.AG.mass + log10.Zostera.shoots.core + log10.Zostera.sheath.width + log10.Zostera.sheath.length + 
                         log10.Zostera.longest.leaf.length + log10.grazer.richness.site + log10.AllelicRichness,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = x)

# Impute missing values
x[is.na(x$log10.periphyton.mass.per.g.zostera), 
  "log10.periphyton.mass.per.g.zostera"] = predict(peri.rf, x[is.na(x$log10.periphyton.mass.per.g.zostera), ])
sum(is.na(x$log10.periphyton.mass.per.g.zostera)) # 0


# Create data frame containing the imputed values and add them to the master data frame
names(y)
imputed.values.y <- y[c("Unique.ID",  "log10.Zostera.shoots.core", "log10.Zostera.AG.mass", "log10.Leaf.PercN", "log10.mesograzer.mass.per.g.plant", 
                        "log10.crustacean.mass.per.g.plant", "log10.gammarid.mass.per.g.plant", "log10.gastropod.mass.per.g.plant")]
names(imputed.values.y)

# Rename imputed values
colnames(imputed.values.y)[2:8] <- c("log10.Zostera.shoots.core.imputed", "log10.Zostera.AG.mass.imputed", "log10.Leaf.PercN.imputed", 
                                     "log10.mesograzer.mass.per.g.plant.imputed", "log10.crustacean.mass.per.g.plant.imputed", "log10.gammarid.mass.per.g.plant.imputed", 
                                     "log10.gastropod.mass.per.g.plant.imputed") 

# Create data frame containing the imputed value for periphyton and add to the 49-site data frame
names(x)
imputed.values <- x[c("Unique.ID", "log10.periphyton.mass.per.g.zostera")]
names(imputed.values)

# Rename iputed values
colnames(imputed.values)[2] <- c("log10.periphyton.mass.per.g.zostera.imputed" ) 
names(imputed.values)


# Integrate the imputed values back into master data set 
ZEN_2014_master_data_imputed <- ZEN_2014_master_data

ZEN_2014_master_data_imputed$log10.Zostera.shoots.core.imputed <- 
  imputed.values.y$log10.Zostera.shoots.core.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.Zostera.AG.mass.imputed <- 
  imputed.values.y$log10.Zostera.AG.mass.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.Leaf.PercN.imputed <- 
  imputed.values.y$log10.Leaf.PercN.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.mesograzer.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.mesograzer.mass.per.g.plant.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.crustacean.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.crustacean.mass.per.g.plant.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.gammarid.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.gammarid.mass.per.g.plant.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.gastropod.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.gastropod.mass.per.g.plant.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values.y$Unique.ID)]

ZEN_2014_master_data_imputed$log10.periphyton.mass.per.g.zostera.imputed <- 
  imputed.values$log10.periphyton.mass.per.g.zostera.imputed[match(ZEN_2014_master_data_imputed$Unique.ID, imputed.values$Unique.ID)]

# Create a data subset excluding SW.A (49 sites) to conduct analyses including periphyton 
ZEN_2014_master_data_49_imputed <- droplevels(subset(ZEN_2014_master_data_imputed, Site != "SW.A"))

names(ZEN_2014_master_data_imputed)
nrow(ZEN_2014_master_data_imputed) # 1000 - good

names(ZEN_2014_master_data_49_imputed)
nrow(ZEN_2014_master_data_49_imputed) # 980 - good


# How many observations are missing for each variable?
sum(is.na(ZEN_2014_master_data_imputed$log10.Zostera.AG.mass)) # 24
sum(is.na(ZEN_2014_master_data_imputed$log10.Zostera.AG.mass.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$log10.Zostera.shoots.core)) # 15
sum(is.na(ZEN_2014_master_data_imputed$log10.Zostera.shoots.core.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$Zostera.longest.leaf.length)) # 0
# sum(is.na(ZEN_2014_master_data_imputed$pct.cover.macroalgae)) # 20
# sum(is.na(ZEN_2014_master_data_imputed$pct.cover.seagrass)) # 20
sum(is.na(ZEN_2014_master_data_imputed$grazer.richness.site)) # 0
sum(is.na(ZEN_2014_master_data_imputed$log10.Leaf.PercN)) # 14
sum(is.na(ZEN_2014_master_data_imputed$log10.Leaf.PercN.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$Temperature.C)) # 0
sum(is.na(ZEN_2014_master_data_imputed$Salinity.ppt)) # 0
sum(is.na(ZEN_2014_master_data_imputed$Day.length.hours)) # 0
sum(is.na(ZEN_2014_master_data_imputed$Perc.Sand)) # 420
sum(is.na(ZEN_2014_master_data_imputed$GenotypicRichness)) # 0
sum(is.na(ZEN_2014_master_data_imputed$AllelicRichness)) # 0 
sum(is.na(ZEN_2014_master_data_imputed$Inbreeding)) # 0
sum(is.na(ZEN_2014_master_data_imputed$log10.periphyton.mass.per.g.zostera.imputed)) # 20 These are for site SW.A
sum(is.na(ZEN_2014_master_data_imputed$log10.mesograzer.mass.per.g.plant.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$log10.crustacean.mass.per.g.plant.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$log10.gastropod.mass.per.g.plant.imputed)) # 0
sum(is.na(ZEN_2014_master_data_imputed$amphipod.survival.24hr)) # 191


# Summary: We can now build and compare models that will have same number of observations
# BUT can't include these predictors because biased by missing from entire site: 
#   periphyton (missing from SW.A),
# Also can't include sediment grain size because missing from nearly half of plots. 

levels(ZEN_2014_master_data_49_imputed$Site) # 49 sites
# write.csv(ZEN_2014_master_data_49_imputed, "ZEN_2014_master_data_imputed_49_2016-05-23.csv", row.names = F)


###################################################################################
# EXPLORATORY PLOTS: PAIRS PANELS                                                 #
###################################################################################

# Environmental variables
pairs.panels(ZEN_2014_site_means[,c("Latitude","Longitude", "Temperature.C","Salinity.ppt", "log10.mean.fetch",
  "pop.density.2015", "Leaf.PercN.site", "log10.macrophytes.total.AG.mass.core.site", "log10.Zostera.shoots.core.site", 
  "log10.periphyton.mass.per.g.zostera.site","log10.crustacean.mass.per.g.plant.site",
  "amphipod.survival.24hr.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Seagrass above-ground biomass
pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
                                    "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "AllelicRichness", 
                                    "log10.periphyton.mass.per.g.zostera.site", "log10.crustacean.mass.per.g.plant.site", "grazer.richness.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Seagrass above-ground biomass (reduced)
pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.AG.mass.site", 
                                    "Latitude", "Temperature.C", "Leaf.PercN.site", 
                                    "log10.periphyton.mass.per.g.zostera.site", "log10.crustacean.mass.per.g.plant.site", "grazer.richness.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 10) 

# Seagrass morphology
pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.AG.mass.site", 
                                    "log10.Zostera.shoots.core.site", "log10.Zostera.longest.leaf.length.cm.site", "log10.Zostera.sheath.length.site",
                                    "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", 
                                    "log10.periphyton.mass.per.g.zostera.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 3) 
# NOTE: Zostera (log) sheath length is a substantially better predictor of AG biomass than is 
# longest leaf length. 

# Seagrass genetics and fitness
pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site", 
                                    "Zostera.sheath.width.site", "log10.Zostera.sheath.length.site", "log10.Zostera.longest.leaf.length.cm.site", 
                                    "Latitude", "Temperature.C", "Leaf.PercN.site", "log10.periphyton.mass.per.g.zostera.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 3) 

# Mesograzer biomass
pairs.panels(ZEN_2014_site_means[,c("log10.mesograzer.mass.per.g.plant.site", "log10.crustacean.mass.per.g.plant.site", 
                                    "log10.gastropod.mass.per.g.plant.site", "Latitude", "Temperature.C","Salinity.ppt", "log10.Zostera.shoots.core.site",
                                    "log10.Zostera.AG.mass.site", "Leaf.PercN.site", "AllelicRichness", 
                                    "log10.periphyton.mass.per.g.zostera.site", "amphipod.survival.24hr.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=1,scale=F, cex.cor = 3) 

# Crustacean mesograzer biomass
pairs.panels(ZEN_2014_site_means[,c("log10.crustacean.mass.per.g.plant.site", "log10.crustacean.abund.per.g.plant.site", 
                                    "log10.gastropod.abund.per.g.plant.site", "Latitude", "Temperature.C","Salinity.ppt", "log10.Zostera.shoots.core.site",
                                    "log10.Zostera.AG.mass.site", "Leaf.PercN.site", "AllelicRichness", 
                                    "log10.periphyton.mass.per.g.zostera.site", "amphipod.survival.24hr.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=1,scale=F, cex.cor = 3) 
# RESULT: crustacean B decreases with latitude, increases with Zostera B, increases with leaf % N.

# Gastropod mesograzer biomass
pairs.panels(ZEN_2014_site_means[,c("log10.gastropod.mass.per.g.plant.site", "Latitude", "Temperature.C","Salinity.ppt", "log10.Zostera.shoots.core.site",
                                    "log10.Zostera.AG.mass.site",  "Leaf.PercN.site", "AllelicRichness", 
                                    "log10.periphyton.mass.per.g.zostera.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=1,scale=F, cex.cor = 1) 
# RESULT: Gastropod B is not strongly related to any predictor except possibly DEcreasing with Zostera biomass.

# Periphyton biomass
pairs.panels(ZEN_2014_site_means[,c("log10.periphyton.mass.per.g.zostera.site", "log10.gastropod.mass.per.g.plant.site",  
                                    "log10.crustacean.mass.per.g.plant.site", "Latitude", "Temperature.C","Salinity.ppt", "log10.Zostera.shoots.core.site",
                                    "log10.Zostera.AG.mass.site",  "Leaf.PercN.site", "AllelicRichness", 
                                    "grazer.richness.site", "amphipod.survival.24hr.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=1,scale=F, cex.cor = 3) 
# RESULT: Hmm. Strongest predictor of periphyton biomass is grazer richness, which is POSITIVE. This is opposite
# what we see in the hierarchical model, where periphyton is negatively related to grazer richnes. However,
# hierarchical model excludes four sites so that may partially explain the difference. 

# Seagrass genetics
pairs.panels(ZEN_2014_site_means[,c("AllelicRichness", "log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
                                    "Leaf.PercN.site", "log10.Zostera.sheath.length.site", "log10.Zostera.sheath.width.site", "log10.Zostera.longest.leaf.length.cm.site",
                                    "Latitude", "Temperature.C", "Salinity.ppt", 
                                    "log10.periphyton.mass.per.g.zostera.site", "log10.crustacean.abund.per.g.plant.site", "grazer.richness.site")],
             smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 10) 

# Seagrass genetics
pairs.panels(ZEN_2014_site_means[,c("AllelicRichness", "GenotypicRichness", "Inbreeding", "log10.Zostera.AG.mass.site", 
                                    "log10.Zostera.shoots.core.site", "log10.Zostera.sheath.length.site", "log10.Zostera.longest.leaf.length.cm.site"
)], smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 10) 
# NOTE: Genotypic richness is highly correlated with bth allelic richness and inbreeding but the latter
# are relatively independent of one another. I suggest this means that if we use only one genetic
# predictor it should be genotypic richness (which relates to the other two) and if we use two 
# predictors we use allelic richness and inbreding, which are relatively independent predictors. One 
# concern is that genotypic richness is highly non-normal and cannot be fixed -- therefore, probably
# best to use allelic richness and inbreeding. 


###################################################################################
# BIVARIATE PLOTS                                                                 #
###################################################################################

# Zostera above-ground biomass x latitude 
ZAG.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.Zostera.AG.mass.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log Zostera above-ground dry mass (g) \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
ZAG.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.Zostera.AG.mass.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.2601, P = 0.00235, beta = -0.53  

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.Zostera.AG.mass.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0.2711, P = 0.01085, beta = -0.56  


# Zostera above-ground biomass x % N
ZAG.N = ggplot(ZEN_2014_site_means, aes(x = Leaf.PercN.site, y = log10.Zostera.AG.mass.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  Leaf % N") +  
  ylab("log Zostera above-ground dry mass (g) \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
ZAG.N + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


# Zostera above-ground biomass x mean fetch 
ZAG.fetch = ggplot(ZEN_2014_site_means, aes(x = log10.mean.fetch, y = log10.Zostera.AG.mass.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  log mean fetch") +  
  ylab("log Zostera above-ground dry mass (g) \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
ZAG.fetch + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


# Zostera shoot density x latitude 
Zshoots.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.Zostera.shoots.core.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log Zostera shoot density \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
# Zshoots.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.Zostera.shoots.core.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0, P = 0.8007, beta =   

# Pacific
y <- lm(ZEN_2014_site_means.Pac$log10.Zostera.shoots.core.site ~ ZEN_2014_site_means.Pac$Latitude)
summary(y) # Adjusted R-squared:  0, P = 0.3476  


# Zostera longest leaf length x latitude 
Zlength.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.Zostera.longest.leaf.length.cm.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log Zostera leaf length (cm) \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zlength.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.Zostera.longest.leaf.length.cm.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.2385, P = 0.003622, b = -0.52  

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.Zostera.longest.leaf.length.cm.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.7973, b = 0.06  


# Zostera allelic richness x latitude 
Zallele.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = AllelicRichness, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("Zostera AllelicRichness \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zallele.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$AllelicRichness) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.2498, P = 0.00289, b = -0.52  

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$AllelicRichness) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.88, b = 0.036  


# Periphyton biomass x latitude 
peri.lat = ggplot(ZEN_2014_site_means_49, aes(x = Latitude, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means_49$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log periphyton dry mass per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.323, P = 0.0007701, b = -0.58

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.97, b = 0.009  


# Periphyton biomass x %N 
peri.N = ggplot(ZEN_2014_site_means_49, aes(x = Leaf.PercN.site, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means_49$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  N availability (Leaf % N)") +  
  ylab("log periphyton dry mass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.N + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Periphyton biomass x leaf length
peri.length = ggplot(ZEN_2014_site_means_49, aes(x = log10.Zostera.longest.leaf.length.cm.site, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means_49$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log leaf length") +  
  ylab("log periphyton dry mass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.length + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


# Periphyton biomass x grazer richness 
peri.grich = ggplot(ZEN_2014_site_means, aes(x = log10.grazer.richness.site, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Grazer richness (site)") +  
  ylab("log periphyton dry mass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.grich + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Atl$log10.grazer.richness.site))
summary(x) # Adjusted R-squared:  0.08175, P = 0.07252, b = 0.33

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Pac$log10.grazer.richness.site))
summary(y) # Adjusted R-squared:  0, P = 0.9589, b = 0.01  



# Periphyton biomass x amphipod predation vulnerability 
peri.pred = ggplot(ZEN_2014_site_means, aes(x = amphipod.survival.24hr.site, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Amphipod predation (site)") +  
  ylab("log periphyton dry mass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.pred + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Atl$amphipod.survival.24hr.site))
summary(x) # Adjusted R-squared:  0.393, P = 0.0004767, b = 0.66

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.periphyton.mass.per.g.zostera.site) ~ scale(ZEN_2014_site_means.Pac$amphipod.survival.24hr.site))
summary(y) # Adjusted R-squared:  0, P = 0.4985, b = -0.19  


# Mesograzer biomass x latitude 
meso.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.mesograzer.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log mesograzer biomass per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
meso.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Gastropod abundance x latitude 
gast.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.gastropod.abund.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log gastropod abundance per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gast.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.gastropod.abund.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.4579, P = 2.422e-05, b = 0.69

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.gastropod.abund.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.8099, b = -0.06  



# Gastropod biomass x latitude 
gastB.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.gastropod.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log gastropod biomass per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gastB.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.gastropod.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.2219, P = 0.00503, b = 0.50

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.gastropod.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.5585, b = -0.01  


# Gastropod biomass x crustacean biomass 
gastB.crustB = ggplot(ZEN_2014_site_means, aes(x = log10.crustacean.mass.per.g.plant.site, y = log10.gastropod.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("log crustacean biomass per g macrophytes \n") +  
  ylab("log gastropod biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gastB.crustB + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Gastropod biomass x predation intensity (on amphipods)
gastB.pred = ggplot(ZEN_2014_site_means, aes(x = amphipod.survival.24hr.site, y = log10.gastropod.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("Predation intensity (on amphipods) \n") +  
  ylab("log gastropod biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gastB.pred + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


# Crustacean biomass x latitude 
crustB.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
crustB.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.crustacean.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.1519, P = 0.01902, b = -0.42

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.crustacean.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0, P = 0.9867, b = 0.004  


# Gammarid abundance x latitude 
gamm.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.gammarid.abund.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log gammarid abundance per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gamm.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.gammarid.abund.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.1604, P = 0.01626, b = -0.43

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.gammarid.abund.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0.1585, P = 0.0463, b = -0.45  


# Gammarid biomass x latitude 
gammB.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.gammarid.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log gammarid biomass per g macrophytes \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gammB.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.gammarid.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.1174, P = 0.03592, b = -0.38

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.gammarid.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0.04422, P = 0.1873, b = -0.31  


# Mesograzer richness (site-level) x latitude 
grich.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = grazer.richness.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("mesograzer species richness (site-level) \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
grich.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$grazer.richness.site) ~ scale(ZEN_2014_site_means.Atl$Latitude))
summary(x) # Adjusted R-squared:  0.1612, P = 0.016, b = -0.43

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$grazer.richness.site) ~ scale(ZEN_2014_site_means.Pac$Latitude))
summary(y) # Adjusted R-squared:  0.04599, P = 0.1832, b = -0.31  


# For some reason the dataframe had a bogus empty level of Ocean ...
ZEN_2014_site_means = droplevels(ZEN_2014_site_means)


# Crustacean biomass x Zostera AG biomass 
crustB.ZAG = ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.AG.mass.site, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Zostera above-ground mass") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,100))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
crustB.ZAG + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Crustacean biomass x Zostera shoot density
crustB.shoots = ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.shoots.core.site, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Zostera shoot density") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,100))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
crustB.shoots + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Zostera allelic richness x latitude 
Zallrich.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = AllelicRichness, col = Ocean)) +
  geom_point(size = 5, shape = 16) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("Zostera allelic richness (site-level) \n") +  
  scale_x_continuous(limits = c(30,72)) +
  #   scale_y_continuous(limits=c(0,3)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zallrich.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Zostera genotypic richness x latitude 
Zgrich.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = GenotypicRichness, col = Ocean)) +
  geom_point(size = 5, shape = 16) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("Zostera genotypic richness (site-level) \n") +  
  scale_x_continuous(limits = c(30,72)) +
  #   scale_y_continuous(limits=c(0,3)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zgrich.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Zostera inbreeding x latitude 
Zgrich.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = GenotypicRichness, col = Ocean)) +
  geom_point(size = 5, shape = 16) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("Zostera inbreeding \n") +  
  scale_x_continuous(limits = c(30,72)) +
  #   scale_y_continuous(limits=c(0,3)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zgrich.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Zostera sheath length vs above-ground biomass
sheath.ZAG = ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.sheath.length.site, y = log10.Zostera.AG.mass.site, col = Ocean)) +
  geom_point(size = 5, shape = 16) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Zostera sheath length") +  
  ylab("log Zostera above-ground biomass \n") +  
  # scale_x_continuous(limits = c(30,72)) +
  # scale_y_continuous(limits=c(0,3)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
sheath.ZAG + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Periphyton biomass x Gastropod abundance 
peri.gast = ggplot(ZEN_2014_site_means, aes(x = log10.gastropod.abund.per.g.plant.site, y = log10.periphyton.mass.per.g.zostera.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log gastropod abundance per g macrophytes") +  
  ylab("log periphyton biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
peri.gast + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Zostera allelic richness x crustacean abundance 
Zall.crust = ggplot(ZEN_2014_site_means, aes(x = AllelicRichness, y = log10.crustacean.abund.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Zostera allelic richness (site)") +  
  ylab("log crustacean abundance per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zall.crust + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)

# Zostera allelic richness x crustacean biomass 
Zall.crustB = ggplot(ZEN_2014_site_means, aes(x = AllelicRichness, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n  log Zostera allelic richness (site)") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zall.crustB + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


# Crustacean biomass x Zostera leaf %N 
pctN.crust = ggplot(ZEN_2014_site_means, aes(x = log10.Leaf.PercN.site, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n   log Zostera leaf % N") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
pctN.crust + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


# Atlantic
x <- lm(scale(ZEN_2014_site_means.Atl$log10.crustacean.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Atl$log10.Leaf.PercN.site))
summary(x) # Adjusted R-squared:  0.03275, P = 0.1702, b = 0.26

# Pacific
y <- lm(scale(ZEN_2014_site_means.Pac$log10.crustacean.mass.per.g.plant.site) ~ scale(ZEN_2014_site_means.Pac$log10.Leaf.PercN.site))
summary(y)# Adjusted R-squared:  0.1206, P = 0.07376, b = 0.41


# Crustacean biomass x predation on amphipods
pred.crust = ggplot(ZEN_2014_site_means, aes(x = amphipod.survival.24hr.site, y = log10.crustacean.mass.per.g.plant.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n   Predation on amphipods") +  
  ylab("log crustacean biomass per g macrophytes \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
pred.crust + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)

# Crustacean mean body mass x latitude
# NOTE: Need to remove sites with zero values before doing this calculation!
crustBN = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.crustacean.mean.body.mass.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n   Latitude") +  
  ylab("log crustacean mean body mass per \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
crustBN + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)

# Gastropod mean body mass x latitude
# NOTE: Need to remove sites with zero values before doing this calculation!
gastBN = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.gastropod.mean.body.mass.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue","forestgreen")) +
  xlab("\n   Latitude") +  
  ylab("log gastropod mean body mass per \n") +  
  # scale_x_continuous(limits = c(30,72))+
  # scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
gastBN + geom_smooth(method = lm, fullrange = F, se = F, lwd = 0, col="black", na.rm=T)


###################################################################################
# IDENTIFY BEST EXPLANATORY VARIABLES USING RANDOM FORESTS - SITE MEANS           #
###################################################################################

# Construct random forest: ZOSTERA ABOVE-GROUND BIOMASS (WITH biogeographic variables)
RF.log10.Zostera.AG.mass.site.1 = randomForest(
  log10.Zostera.AG.mass.site ~ 
    Ocean
  + Coast
  + Latitude
  + Longitude
  + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  #   + log10.Zostera.AG.mass
  #   + log10.macrophytes.total.AG.mass.core
  #   + log10.Zostera.shoots.core
  #   + Zostera.sheath.width
  #   + Zostera.sheath.length
  #   + Zostera.longest.leaf.length
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site  
  + log10.mesograzer.mass.per.g.plant.site
  + log10.crustacean.mass.per.g.plant.site
  #   + log10.amphipod.mass.per.g.plant
  + log10.gammarid.mass.per.g.plant.site
  + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.log10.Zostera.AG.mass.site.1 # % Var explained: 50.09

# Plot error as a function of # of trees
plot(RF.log10.Zostera.AG.mass.site.1) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(RF.log10.Zostera.AG.mass.site.1) # Look at % Inc MSE (first plot) not Inc Node Purity

# RESULT: Latitude is by far the best predictor. Then coast, day length, salinity, depth. 


# Construct random forest: ZOSTERA ABOVE-GROUND BIOMASS (biogeographic variables omitted)
RF.log10.Zostera.AG.mass.site.2 = randomForest(
  log10.Zostera.AG.mass.site ~ 
    # Ocean
    # + Coast
    # + Latitude
    # + Longitude
    + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  #   + log10.Zostera.AG.mass
  #   + log10.macrophytes.total.AG.mass.core
  #   + log10.Zostera.shoots.core
  #   + Zostera.sheath.width
  #   + Zostera.sheath.length
  #   + Zostera.longest.leaf.length
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site  
  + log10.mesograzer.mass.per.g.plant.site
  + log10.crustacean.mass.per.g.plant.site
  #   + log10.amphipod.mass.per.g.plant
  + log10.gammarid.mass.per.g.plant.site
  + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.log10.Zostera.AG.mass.site.2 # % Var explained: 34.65 - much less than when biogeographic variables are included

# Plot error as a function of # of trees
plot(RF.log10.Zostera.AG.mass.site.2) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(RF.log10.Zostera.AG.mass.site.2) # Look at % Inc MSE (first plot) not Inc Node Purity

# RESULT: When biogeographic variables (Ocean, Latotude, Longitude, Coast) are omitted, the highest scores
# are day length (proxy for latitude), crustacean mass, and allelic richness. 

plot(ZEN_2014_site_means$log10.Zostera.AG.mass.site ~ ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site)
x <- lm(ZEN_2014_site_means$log10.Zostera.AG.mass.site ~ ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site)
summary(x) # Adjusted R-squared:  0.166, P = 0.00214 **  

plot(ZEN_2014_site_means$log10.Zostera.AG.mass.site ~ ZEN_2014_site_means$AllelicRichness)
y <- lm(ZEN_2014_site_means$log10.Zostera.AG.mass.site ~ ZEN_2014_site_means$AllelicRichness)
summary(y) # Adjusted R-squared:  0.1459, P = 0.00361 ** 

z <- lm(log10.Zostera.AG.mass.site ~ AllelicRichness + log10.periphyton.mass.per.g.zostera.site, data = ZEN_2014_site_means)
summary(z) # Adjusted R-squared:  0.2319, Allelic P = 0.0298 *, periphyton P = 0.0171 * 


# Construct random forest: PERIPHYTON BIOMASS
RF.periphyton.mass.1 = randomForest(
  log10.periphyton.mass.per.g.zostera.site ~ 
    Ocean
  + Coast
  + Latitude
  + Longitude
  + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass.site
  + log10.macrophytes.total.AG.mass.core.site
  + log10.Zostera.shoots.core.site
  #   + Zostera.sheath.width
  #   + Zostera.sheath.length
  + log10.Zostera.longest.leaf.length.cm.site
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  # + log10.periphyton.mass.per.g.zostera.site  
  + log10.mesograzer.mass.per.g.plant.site
  + log10.crustacean.mass.per.g.plant.site
  + log10.amphipod.mass.per.g.plant.site
  + log10.gammarid.mass.per.g.plant.site
  + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.periphyton.mass.1 # % Var explained: 39.34

# Plot error as a function of # of trees
plot(RF.periphyton.mass.1) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(RF.periphyton.mass.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULT: Strongest predictors of periphyton, even with biogeographic variables in model, are 
# Zostera leaf length and predation rate on amphipods. Trophic cascade?

plot(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
pred.peri <- lm(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
summary(pred.peri) # Adjusted R-squared:  0.1662, P = 0.00474 ** 
# Positive relationship between predation pressure on amphipods and periphyton. 

pred.leaf.peri <- lm(log10.periphyton.mass.per.g.zostera.site ~ amphipod.survival.24hr.site + log10.Zostera.longest.leaf.length.cm.site, data = ZEN_2014_site_means)
summary(pred.leaf.peri) # Adjusted R-squared:  0.4072, predation P = 0.016819 * 
# Positive effect of predation remians when Zostera leaf length is in model. 

pred.lat.peri <- lm(log10.periphyton.mass.per.g.zostera.site ~ amphipod.survival.24hr.site + Latitude, data = ZEN_2014_site_means)
summary(pred.lat.peri) 
#                             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                 -0.31068    0.84513  -0.368   0.7152  
# amphipod.survival.24hr.site  0.69960    0.50601   1.383   0.1749  
# Latitude                    -0.02621    0.01284  -2.041   0.0482 *
  

# Construct random forest: CRUSTACEAN BIOMASS (WITH biogeographic variables)
RF.crustacean.mass.1 = randomForest(
  log10.crustacean.mass.per.g.plant.site ~ 
    Ocean
  + Coast
  + Latitude
  + Longitude
  + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass.site
  + log10.macrophytes.total.AG.mass.core.site
  + log10.Zostera.shoots.core.site
  #   + Zostera.sheath.width
  #   + Zostera.sheath.length
  + log10.Zostera.longest.leaf.length.cm.site
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site  
  # + log10.mesograzer.mass.per.g.plant.site
  # + log10.crustacean.mass.per.g.plant.site
  # + log10.amphipod.mass.per.g.plant.site
  # + log10.gammarid.mass.per.g.plant.site
  + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  # + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.crustacean.mass.1 # % Var explained: 37.06

# Plot error as a function of # of trees
plot(RF.crustacean.mass.1) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(RF.crustacean.mass.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULT: Crustacean biomass is basically driven entirely by geographic variation (long, coast)


# Construct random forest: CRUSTACEAN BIOMASS (omitting biogeographic variables)
RF.crustacean.mass.2 = randomForest(
  log10.crustacean.mass.per.g.plant.site ~ 
    # Ocean
    # + Coast
    # + Latitude
    # + Longitude
    + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass.site
  + log10.macrophytes.total.AG.mass.core.site
  + log10.Zostera.shoots.core.site
  #   + Zostera.sheath.width
  #   + Zostera.sheath.length
  + log10.Zostera.longest.leaf.length.cm.site
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site  
  # + log10.mesograzer.mass.per.g.plant.site
  # + log10.crustacean.mass.per.g.plant.site
  # + log10.amphipod.mass.per.g.plant.site
  # + log10.gammarid.mass.per.g.plant.site
  + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  # + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.crustacean.mass.2 # % Var explained: 17.87 - much less

# Plot error as a function of # of trees
plot(RF.crustacean.mass.2) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(RF.crustacean.mass.2) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULT: When geographic variables are omitted, leaf % N comes in just after day length. Fits with
# results of ZEN 1. No effect of allelic richness, however. 

plot(log10.crustacean.mass.per.g.plant.site ~ log10.Leaf.PercN.site, data = ZEN_2014_site_means)
crust.N <- lm(log10.crustacean.mass.per.g.plant.site ~ log10.Leaf.PercN.site, data = ZEN_2014_site_means)
summary(crust.N) # Adjusted R-squared:  0.1469, P = 0.00349 ** 

crust.lat.N <- lm(log10.crustacean.mass.per.g.plant.site ~ log10.Leaf.PercN.site + Latitude, data = ZEN_2014_site_means)
summary(crust.lat.N) # Adjusted R-squared:  0.3273, leaf N: P = 0.001073 ** 
# RESULT: Leaf % N still strongly significant predictor of crustacean B even when latitude is in model.



plot(ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
pred.gast <- lm(log10.gastropod.mass.per.g.plant.site ~ amphipod.survival.24hr.site, data = ZEN_2014_site_means)
summary(pred.gast) # Adjusted R-squared:  0.1168, P = 0.01527

plot(ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site ~ ZEN_2014_site_means$gastropod.survival.24hr.site)
pred.gast.2 <- lm(ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site ~ ZEN_2014_site_means$gastropod.survival.24hr.site)
summary(pred.gast.2) # P = 0.068964 .  Too little variation in gastropod predation (mostly zero)

plot(ZEN_2014_site_means$log10.gammarid.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
pred.gam <- lm(ZEN_2014_site_means$log10.gammarid.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
summary(pred.gam)

plot(ZEN_2014_site_means$log10.amphipod.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
pred.amph <- lm(ZEN_2014_site_means$log10.amphipod.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
summary(pred.amph)

plot(ZEN_2014_site_means$log10.crustacean.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
pred.crust <- lm(ZEN_2014_site_means$log10.crustacean.mass.per.g.plant.site ~ ZEN_2014_site_means$amphipod.survival.24hr.site)
summary(pred.crust)

plot(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site ~ ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site)
peri.gast <- lm(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site ~ ZEN_2014_site_means$log10.gastropod.mass.per.g.plant.site)
summary(peri.gast) # P = 0.035


# Construct random forest: GRAZER RICHNESS
RF.grich.1 = randomForest(
  log10.grazer.richness.site ~ 
    Ocean
  + Coast
  + Latitude
  + Longitude
  + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass.site
  + log10.macrophytes.total.AG.mass.core.site
  + log10.Zostera.shoots.core.site
  + Zostera.sheath.width.site
  + Zostera.sheath.length.site
  + log10.Zostera.longest.leaf.length.cm.site
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site 
  # + log10.mesograzer.abund.per.g.plant.site
  # + log10.mesograzer.mass.per.g.plant.site
  # + log10.crustacean.mass.per.g.plant.site
  # + log10.amphipod.mass.per.g.plant.site
  # + log10.gammarid.mass.per.g.plant.site
  # + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  # + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.grich.1 # % Var explained: 25.79

# Plot error as a function of # of trees
plot(RF.grich.1) # Good.

# Plot variable importance
varImpPlot(RF.grich.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULT: Longitude dominates.

plot(log10.grazer.richness.site ~ log10.periphyton.mass.per.g.zostera.site, data = ZEN_2014_site_means)
grich.peri <- lm(log10.grazer.richness.site ~ log10.periphyton.mass.per.g.zostera.site, data = ZEN_2014_site_means)
summary(grich.peri) # Adjusted R-squared:  0.1544, P =   0.00305 **

grich.lat.peri <- lm(log10.grazer.richness.site ~ log10.periphyton.mass.per.g.zostera.site + Latitude, data = ZEN_2014_site_means)
summary(grich.lat.peri) # Adjusted R-squared:  0.2705, periphyton P =   0.16278


# Construct random forest: ALLELIC RICHNESS
RF.zrich.1 = randomForest(
  AllelicRichness ~ 
    Ocean
  + Coast
  + Latitude
  + Longitude
  + log10.mean.fetch
  + Temperature.C
  + Salinity.ppt
  + Day.length.hours
  + Depth.Categorical
  + pop.density.2015
  #  + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  # + AllelicRichness
  # + Inbreeding
  # + log10.Zostera.AG.mass.site
  # + log10.macrophytes.total.AG.mass.core.site
  # + log10.Zostera.shoots.core.site
  # + Zostera.sheath.width.site
  # + Zostera.sheath.length.site
  # + log10.Zostera.longest.leaf.length.cm.site
  #   + log10.pct.cover.macroalgae
  #   + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera.site 
  # + log10.mesograzer.abund.per.g.plant.site
  # + log10.mesograzer.mass.per.g.plant.site
  # + log10.crustacean.mass.per.g.plant.site
  # + log10.amphipod.mass.per.g.plant.site
  # + log10.gammarid.mass.per.g.plant.site
  # + log10.gastropod.mass.per.g.plant.site
  # + log10.grazer.richness.plot
  + log10.grazer.richness.site
  + amphipod.survival.24hr.site
  + log10.Leaf.PercN.site
  ,
  data = ZEN_2014_site_means,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.zrich.1 # % Var explained: 20.88

# Plot error as a function of # of trees
plot(RF.zrich.1) # Hmm

# Plot variable importance
varImpPlot(RF.zrich.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULT: Best predictors are human population density (!) and temperature 

plot(AllelicRichness ~ Temperature.C, data = ZEN_2014_site_means)
zrich.temp <- lm(AllelicRichness ~ Temperature.C , data = ZEN_2014_site_means)
summary(zrich.temp) # Adjusted R-squared:  0.201, Temperature P = 0.000644 ***

plot(AllelicRichness ~ amphipod.survival.24hr.site, data = ZEN_2014_site_means)
zrich.temp.pred <- lm(AllelicRichness ~ Temperature.C + amphipod.survival.24hr.site, data = ZEN_2014_site_means)
summary(zrich.temp.pred) # Adjusted R-squared:  0.1402, Temperature P = 0.0117 * (no predation effect)


###################################################################################
# PERFORM VIF ANALYSIS: ZOSTERA ABOVE-GROUND MASS (NO OTHER ZOSTERA PREDICTORS)   #
###################################################################################

# NOTE: After exploring this a bit, I have found that too severe pruning based on VIF
# can remove important variables. For example, in third stage of this analysis for 
# Zostera AG biomass, Inbreeding comes up with VIF = 2.940177. Using a strict threshold
# of VIF > 2 to delete eliminates this. But leaving it in shows that inbreeding is among the 
# most highly significant predictors of Zostera biomass. So there seems to be no 
# substitute in this process (as in modeling these complex data sets generally) for a 
# certain amount of judgment . . . 


vifZAG.1 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                + log10.mean.fetch
                + pop.density.2015
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                # + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                + Inbreeding
                + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                # + Perc.Silt
                # + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.1)
# 

vifZAG.2 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                + Inbreeding
                + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                + Perc.Silt
                + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.2)
# Temperature.C = 14.356929 But this is high in RF and an important variable mechanistically. 
# Let's remove next highest: Perc.Gravel = 14.132502 


vifZAG.3 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                + Inbreeding
                + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                + Perc.Silt
                # + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.3)
# Temperature.C = 5.943901 But this is high in RF and a major environmental variable. 
# Let's remove next highest: GenotypicRichness = 4.047051 


vifZAG.4 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                + Inbreeding
                # + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                + Perc.Silt
                # + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.4)
# Perc.Silt = 2.332913


vifZAG.5 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                + Inbreeding
                # + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                # + Perc.Silt
                # + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.5)
# AllelicRichness = 3.181334. But this is a more informative and statistically well-behaved
# variable than inbreeding, so let's remove the latter. 

vifZAG.6 <- lme(log10.Zostera.AG.mass ~ 
                  + Temperature.C
                + Salinity.ppt
                + Day.length.hours
                # + log10.Zostera.AG.mass
                # + log10.pct.cover.seagrass
                + log10.pct.cover.macroalgae
                # + log10.Zostera.shoots.core
                # + Zostera.sheath.width
                # + Zostera.sheath.length
                # + Zostera.longest.leaf.length
                + AllelicRichness
                # + Inbreeding
                # + GenotypicRichness
                + log10.periphyton.mass.per.g.zostera
                # + Perc.Sand
                # + Perc.Silt
                # + Perc.Gravel
                + Leaf.PercN
                + log10.crustacean.abund.per.g.plant
                + log10.gastropod.abund.per.g.plant
                + grazer.richness.site
                , random = ~1 | Ocean/Site.Code
                , na.action = na.omit
                , data = ZEN_2014_master_data)
vif(vifZAG.6)
# All variables now have VIF < 2.


###################################################################################
# PERFORM VIF ANALYSIS: PERIPHYTON BIOMASS                                        #
###################################################################################

# NOTE that Sweden subsite A did not collect info on periphyton and should be omitted from 
# analyses that consider periphyton:
levels(ZEN_2014_master_data_NO_SW.A$Site)

vifperi.1 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 + Perc.Silt
                 + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.1)
# Perc.Silt = 591.904801 


vifperi.2 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 # + Perc.Silt
                 + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.2)
# Temperature.C = 23.588716. Remove next highest: Perc.Gravel = 17.452378


vifperi.3 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.3)
# Temperature.C = 24.138831. Next highest: GenotypicRichness = 14.652632


vifperi.4 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 # + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.4)
# Salinity.ppt = 4.469234


vifperi.5 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 # + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 # + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.5)
# Temperature.C = 3.403308. Next highest: Zostera.sheath.width = 3.062357


vifperi.6 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 # + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 # + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 # + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.6)
# Perc.Sand = 2.950242


vifperi.7 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 # + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 # + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 + Inbreeding
                 # + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 # + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.7)
# AllelicRichness = 3.371488. Next highest: Inbreeding = 2.437403


vifperi.8 <- lme(log10.periphyton.mass.per.g.zostera ~ 
                   + Temperature.C
                 # + Salinity.ppt
                 + Day.length.hours
                 + log10.Zostera.AG.mass
                 + log10.pct.cover.seagrass
                 + log10.pct.cover.macroalgae
                 + log10.Zostera.shoots.core
                 # + Zostera.sheath.width
                 + Zostera.sheath.length
                 + Zostera.longest.leaf.length
                 + AllelicRichness
                 # + Inbreeding
                 # + GenotypicRichness
                 # + log10.periphyton.mass.per.g.zostera
                 # + Perc.Sand
                 # + Perc.Silt
                 # + Perc.Gravel
                 + Leaf.PercN
                 + log10.crustacean.abund.per.g.plant
                 + log10.gastropod.abund.per.g.plant
                 + grazer.richness.site
                 , random = ~1 | Ocean/Site.Code
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_NO_SW.A)
vif(vifperi.8)
# All remaining variables have VIF < 2.


###################################################################################
# PERFORM VIF ANALYSIS: CRUSTACEAN ABUNDANCE                                      #
###################################################################################

vifcrust.1 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  + Perc.Sand
                  + Perc.Silt
                  + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.1)
# Perc.Sand = 317.796583


vifcrust.2 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  # + Perc.Sand
                  + Perc.Silt
                  + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.2)
# Temperature.C = 14.265183 Next in line is Day.length.hours = 11.153726. But close third, 
# which we delete, is: Perc.Gravel = 10.864142


vifcrust.3 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  # + Perc.Sand
                  + Perc.Silt
                  # + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.3)
# Temperature.C = 12.022051. Next: GenotypicRichness = 6.992560


vifcrust.4 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  # + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  # + Perc.Sand
                  + Perc.Silt
                  # + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.4)
# Temperature.C = 3.180547. Next: Zostera.sheath.width = 2.895072


vifcrust.5 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  # + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  # + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  # + Perc.Sand
                  + Perc.Silt
                  # + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.5)
# log10.pct.cover.macroalgae = 2.507968.  Temperature.C = 2.498815. Delete Day.length.hours = 2.477936


vifcrust.6 <- lme(log10.crustacean.abund.per.g.plant ~ 
                    + Temperature.C
                  + Salinity.ppt
                  # + Day.length.hours
                  + log10.Zostera.AG.mass
                  + log10.pct.cover.seagrass
                  + log10.pct.cover.macroalgae
                  + log10.Zostera.shoots.core
                  # + Zostera.sheath.width
                  + Zostera.sheath.length
                  + Zostera.longest.leaf.length
                  + AllelicRichness
                  + Inbreeding
                  # + GenotypicRichness
                  + log10.periphyton.mass.per.g.zostera
                  # + Perc.Sand
                  + Perc.Silt
                  # + Perc.Gravel
                  + Leaf.PercN
                  # + log10.crustacean.abund.per.g.plant
                  # + log10.gastropod.abund.per.g.plant
                  # + grazer.richness.site
                  , random = ~1 | Ocean/Site.Code
                  , na.action = na.omit
                  , data = ZEN_2014_master_data)
vif(vifcrust.6)
# Highest now is Salinity.ppt = 2.189858. Let's keep the rest in 


###################################################################################
# HIERARCHICAL MODELING APPROACH: RATIONALE                                       #
###################################################################################

# Some notes on modeling approach: 

# 1) The random forest analyses above confirm that latitude and
# longitude have big, overriding effects on all variables examined. This is not surprising
# but also not very helpful since these are  non-mechanistic predictors. That is, latitude 
# is a proxy for a whole bunch of covarying factors (many of which we actually have 
# measured) and longitude is arguably a proxy for the effect of evolutionary history 
# resulting in phylogenetically differentiated biotas in different oceans. Therefore in the following
# models (a) I ignore latitude and instead include the environmental variables covarying 
# with it that have mechanistic bases, and (b) to account for longitude I use Ocean as
# a random grouping variable in the hierarchical model. 

# 2) In contrast to previous versions of this script, I now use "Site" (n = 50) as the grouping variable
# rather than "Site.Code" (n = 25) because  Site is the level at which "top-level" variables are measured,
# including temperature, salinity, etc. as well as alleic and grazer richness. 

# 3) Remember that, as in the ZEN 1 analysis, submodels of site-level variables as the response must 
# use standard linear model (lm) rather than hierarchical (lme) since we only have a single value per site 
# for those variables (e.g., allelic richness, grazer richness). Morever, these lm's must use the site mean
# data frame (N = 50) rather than the mian data frame (N = 1000) since there are only 50 independent observations.

# IN FACT, it sems reasonable to use simple linear models iof site means for most analyses since 
# most analyses are using mostly or excvlusively site-level predictors. 

# 4) For the purposes of model comparisons using AIC it is essential that all models use the 
# same data set. This in turn means that it is essential to impute missing values since, if this is 
# not done,  fitting a model that omits a particular predictor variable can change the total number 
# of observations (plots) in the analysis if data were missing for that variable, and hence the DF of
# the model, which will change the AIC for reasons unrelated to model fit.

# 5) The starting approach, after eliminating latitude and longitude, was to begin with all 
# variables in the data set and use VIF to prune out highly collinear variables, with one 
# caveat: variables scoring high in the RF analyses, and those with clear mechanistic
# links to the response (e.g., temperature) were retained even if VIF analysis 
# suggested collinearity. In these cases I would move to the variable with next highest VIF score 
# and delete that one.  

# 6) I looked at sets of closely related variables to try to minimize redundancy. Thus, 
# (a) I eliminated genotypic richness from predictors because it was correlated with each of 
# allelic richness and inbreeding (which were relatively uncorrelated with one another) and it is
# highly non-normally distributed. After consultation wiht Jay I also eliminated inbreeding
# since it is highly non-normal and allelic richness is a better metric of genetic diversity. 
# Similarly, I eliminated (b) seagrass sheath width and length because they are correlated with 
# longest leaf length and the latter seems more clearly mechanistically
# related to most response variables. (c) I eliminated % silt and % gravel and retained % sand, 
# since these are inherently correlated (they add up to one) and sand was the largest fraction. 
# IN FACT, I ELIMINATED ALL SEDIMENT VARIABLES SINCE NEARLY HAFL THE PLOTS HAD NO MEASUREMENTS. 

# 7) I eliminated as predictors variables that seem likely to be correlated errors rather
# than causal predictors. So, for example, when modeling crustacean biomass I eliminated 
# gastropod biomass as a predictor. I also eliminated grazer richness as a predictor UNTIL we 
# can get a Chao-Jost estimate of richness that corrects for number of individuls sampled. 

# 8) Finally, I eliminated as predictors those varables with highly non-normal distributions, like 
# % cover of macroalgae.


###################################################################################
# BUILD HIERARCHICAL MODELS: ZOSTERA ABOVE-GROUND BIOMASS                         #
###################################################################################

# Using data for 49 sites to include periphyton as predictor.

names(ZEN_2014_master_data_49_imputed)

ZosteraAG1.49 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    + Latitude
                  + Temperature.C
                  + Salinity.ppt
                  + glaciation
                  + log10.mean.fetch 
                  + pop.density.2015 
                  + Day.length.hours
                  # + log10.Zostera.AG.mass.imputed
                  # + log10.Zostera.shoots.core.imputed
                  # + log10.Zostera.sheath.length
                  # + log10.Zostera.longest.leaf.length
                  + AllelicRichness
                  + log10.periphyton.mass.per.g.zostera.imputed 
                  + log10.Leaf.PercN.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  + log10.crustacean.mass.per.g.plant.imputed
                  + log10.gastropod.mass.per.g.plant.imputed
                  + grazer.richness.site
                  , random = ~1 | Ocean/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZosteraAG1.49)
# N = 980 observations
# AIC = -22.09148

sem.coefs(list(ZosteraAG1.49), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                         response                                   predictor     estimate  std.error p.value
# 1  log10.Zostera.AG.mass.imputed                                    Latitude -0.553178142 0.22051373  0.0166
# 12 log10.Zostera.AG.mass.imputed    log10.gastropod.mass.per.g.plant.imputed -0.056529067 0.03227852  0.0802
# 3  log10.Zostera.AG.mass.imputed                                Salinity.ppt -0.222226472 0.13757367  0.1147
# 2  log10.Zostera.AG.mass.imputed                               Temperature.C  0.127853880 0.13872587  0.3627
# 8  log10.Zostera.AG.mass.imputed                             AllelicRichness  0.061895762 0.12130290  0.6129
# 9  log10.Zostera.AG.mass.imputed log10.periphyton.mass.per.g.zostera.imputed  0.023236785 0.04722421  0.6228
# 10 log10.Zostera.AG.mass.imputed                    log10.Leaf.PercN.imputed -0.017342170 0.03688868  0.6384
# 7  log10.Zostera.AG.mass.imputed                            Day.length.hours  0.069515573 0.16123473  0.6689
# 11 log10.Zostera.AG.mass.imputed   log10.crustacean.mass.per.g.plant.imputed  0.014764693 0.03552011  0.6777
# 5  log10.Zostera.AG.mass.imputed                            log10.mean.fetch  0.027158048 0.11865727  0.8202
# 13 log10.Zostera.AG.mass.imputed                        grazer.richness.site  0.021546266 0.11049867  0.8465
# 6  log10.Zostera.AG.mass.imputed                            pop.density.2015 -0.016242407 0.11383807  0.8873
# 4  log10.Zostera.AG.mass.imputed                                  glaciation  0.007867825 0.18165498  0.9657


# Since there is no hint of periphyton or human population density effects let's go back to 
# the complete data set (50 sites) and take these out.

ZosteraAG1 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    + Latitude
                  + Temperature.C
                  + Salinity.ppt
                  + glaciation
                  + log10.mean.fetch 
                  # + pop.density.2015 # one site missing: IR.B
                  + Day.length.hours
                  # + log10.Zostera.AG.mass.imputed
                  # + log10.Zostera.shoots.core.imputed
                  # + log10.Zostera.sheath.length
                  # + log10.Zostera.longest.leaf.length
                  + AllelicRichness
                  # + log10.periphyton.mass.per.g.zostera.imputed 
                  # + log10.Leaf.PercN.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  + log10.crustacean.mass.per.g.plant.imputed
                  + log10.gastropod.mass.per.g.plant.imputed
                  + grazer.richness.site
                  , random = ~1 | Ocean/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_imputed)
summary(ZosteraAG1)

sem.coefs(list(ZosteraAG1), ZEN_2014_master_data_imputed, standardize = "scale")
#                         response                                 predictor    estimate  std.error p.value
# 9  log10.Zostera.AG.mass.imputed  log10.gastropod.mass.per.g.plant.imputed -0.06850580 0.03131236  0.0289
# 1  log10.Zostera.AG.mass.imputed                                  Latitude -0.47601771 0.21131106  0.0298

# No geography
ZosteraAG2 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    # + Latitude
                  + Temperature.C
                  + Salinity.ppt
                  + glaciation
                  + log10.mean.fetch 
                  # + pop.density.2015 # one site missing: IR.B
                  # + Day.length.hours
                  # + log10.Zostera.AG.mass.imputed
                  # + log10.Zostera.shoots.core.imputed
                  # + log10.Zostera.sheath.length
                  # + log10.Zostera.longest.leaf.length
                  + AllelicRichness
                  # + log10.periphyton.mass.per.g.zostera.imputed 
                  # + log10.Leaf.PercN.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  + log10.crustacean.mass.per.g.plant.imputed
                  + log10.gastropod.mass.per.g.plant.imputed
                  + grazer.richness.site
                  , random = ~1 | Ocean/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_imputed)
summary(ZosteraAG2)
sem.coefs(list(ZosteraAG2), ZEN_2014_master_data_imputed, standardize = "scale")


# Try an equivalent model using site means
names(ZEN_2014_site_means)

ZAG.site <- lm(log10.Zostera.AG.mass.site ~ Latitude, data = ZEN_2014_site_means)
summary(ZAG.site)

#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.95759    0.17006  17.391  < 2e-16 ***
#   Latitude    -0.02071    0.00371  -5.583 1.08e-06 ***
# 
# Residual standard error: 0.2638 on 48 degrees of freedom
# Multiple R-squared:  0.3937,	Adjusted R-squared:  0.3811 
# F-statistic: 31.17 on 1 and 48 DF,  p-value: 1.08e-06


###################################################################################
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: RANDOM FORESTS                         #
###################################################################################

# Explore graphically

pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.AG.mass.site", "log10.Zostera.proxy.production.rate.site", "log10.crustacean.mass.per.g.plant.site", 
  "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "log10.mean.fetch", "pop.density.2015",
  "log10.periphyton.mass.per.g.zostera.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Construct random forest: EELGRASS ABOVE-GROUND BIOMASS (with geographic variables)
RF.ZAG.1 = randomForest(log10.Zostera.AG.mass ~ 
                                  Ocean
                                + Coast
                                + basin
                                + Latitude
                                + Longitude
                                + Temperature.C
                                + Salinity.ppt
                                # + Day.length.hours
                                + log10.mean.fetch
                                + pop.density.2015
                                + Depth.Categorical
                                # + Depth.m
                                # + GenotypicRichness
                                + AllelicRichness
                                # + Inbreeding
                                # + log10.Zostera.AG.mass
                                # + log10.macrophytes.total.AG.mass.core
                                # + log10.Zostera.shoots.core
                                # + Zostera.sheath.width
                                # + Zostera.sheath.length
                                # + Zostera.longest.leaf.length
                                # + log10.pct.cover.macroalgae
                                # + pct.cover.seagrass
                                + log10.periphyton.mass.per.g.zostera  
                                + log10.mesograzer.mass.per.g.plant
                                + log10.crustacean.mass.per.g.plant
                                + log10.amphipod.mass.per.g.plant
                                + log10.gammarid.mass.per.g.plant
                                + amphipod.survival.24hr
                                + log10.Leaf.PercN
                                # + log10.GenotypicRichness
                                ,
                                data = ZEN_2014_master_data,
                                na.action = na.roughfix,
                                ntrees = 500,
                                importance = T
)

# Examine summary output
RF.ZAG.1 # % Var explained: 65.23

# Plot variable importance
varImpPlot(RF.ZAG.1) 
# RESULTS: Temperature, then (latitude, fetch, salinity)


# Construct random forest: EELGRASS ABOVE-GROUND BIOMASS (omitting geographic variables)
RF.ZAG.2 = randomForest(log10.Zostera.AG.mass ~ 
                        #   Ocean
                        # + Coast
                        # + basin
                        # + Latitude
                        # + Longitude
                        + Temperature.C
                        + Salinity.ppt
                        # + Day.length.hours
                        + log10.mean.fetch
                        + pop.density.2015
                        + Depth.Categorical
                        # + Depth.m
                        # + GenotypicRichness
                        + AllelicRichness
                        # + Inbreeding
                        # + log10.Zostera.AG.mass
                        # + log10.macrophytes.total.AG.mass.core
                        # + log10.Zostera.shoots.core
                        # + Zostera.sheath.width
                        # + Zostera.sheath.length
                        # + Zostera.longest.leaf.length
                        # + log10.pct.cover.macroalgae
                        # + pct.cover.seagrass
                        + log10.periphyton.mass.per.g.zostera  
                        + log10.mesograzer.mass.per.g.plant
                        + log10.crustacean.mass.per.g.plant
                        + log10.amphipod.mass.per.g.plant
                        + log10.gammarid.mass.per.g.plant
                        + amphipod.survival.24hr
                        + log10.Leaf.PercN
                        # + log10.GenotypicRichness
                        ,
                        data = ZEN_2014_master_data,
                        na.action = na.roughfix,
                        ntrees = 500,
                        importance = T
)

# Examine summary output
RF.ZAG.2 # % Var explained: 65.23

# Plot variable importance
varImpPlot(RF.ZAG.2) 
# RESULTS: Salinity, fetch, (periphyton, allelic richness)


###################################################################################
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: PLOT LEVEL (HIERARCHICAL)              #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable.
ZAG.plot.1 <- lme(log10.Zostera.AG.mass.imputed ~ 
                            Latitude
                          + Temperature.C
                          + Salinity.ppt
                          # + Day.length.hours
                          + log10.mean.fetch
                          # + Depth.Categorical
                          + pop.density.2015
                          # + log10.Zostera.AG.mass.imputed
                          # + log10.Zostera.shoots.core.imputed
                          # + log10.Zostera.sheath.length
                          # + log10.Zostera.longest.leaf.length
                          + AllelicRichness
                          + log10.Leaf.PercN.imputed
                          + log10.periphyton.mass.per.g.zostera.imputed
                          + log10.mesograzer.mass.per.g.plant.imputed
                          + log10.crustacean.mass.per.g.plant.imputed
                          + log10.gastropod.mass.per.g.plant.imputed
                          + grazer.richness.site
                          # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.1)
# N = 980
# AIC = -26.96025
sem.coefs(list(ZAG.plot.1), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                         response                                   predictor     estimate  std.error p.value
# 1  log10.Zostera.AG.mass.imputed                                    Latitude -0.311974011 0.16821339  0.0714
# 2  log10.Zostera.AG.mass.imputed                               Temperature.C  0.171451349 0.11638829  0.1490


# Hierarchical model (plot level): "COMPETITION" 
ZAG.plot.2 <- lme(log10.Zostera.AG.mass.imputed ~ 
                  #   Latitude
                  # + Temperature.C
                  # + Salinity.ppt
                  # # + Day.length.hours
                  # + log10.mean.fetch
                  # # + Depth.Categorical
                  # + pop.density.2015
                  # # + log10.Zostera.AG.mass.imputed
                  # # + log10.Zostera.shoots.core.imputed
                  # # + log10.Zostera.sheath.length
                  # # + log10.Zostera.longest.leaf.length
                  # + AllelicRichness
                  # + log10.Leaf.PercN.imputed
                  + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + grazer.richness.site
                  # + squid.survival.24hr.imputed
                  # + amphipod.survival.24hr # many missing value because only done at some sites
                  , random = ~1 | Coast/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.2)
# N = 980
# AIC = -108.0766
sem.coefs(list(ZAG.plot.2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                        response                                   predictor   estimate  std.error p.value
# 1 log10.Zostera.AG.mass.imputed log10.periphyton.mass.per.g.zostera.imputed 0.05051419 0.04601804  0.2726


# Hierarchical model (plot level): PRODUCTIVITY
ZAG.plot.3 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    #   Latitude
                    # + Temperature.C
                    # + Salinity.ppt
                    # # + Day.length.hours
                    # + log10.mean.fetch
                    # # + Depth.Categorical
                    # + pop.density.2015
                    # # + log10.Zostera.AG.mass.imputed
                    # # + log10.Zostera.shoots.core.imputed
                    # # + log10.Zostera.sheath.length
                    # # + log10.Zostera.longest.leaf.length
                  # + AllelicRichness
                  + log10.Leaf.PercN.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + grazer.richness.site
                  # + squid.survival.24hr.imputed
                  # + amphipod.survival.24hr # many missing value because only done at some sites
                  , random = ~1 | Coast/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.3)
# N = 980
# AIC = -111.1597
sem.coefs(list(ZAG.plot.3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                        response                predictor    estimate  std.error p.value
# 1 log10.Zostera.AG.mass.imputed log10.Leaf.PercN.imputed -0.03995435 0.03624634  0.2706


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
ZAG.plot.4 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    #   Latitude
                    + Temperature.C
                    + Salinity.ppt
                    + Day.length.hours
                    + log10.mean.fetch
                    # # + Depth.Categorical
                    # + pop.density.2015
                    # # + log10.Zostera.AG.mass.imputed
                    # # + log10.Zostera.shoots.core.imputed
                    # # + log10.Zostera.sheath.length
                    # # + log10.Zostera.longest.leaf.length
                  # + AllelicRichness
                  # + log10.Leaf.PercN.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + grazer.richness.site
                  # + squid.survival.24hr.imputed
                  # + amphipod.survival.24hr # many missing value because only done at some sites
                  , random = ~1 | Coast/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.4)
# N = 980
# AIC = -90.11734
sem.coefs(list(ZAG.plot.4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                        response        predictor    estimate std.error p.value
# 1 log10.Zostera.AG.mass.imputed    Temperature.C  0.26534373 0.1082524  0.0186
# 4 log10.Zostera.AG.mass.imputed log10.mean.fetch -0.11294057 0.1062914  0.2942
# 2 log10.Zostera.AG.mass.imputed     Salinity.ppt -0.07617149 0.1100448  0.4927
# 3 log10.Zostera.AG.mass.imputed Day.length.hours -0.04534299 0.1262972  0.7214


# Hierarchical model (plot level): ZOSTERA GENETICS
ZAG.plot.5 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    #   Latitude
                  #   + Temperature.C
                  # + Salinity.ppt
                  # + Day.length.hours
                  # + log10.mean.fetch
                  # # + Depth.Categorical
                  # + pop.density.2015
                  # # + log10.Zostera.AG.mass.imputed
                  # # + log10.Zostera.shoots.core.imputed
                  # # + log10.Zostera.sheath.length
                  # # + log10.Zostera.longest.leaf.length
                  + AllelicRichness
                  # + log10.Leaf.PercN.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + grazer.richness.site
                  # + squid.survival.24hr.imputed
                  # + amphipod.survival.24hr # many missing value because only done at some sites
                  , random = ~1 | Coast/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.5)
# N = 980
# AIC = -109.2214
sem.coefs(list(ZAG.plot.5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                        response       predictor  estimate std.error p.value
# 1 log10.Zostera.AG.mass.imputed AllelicRichness 0.1423411 0.1021208  0.1704


# Hierarchical model (plot level): TOP-DOWN CONTROL
ZAG.plot.6 <- lme(log10.Zostera.AG.mass.imputed ~ 
                    #   Latitude
                    #   + Temperature.C
                    # + Salinity.ppt
                    # + Day.length.hours
                    # + log10.mean.fetch
                    # # + Depth.Categorical
                    # + pop.density.2015
                    # # + log10.Zostera.AG.mass.imputed
                    # # + log10.Zostera.shoots.core.imputed
                    # # + log10.Zostera.sheath.length
                    # # + log10.Zostera.longest.leaf.length
                  # + AllelicRichness
                  # + log10.Leaf.PercN.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  + log10.mesograzer.mass.per.g.plant.imputed
                  + log10.crustacean.mass.per.g.plant.imputed
                  + log10.gastropod.mass.per.g.plant.imputed
                  + grazer.richness.site
                  # + squid.survival.24hr.imputed
                  # + amphipod.survival.24hr # many missing value because only done at some sites
                  , random = ~1 | Coast/Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed)
summary(ZAG.plot.6)
# N = 980
# AIC = -84.03823
sem.coefs(list(ZAG.plot.6), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                        response                                 predictor    estimate  std.error p.value
# 3 log10.Zostera.AG.mass.imputed  log10.gastropod.mass.per.g.plant.imputed -0.04810919 0.04588568  0.2947
# 4 log10.Zostera.AG.mass.imputed                      grazer.richness.site  0.10022777 0.10100778  0.3265
# 2 log10.Zostera.AG.mass.imputed log10.crustacean.mass.per.g.plant.imputed  0.03075345 0.04559190  0.5001
# 1 log10.Zostera.AG.mass.imputed log10.mesograzer.mass.per.g.plant.imputed -0.01812772 0.05582886  0.7455 

# Best model based on AIC is "productivity" based on leaf % N, but the single predictor is non-significant. 
# It appears that Zostera biomass is too variable on a plot scale to be predicted by the coarse 
# predictor variables we have. Interestingly, Zostera site-sale productivity is much better correlated
# with pedictors (see below). 


###################################################################################
# MODELS OF EELGRASS ABOVE-GROUND BIOMASS: SITE LEVEL (SIMPLE)                    #
###################################################################################

# Linear model (site level): PRODUCTIVITY
ZAG.lm.1 <- lm(log10.Zostera.AG.mass.site ~ 
                         log10.periphyton.mass.per.g.zostera.site 
                       + Leaf.PercN.site
                       , data = ZEN_2014_site_means_49)
summary(ZAG.lm.1)
#                                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                               2.35851    0.19765  11.933 1.11e-15 ***
#   log10.periphyton.mass.per.g.zostera.site  0.20394    0.06141   3.321  0.00176 ** 
#   Leaf.PercN.site                          -0.06438    0.08428  -0.764  0.44884    
# 
# Residual standard error: 0.3098 on 46 degrees of freedom
# Multiple R-squared:  0.1936,	Adjusted R-squared:  0.1586 
# F-statistic: 5.524 on 2 and 46 DF,  p-value: 0.00708
AIC(ZAG.lm.1) # 29.10888


# Linear model (site level): ABIOTIC ENVIRONMENT
ZAG.lm.2 <- lm(log10.Zostera.AG.mass.site ~ 
                         Temperature.C 
                       + Salinity.ppt
                       + Day.length.hours 
                       + log10.mean.fetch
                       , data = ZEN_2014_site_means_49)
summary(ZAG.lm.2)
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       2.9406785  0.5083713   5.785 6.99e-07 ***
# Temperature.C     0.0066434  0.0093762   0.709  0.48234    
# Salinity.ppt      0.0005114  0.0064431   0.079  0.93710    
# Day.length.hours -0.0682093  0.0220202  -3.098  0.00339 ** 
# log10.mean.fetch -0.0193782  0.1069347  -0.181  0.85703    
# 
# Residual standard error: 0.2992 on 44 degrees of freedom
# Multiple R-squared:  0.2806,	Adjusted R-squared:  0.2151 
# F-statistic: 4.289 on 4 and 44 DF,  p-value: 0.005125
AIC(ZAG.lm.2) # 27.52131


# Linear model (site level): GEOGRAPHY
ZAG.lm.3 <- lm(log10.Zostera.AG.mass.site ~ 
                         Latitude
                       + Ocean 
                       + Coast
                       , data = ZEN_2014_site_means_49)
summary(ZAG.lm.3)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         2.406928   0.256637   9.379 4.65e-12 ***
# Latitude           -0.012797   0.004667  -2.742   0.0088 ** 
# OceanPacific        0.306420   0.137344   2.231   0.0308 *  
# CoastEast.Pacific   0.001234   0.118642   0.010   0.9918    
# CoastWest.Atlantic  0.260591   0.110185   2.365   0.0225 *  
# CoastWest.Pacific         NA         NA      NA       NA    
# 
# Residual standard error: 0.2418 on 44 degrees of freedom
# Multiple R-squared:  0.5299,	Adjusted R-squared:  0.4871 
# F-statistic:  12.4 on 4 and 44 DF,  p-value: 7.787e-07
AIC(ZAG.lm.3) # 6.674112


# Linear model (site level): ZOSTERA GENETICS
ZAG.lm.4 <- lm(log10.Zostera.AG.mass.site ~ 
                         AllelicRichness
                       , data = ZEN_2014_site_means_49)
summary(ZAG.lm.4)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.55285    0.16146   9.618 1.11e-12 ***
# AllelicRichness  0.09967    0.03260   3.058  0.00368 ** 
# 
# Residual standard error: 0.3117 on 47 degrees of freedom
# Multiple R-squared:  0.1659,	Adjusted R-squared:  0.1482 
# F-statistic: 9.348 on 1 and 47 DF,  p-value: 0.003676
AIC(ZAG.lm.4) # 28.76661


# SUMMARY OF RESULTS: SITE level
# Geography is the story for site-level Zostera biomass. Geography model is far better than 
# any other and suggests that Zostera biomass declines with latitude and is generally
# greater in West Atlantic than on other coasts. 


###################################################################################
# MODELS OF EELGRASS PRODUCTIVITY: RANDOM FORESTS                                 #
###################################################################################

# First, create an estimate of eelgrass biomass producivity rate. Here we use sheath length
# as a proxy for leaf elongation rate, based on Jen Ruesink's analysis of ZEN 1 (2011) data, 
# and multiply that estimate of leaf elongation rate by shoot density to get rough proxy 
# for biomass production rate.

ZEN_2014_master_data$Zostera.proxy.production.rate <- ZEN_2014_master_data$Zostera.sheath.length * ZEN_2014_master_data$Zostera.shoots.per.m2.core
ZEN_2014_master_data$log10.Zostera.proxy.production.rate <- log10(ZEN_2014_master_data$Zostera.proxy.production.rate)
hist(ZEN_2014_master_data$Zostera.proxy.production.rate)
hist(ZEN_2014_master_data$log10.Zostera.proxy.production.rate) # Nice. Use log10-transformed data

ZEN_2014_master_data_imputed$log10.Zostera.proxy.production.rate <- ZEN_2014_master_data$log10.Zostera.proxy.production.rate[match(ZEN_2014_master_data_imputed$Unique.ID, ZEN_2014_master_data$Unique.ID)]
ZEN_2014_master_data_49_imputed$log10.Zostera.proxy.production.rate <- ZEN_2014_master_data$log10.Zostera.proxy.production.rate[match(ZEN_2014_master_data_49_imputed$Unique.ID, ZEN_2014_master_data$Unique.ID)]

x <- ddply(ZEN_2014_master_data, c("Site"), summarize, log10.Zostera.proxy.production.rate.site = mean(log10.Zostera.proxy.production.rate, na.rm = T))
nrow(x)

ZEN_2014_site_means$log10.Zostera.proxy.production.rate.site <- x$log10.Zostera.proxy.production.rate.site
  [match(ZEN_2014_site_means$Site, x$Site)]
names(ZEN_2014_site_means)

ZEN_2014_site_means_49$log10.Zostera.proxy.production.rate.site <- x$log10.Zostera.proxy.production.rate.site[match(ZEN_2014_site_means_49$Site, x$Site)]
names(ZEN_2014_site_means)


# Explore graphically

pairs.panels(ZEN_2014_site_means[,c("log10.Zostera.proxy.production.rate.site", "log10.crustacean.mass.per.g.plant.site", "log10.Zostera.AG.mass.site",
  "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "log10.mean.fetch", "pop.density.2015",
  "log10.periphyton.mass.per.g.zostera.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Construct random forest: ZOSTERA PRODUCTION RATE (with geographic variables)
RF.Zproduction.1 = randomForest(log10.Zostera.proxy.production.rate ~ 
    Ocean
  + Coast
  + basin
  + Latitude
  + Longitude
  + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  # + log10.Zostera.AG.mass
  # + log10.macrophytes.total.AG.mass.core
  # + log10.Zostera.shoots.core
  # + Zostera.sheath.width
  # + Zostera.sheath.length
  # + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera  
  + log10.mesograzer.mass.per.g.plant
  + log10.crustacean.mass.per.g.plant
  + log10.amphipod.mass.per.g.plant
  + log10.gammarid.mass.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.Zproduction.1 # % Var explained: 64.41

# Plot variable importance
varImpPlot(RF.Zproduction.1) 
# RESULTS: temperature stands out


# Construct random forest: ZOSTERA PRODUCTION RATE (omitting geographic variables)
RF.Zproduction.2 = randomForest(log10.Zostera.proxy.production.rate ~ 
                               #   Ocean
                               # + Coast
                               # + basin
                               # + Latitude
                               # + Longitude
                               + Temperature.C
                               + Salinity.ppt
                               # + Day.length.hours
                               + log10.mean.fetch
                               + pop.density.2015
                               + Depth.Categorical
                               # + Depth.m
                               # + Perc.Silt
                               # + Perc.Sand
                               # + Perc.Gravel
                               # + GenotypicRichness
                               + AllelicRichness
                               # + Inbreeding
                               # + log10.Zostera.AG.mass
                               # + log10.macrophytes.total.AG.mass.core
                               # + log10.Zostera.shoots.core
                               # + Zostera.sheath.width
                               # + Zostera.sheath.length
                               # + Zostera.longest.leaf.length
                               # + log10.pct.cover.macroalgae
                               # + pct.cover.seagrass
                               + log10.periphyton.mass.per.g.zostera  
                               + log10.mesograzer.mass.per.g.plant
                               + log10.crustacean.mass.per.g.plant
                               + log10.amphipod.mass.per.g.plant
                               + log10.gammarid.mass.per.g.plant
                               #   + log10.gastropod.abund.per.g.plant
                               #   + log10.grazer.richness.plot
                               #   + log10.grazer.richness.site
                               + amphipod.survival.24hr
                               + log10.Leaf.PercN
                               # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.Zproduction.2 # % Var explained: 63.34

# Plot variable importance
varImpPlot(RF.Zproduction.2) 
# RESULT: temperature. Inteesting - geography has no discrenable effect on Zostera 
# production, beyond the effect of temperature. 


###################################################################################
# MODELS OF EELGRASS PRODUCTIVITY: PLOT LEVEL (HIERARCHICAL)                      #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable.
Zproduction.plot.1 <- lme(log10.Zostera.proxy.production.rate ~ 
                           Latitude
                         + Temperature.C
                         + Salinity.ppt
                         # + Day.length.hours
                         + log10.mean.fetch
                         # + Depth.Categorical
                         + pop.density.2015
                         # + log10.Zostera.AG.mass.imputed
                         # + log10.Zostera.shoots.core.imputed
                         # + log10.Zostera.sheath.length
                         # + log10.Zostera.longest.leaf.length
                         + AllelicRichness
                         + log10.Leaf.PercN.imputed
                         + log10.periphyton.mass.per.g.zostera.imputed
                         + log10.mesograzer.mass.per.g.plant.imputed
                         + log10.crustacean.mass.per.g.plant.imputed
                         + log10.gastropod.mass.per.g.plant.imputed
                         + grazer.richness.site
                         # + squid.survival.24hr.imputed
                         # + amphipod.survival.24hr # many missing value because only done at some sites
                         , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.1)
# N = 965
# AIC = 110.5114
sem.coefs(list(Zproduction.plot.1), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                               response                                   predictor     estimate  std.error p.value
# 8  log10.Zostera.proxy.production.rate log10.periphyton.mass.per.g.zostera.imputed  0.169684614 0.04832652  0.0005
# 1  log10.Zostera.proxy.production.rate                                    Latitude -0.456668676 0.17651770  0.0136
# 3  log10.Zostera.proxy.production.rate                                Salinity.ppt -0.371676778 0.14624499  0.0152
# 5  log10.Zostera.proxy.production.rate                            pop.density.2015 -0.317052904 0.12672012  0.0168
# 9  log10.Zostera.proxy.production.rate   log10.mesograzer.mass.per.g.plant.imputed -0.133391755 0.05711571  0.0197
# 2  log10.Zostera.proxy.production.rate                               Temperature.C  0.146481611 0.13497815  0.2847


# Hierarchical model (plot level): "COMPETITION" (but relation turns out to be positive ...)
Zproduction.plot.2 <- lme(log10.Zostera.proxy.production.rate ~ 
                          #   Latitude
                          # + Temperature.C
                          # + Salinity.ppt
                          # # + Day.length.hours
                          # + log10.mean.fetch
                          # # + Depth.Categorical
                          # + pop.density.2015
                          # # + log10.Zostera.AG.mass.imputed
                          # # + log10.Zostera.shoots.core.imputed
                          # # + log10.Zostera.sheath.length
                          # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.2)
# N = 965
# AIC = 40.20471
sem.coefs(list(Zproduction.plot.2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                              response                                   predictor  estimate std.error p.value
# 1 log10.Zostera.proxy.production.rate log10.periphyton.mass.per.g.zostera.imputed 0.1825507 0.0478227   1e-04


# Hierarchical model (plot level): PRODUCTIVITY
Zproduction.plot.3 <- lme(log10.Zostera.proxy.production.rate ~ 
                          #   Latitude
                          # + Temperature.C
                          # + Salinity.ppt
                          # # + Day.length.hours
                          # + log10.mean.fetch
                          # # + Depth.Categorical
                          # + pop.density.2015
                          # # + log10.Zostera.AG.mass.imputed
                          # # + log10.Zostera.shoots.core.imputed
                          # # + log10.Zostera.sheath.length
                          # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.3)
# N = 965
# AIC = 51.43112
sem.coefs(list(Zproduction.plot.3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                              response                predictor   estimate  std.error p.value
# 1 log10.Zostera.proxy.production.rate log10.Leaf.PercN.imputed 0.01259627 0.03793142  0.7399


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
Zproduction.plot.4 <- lme(log10.Zostera.proxy.production.rate ~ 
                            #   Latitude
                            + Temperature.C
                          + Salinity.ppt
                          + Day.length.hours
                          + log10.mean.fetch
                          # + Depth.Categorical
                          # + pop.density.2015
                          # + log10.Zostera.AG.mass.imputed
                          # + log10.Zostera.shoots.core.imputed
                          # + log10.Zostera.sheath.length
                          # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.Leaf.PercN.imputed
                          # # + log10.mesograzer.mass.per.g.plant.imputed
                          # # + log10.crustacean.mass.per.g.plant.imputed
                          # # + log10.gastropod.mass.per.g.plant.imputed
                          # # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.4)
# N = 965
# AIC = 69.05987
sem.coefs(list(Zproduction.plot.4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                              response        predictor   estimate std.error p.value
# 3 log10.Zostera.proxy.production.rate log10.mean.fetch -0.1819258 0.1306524  0.1711
# 2 log10.Zostera.proxy.production.rate     Salinity.ppt -0.1816213 0.1353363  0.1868
# 1 log10.Zostera.proxy.production.rate    Temperature.C  0.1226584 0.1163773  0.2979


# Hierarchical model (plot level): ZOSTERA GENETICS
Zproduction.plot.5 <- lme(log10.Zostera.proxy.production.rate ~ 
                          # #   Latitude
                          # # + Temperature.C
                          # # + Salinity.ppt
                          # # # + Day.length.hours
                          # # + log10.mean.fetch
                          # # # + Depth.Categorical
                          # + pop.density.2015
                          # # + log10.Zostera.AG.mass.imputed
                          # # + log10.Zostera.shoots.core.imputed
                          # # + log10.Zostera.sheath.length
                          # # + log10.Zostera.longest.leaf.length
                          + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.5)
# N = 965
# AIC = 53.93175
sem.coefs(list(Zproduction.plot.5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                              response       predictor     estimate std.error p.value
# 1 log10.Zostera.proxy.production.rate AllelicRichness -0.001388548 0.1198613  0.9908


# Hierarchical model (plot level): TOP-DOWN CONTROL
Zproduction.plot.6 <- lme(log10.Zostera.proxy.production.rate ~ 
                          #   Latitude
                          # + Temperature.C
                          # + Salinity.ppt
                          # # + Day.length.hours
                          # + log10.mean.fetch
                          # # + Depth.Categorical
                          # + pop.density.2015
                          # # + log10.Zostera.AG.mass.imputed
                          # # + log10.Zostera.shoots.core.imputed
                          # # + log10.Zostera.sheath.length
                          # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          + log10.mesograzer.mass.per.g.plant.imputed
                          + log10.crustacean.mass.per.g.plant.imputed
                          + log10.gastropod.mass.per.g.plant.imputed
                          + grazer.richness.site
                          # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zproduction.plot.6)
# N = 965
# AIC = 66.52063
sem.coefs(list(Zproduction.plot.6), ZEN_2014_master_data_49_imputed, standardize = "scale")
# 1 log10.Zostera.proxy.production.rate log10.mesograzer.mass.per.g.plant.imputed -0.137960425 0.05729216  0.0162
# 4 log10.Zostera.proxy.production.rate                      grazer.richness.site  0.193770975 0.11803229  0.1078
# 3 log10.Zostera.proxy.production.rate  log10.gastropod.mass.per.g.plant.imputed  0.025695400 0.04712337  0.5857
# 2 log10.Zostera.proxy.production.rate log10.crustacean.mass.per.g.plant.imputed  0.003211136 0.04692921  0.9455


# SUMMARY OF RESULTS: PLOT level
# By far the best model is the "competition" model including only periphyton biomass. 
# Since the relationships with periphyton biomass is positive, however, this
# is better thought of as a productivity model that integrates factors promoting plant
# growth (e.g., light, other nutrients, perhaps flow, etc.). What's good for periphyton 
# is good for eelgrass.


###################################################################################
# MODELS OF EELGRASS PRODUCTIVITY: SITE LEVEL (SIMPLE)                            #
###################################################################################

# Linear model (site level): PRODUCTIVITY
Zproduction.lm.1 <- lm(log10.Zostera.proxy.production.rate.site ~ 
                         log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(Zproduction.lm.1)
#                                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                               3.42376    0.21510  15.917   <2e-16 ***
# log10.periphyton.mass.per.g.zostera.site  0.12984    0.06684   1.943   0.0582 .  
# Leaf.PercN.site                           0.08343    0.09172   0.910   0.3678    
# 
# Residual standard error: 0.3371 on 46 degrees of freedom
# Multiple R-squared:  0.1061,	Adjusted R-squared:  0.06721 
# F-statistic: 2.729 on 2 and 46 DF,  p-value: 0.07584
AIC(Zproduction.lm.1) # 37.39819


# Linear model (site level): ABIOTIC ENVIRONMENT
Zproduction.lm.2 <- lm(log10.Zostera.proxy.production.rate.site ~ 
                         Temperature.C 
                       + Salinity.ppt
                       + Day.length.hours 
                       + log10.mean.fetch
                       , data = ZEN_2014_site_means_49)
summary(Zproduction.lm.2)
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       4.6466078  0.5742101   8.092 2.98e-10 ***
# Temperature.C    -0.0007009  0.0105905  -0.066   0.9475    
# Salinity.ppt     -0.0071194  0.0072775  -0.978   0.3333    
# Day.length.hours -0.0615322  0.0248721  -2.474   0.0173 *  
# log10.mean.fetch -0.0603455  0.1207838  -0.500   0.6198    
# 
# Residual standard error: 0.3379 on 44 degrees of freedom
# Multiple R-squared:  0.1408,	Adjusted R-squared:  0.06271 
# F-statistic: 1.803 on 4 and 44 DF,  p-value: 0.1454
AIC(Zproduction.lm.2) # 39.45608


# Linear model (site level): GEOGRAPHY
Zproduction.lm.3 <- lm(log10.Zostera.proxy.production.rate.site ~ 
                         Latitude
                       + Ocean 
                       + Coast
                       , data = ZEN_2014_site_means_49)
summary(Zproduction.lm.3)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         3.491516   0.328043  10.643 9.43e-14 ***
#   Latitude           -0.004765   0.005966  -0.799  0.42879    
# OceanPacific        0.186335   0.175558   1.061  0.29431    
# CoastEast.Pacific   0.054907   0.151652   0.362  0.71904    
# CoastWest.Atlantic  0.389707   0.140842   2.767  0.00824 ** 
#   CoastWest.Pacific         NA         NA      NA       NA    
# 
# Residual standard error: 0.3091 on 44 degrees of freedom
# Multiple R-squared:  0.281,	Adjusted R-squared:  0.2156 
# F-statistic: 4.298 on 4 and 44 DF,  p-value: 0.005068
AIC(Zproduction.lm.3) # 30.73131


# Linear model (site level): ZOSTERA GENETICS
Zproduction.lm.4 <- lm(log10.Zostera.proxy.production.rate.site ~ 
                         AllelicRichness
                       , data = ZEN_2014_site_means_49)
summary(Zproduction.lm.4)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.24733    0.17966  18.075   <2e-16 ***
#   AllelicRichness  0.04612    0.03627   1.272     0.21    
# 
# Residual standard error: 0.3468 on 47 degrees of freedom
# Multiple R-squared:  0.03326,	Adjusted R-squared:  0.01269 
# F-statistic: 1.617 on 1 and 47 DF,  p-value: 0.2098
AIC(Zproduction.lm.4) # 39.2355


# SUMMARY OF RESULTS: SITE level
# Looks like it's all geography. West Atlantic has more productive 
# eelgrass than anywhere else. Who knew?  Reasons remain obscure since there 
# is little or no correlation with leaf N or the usual abiotic factors. 

# Zostera proxy production x latitude 
Zprodn.lat = ggplot(ZEN_2014_site_means, aes(x = Latitude, y = log10.Zostera.proxy.production.rate.site, col = Ocean)) +
  geom_point(size = 5) +
  geom_text(aes(label = unique(ZEN_2014_site_means$Site)), hjust = -0.25, vjust = 0) +
  scale_color_manual(values = c("blue",  "forestgreen")) +
  xlab("\n  Latitude") +  
  ylab("log Zostera production (proxy) \n") +  
  scale_x_continuous(limits = c(30,72))+
  #   scale_y_continuous(limits=c(0,3))+
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),axis.ticks = element_line(size = 1.5))
Zprodn.lat + geom_smooth(method = lm, fullrange = F, se = F, lwd = 1.5, col="black", na.rm=T)


###################################################################################
# MODELS OF EELGRASS SHOOT DENSITY: PLOT LEVEL (HIERARCHICAL)                     #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable.
Zshoot.plot.1 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            Latitude
                          + Temperature.C
                          + Salinity.ppt
                          # + Day.length.hours
                          + log10.mean.fetch
                          # + Depth.Categorical
                          + pop.density.2015
                          # + log10.Zostera.AG.mass.imputed
                          # + log10.Zostera.shoots.core.imputed
                          # + log10.Zostera.sheath.length
                          # + log10.Zostera.longest.leaf.length
                          + AllelicRichness
                          + log10.Leaf.PercN.imputed
                          + log10.periphyton.mass.per.g.zostera.imputed
                          + log10.mesograzer.mass.per.g.plant.imputed
                          + log10.crustacean.mass.per.g.plant.imputed
                          + log10.gastropod.mass.per.g.plant.imputed
                          + grazer.richness.site
                          # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.1)
# N = 980
# AIC = -110.2985
sem.coefs(list(Zshoot.plot.1), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                             response                                   predictor     estimate  std.error p.value
# 3  log10.Zostera.shoots.core.imputed                                Salinity.ppt -0.348540351 0.16662691  0.0432
# 9  log10.Zostera.shoots.core.imputed   log10.mesograzer.mass.per.g.plant.imputed -0.083074614 0.04837212  0.0862
# 8  log10.Zostera.shoots.core.imputed log10.periphyton.mass.per.g.zostera.imputed  0.059652011 0.04147171  0.1507


# Hierarchical model (plot level): "COMPETITION" (but relation turns out to be positive ...)
Zshoot.plot.2 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            #   Latitude
                            # + Temperature.C
                            # + Salinity.ppt
                            # # + Day.length.hours
                            # + log10.mean.fetch
                            # # + Depth.Categorical
                            # + pop.density.2015
                            # # + log10.Zostera.AG.mass.imputed
                            # # + log10.Zostera.shoots.core.imputed
                            # # + log10.Zostera.sheath.length
                            # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.2)
# N = 980
# AIC = -187.1107
sem.coefs(list(Zshoot.plot.2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                            response                                   predictor   estimate  std.error p.value
# 1 log10.Zostera.shoots.core.imputed log10.periphyton.mass.per.g.zostera.imputed 0.05918096 0.04111075  0.1503

# Hierarchical model (plot level): PRODUCTIVITY
Zshoot.plot.3 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            #   Latitude
                            # + Temperature.C
                            # + Salinity.ppt
                            # # + Day.length.hours
                            # + log10.mean.fetch
                            # # + Depth.Categorical
                            # + pop.density.2015
                            # # + log10.Zostera.AG.mass.imputed
                            # # + log10.Zostera.shoots.core.imputed
                            # # + log10.Zostera.sheath.length
                            # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.3)
# N = 980
# AIC = -188.1395
sem.coefs(list(Zshoot.plot.3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                            response                predictor   estimate  std.error p.value
# 1 log10.Zostera.shoots.core.imputed log10.Leaf.PercN.imputed 0.00491861 0.03195534  0.8777


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
Zshoot.plot.4 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            #   Latitude
                            + Temperature.C
                          + Salinity.ppt
                          + Day.length.hours
                          + log10.mean.fetch
                          # + Depth.Categorical
                          # + pop.density.2015
                          # + log10.Zostera.AG.mass.imputed
                          # + log10.Zostera.shoots.core.imputed
                          # + log10.Zostera.sheath.length
                          # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.Leaf.PercN.imputed
                          # # + log10.mesograzer.mass.per.g.plant.imputed
                          # # + log10.crustacean.mass.per.g.plant.imputed
                          # # + log10.gastropod.mass.per.g.plant.imputed
                          # # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.4)
# N = 980
# AIC = -166.7987
sem.coefs(list(Zshoot.plot.4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                            response        predictor   estimate std.error p.value
# 2 log10.Zostera.shoots.core.imputed     Salinity.ppt -0.2621729 0.1437958  0.0756
# 1 log10.Zostera.shoots.core.imputed    Temperature.C  0.2179574 0.1408351  0.1294
# 4 log10.Zostera.shoots.core.imputed log10.mean.fetch -0.1939976 0.1385052  0.1688
# 3 log10.Zostera.shoots.core.imputed Day.length.hours  0.1720676 0.1623290  0.2954


# Hierarchical model (plot level): ZOSTERA GENETICS
Zshoot.plot.5 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            # #   Latitude
                            # # + Temperature.C
                            # # + Salinity.ppt
                            # # # + Day.length.hours
                            # # + log10.mean.fetch
                            # # # + Depth.Categorical
                            # + pop.density.2015
                            # # + log10.Zostera.AG.mass.imputed
                            # # + log10.Zostera.shoots.core.imputed
                            # # + log10.Zostera.sheath.length
                            # # + log10.Zostera.longest.leaf.length
                          + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          # + log10.mesograzer.mass.per.g.plant.imputed
                          # + log10.crustacean.mass.per.g.plant.imputed
                          # + log10.gastropod.mass.per.g.plant.imputed
                          # + grazer.richness.site
                          # # + squid.survival.24hr.imputed
                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.5)
# N = 980
# AIC = -186.2219
sem.coefs(list(Zshoot.plot.5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                            response       predictor   estimate std.error p.value
# 1 log10.Zostera.shoots.core.imputed AllelicRichness 0.01488269 0.1285667  0.9084

# Hierarchical model (plot level): TOP-DOWN CONTROL
Zshoot.plot.6 <- lme(log10.Zostera.shoots.core.imputed ~ 
                            #   Latitude
                            # + Temperature.C
                            # + Salinity.ppt
                            # # + Day.length.hours
                            # + log10.mean.fetch
                            # # + Depth.Categorical
                            # + pop.density.2015
                            # # + log10.Zostera.AG.mass.imputed
                            # # + log10.Zostera.shoots.core.imputed
                            # # + log10.Zostera.sheath.length
                            # # + log10.Zostera.longest.leaf.length
                          # + AllelicRichness
                          # + log10.Leaf.PercN.imputed
                          # + log10.periphyton.mass.per.g.zostera.imputed
                          + log10.mesograzer.mass.per.g.plant.imputed
                          + log10.crustacean.mass.per.g.plant.imputed
                          + log10.gastropod.mass.per.g.plant.imputed
                          + grazer.richness.site
                          # + squid.survival.24hr.imputed
                          # + amphipod.survival.24hr # many missing value because only done at some sites
                          , random = ~1 | Coast/Site
                          , na.action = na.omit
                          , data = ZEN_2014_master_data_49_imputed)
summary(Zshoot.plot.6)
# N = 980
# AIC = -168.6221
sem.coefs(list(Zshoot.plot.6), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                            response                                 predictor    estimate  std.error p.value
# 1 log10.Zostera.shoots.core.imputed log10.mesograzer.mass.per.g.plant.imputed -0.08472008 0.04820793  0.0792
# 4 log10.Zostera.shoots.core.imputed                      grazer.richness.site  0.20400858 0.12623338  0.1132
# 3 log10.Zostera.shoots.core.imputed  log10.gastropod.mass.per.g.plant.imputed  0.01814601 0.03969156  0.6477
# 2 log10.Zostera.shoots.core.imputed log10.crustacean.mass.per.g.plant.imputed -0.01495811 0.03924397  0.7032


# SUMMARY OF RESULTS: PLOT level


###################################################################################
# MODELS OF EELGRASS SHOOT DENSITY: SITE LEVEL (SIMPLE)                           #
###################################################################################

# Linear model (site level): PRODUCTIVITY
Zshoot.lm.1 <- lm(log10.Zostera.shoots.core.site ~ 
                         log10.periphyton.mass.per.g.zostera.site 
                       + Leaf.PercN.site
                       , data = ZEN_2014_site_means_49)
summary(Zshoot.lm.1)
#                                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                               2.37476    0.24258   9.789 7.97e-13 ***
#   log10.periphyton.mass.per.g.zostera.site -0.14861    0.07538  -1.972   0.0547 .  
# Leaf.PercN.site                          -0.02632    0.10344  -0.254   0.8003    
# 
# Residual standard error: 0.3802 on 46 degrees of freedom
# Multiple R-squared:  0.08534,	Adjusted R-squared:  0.04557 
# F-statistic: 2.146 on 2 and 46 DF,  p-value: 0.1285
AIC(Zshoot.lm.1) # 49.18255


# Linear model (site level): ABIOTIC ENVIRONMENT
Zshoot.lm.2 <- lm(log10.Zostera.shoots.core.site ~ 
                         Temperature.C 
                       + Salinity.ppt
                       + Day.length.hours 
                       + log10.mean.fetch
                       , data = ZEN_2014_site_means_49)
summary(Zshoot.lm.2)
#                   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)       2.122420   0.653504   3.248  0.00223 **
# Temperature.C     0.017211   0.012053   1.428  0.16036   
# Salinity.ppt     -0.014010   0.008282  -1.691  0.09782 . 
# Day.length.hours  0.027323   0.028307   0.965  0.33970   
# log10.mean.fetch -0.111015   0.137463  -0.808  0.42367   
# 
# Residual standard error: 0.3846 on 44 degrees of freedom
# Multiple R-squared:  0.1047,	Adjusted R-squared:  0.02334 
# F-statistic: 1.287 on 4 and 44 DF,  p-value: 0.2898
AIC(Zshoot.lm.2) # 52.13272


# Linear model (site level): GEOGRAPHY
Zshoot.lm.3 <- lm(log10.Zostera.shoots.core.site ~ 
                         Latitude
                       + Ocean 
                       + Coast
                       , data = ZEN_2014_site_means_49)
summary(Zshoot.lm.3)
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         2.5421395  0.3736053   6.804 2.21e-08 ***
# Latitude           -0.0001238  0.0067946  -0.018   0.9855    
# OceanPacific       -0.4570895  0.1999415  -2.286   0.0271 *  
# CoastEast.Pacific   0.2692044  0.1727158   1.559   0.1262    
# CoastWest.Atlantic  0.1616185  0.1604041   1.008   0.3192    
# CoastWest.Pacific          NA         NA      NA       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3521 on 44 degrees of freedom
# Multiple R-squared:  0.2497,	Adjusted R-squared:  0.1815 
# F-statistic: 3.661 on 4 and 44 DF,  p-value: 0.01169
AIC(Zshoot.lm.3) # 43.47681


# Linear model (site level): ZOSTERA GENETICS
Zshoot.lm.4 <- lm(log10.Zostera.shoots.core.site ~ 
                         AllelicRichness
                       , data = ZEN_2014_site_means_49)
summary(Zshoot.lm.4)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       2.5175     0.2036   12.37   <2e-16 ***
# AllelicRichness  -0.0107     0.0411   -0.26    0.796    
# 
# Residual standard error: 0.393 on 47 degrees of freedom
# Multiple R-squared:  0.00144,	Adjusted R-squared:  -0.01981 
# F-statistic: 0.06777 on 1 and 47 DF,  p-value: 0.7958
AIC(Zshoot.lm.4) # 51.48302


# SUMMARY OF RESULTS: SITE level
# Shoot density is lower in Pacific. That's it. 



###################################################################################
# MODELS OF PERIPHYTON BIOMASS: RANDOM FORESTS                                    #
###################################################################################

# Explore graphically
pairs.panels(ZEN_2014_site_means[,c("log10.periphyton.mass.per.g.zostera.site", "log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
  "log10.Zostera.longest.leaf.length.cm.site", "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", 
  "log10.mean.fetch", "pop.density.2015", "log10.crustacean.mass.per.g.plant.site", "log10.gastropod.mass.per.g.plant.site", 
  "amphipod.survival.24hr.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Construct random forest: PERIPHYTON BIOMASS (with geographic variables)
RF.periphyton.1 = randomForest(log10.periphyton.mass.per.g.zostera ~
    Ocean
  + Coast
  + basin
  + Latitude
  + Longitude
  + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.proxy.production.rate
  + log10.Zostera.AG.mass
  + log10.macrophytes.total.AG.mass.core
  + log10.Zostera.shoots.core
  + Zostera.sheath.width
  + Zostera.sheath.length
  + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  # + log10.mesograzer.abund.per.g.plant
  + log10.crustacean.mass.per.g.plant
  + log10.amphipod.mass.per.g.plant
  + log10.gammarid.mass.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.periphyton.1 # % Var explained: 85.99

# Plot variable importance
varImpPlot(RF.periphyton.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULTS: Latitude and fetch highest, followed closely by (temperature, longitude, salinity) 


# Construct random forest: PERIPHYTON BIOMASS (omitting geographic variables)
RF.periphyton.2 = randomForest(log10.periphyton.mass.per.g.zostera ~
                               #   Ocean
                               # + Coast
                               # + basin
                               # + Latitude
                               # + Longitude
                               + Temperature.C
                               + Salinity.ppt
                               # + Day.length.hours
                               + log10.mean.fetch
                               + pop.density.2015
                               + Depth.Categorical
                               # + Depth.m
                               # + GenotypicRichness
                               + AllelicRichness
                               # + Inbreeding
                               + log10.Zostera.proxy.production.rate
                               + log10.Zostera.AG.mass
                               + log10.macrophytes.total.AG.mass.core
                               + log10.Zostera.shoots.core
                               + Zostera.sheath.width
                               + Zostera.sheath.length
                               + Zostera.longest.leaf.length
                               # + log10.pct.cover.macroalgae
                               # + pct.cover.seagrass
                               # + log10.mesograzer.abund.per.g.plant
                               + log10.crustacean.mass.per.g.plant
                               + log10.amphipod.mass.per.g.plant
                               + log10.gammarid.mass.per.g.plant
                               #   + log10.gastropod.abund.per.g.plant
                               #   + log10.grazer.richness.plot
                               #   + log10.grazer.richness.site
                               + amphipod.survival.24hr
                               + log10.Leaf.PercN
                               # + log10.GenotypicRichness
                               ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.periphyton.2 # % Var explained: 85.61

# Plot variable importance
varImpPlot(RF.periphyton.2) 
# RESULT: Strongest predictor is fetch, followed closely by (temperature, salinity, allelic richness ...)


###################################################################################
# MODELS OF PERIPHYTON BIOMASS: PLOT LEVEL (HIERARCHICAL)                         #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable.
# We remove predation (amphipod.survival.24hr) as predictor because many sites did not do  
# this so we lose plots, hence DF, and AICs can no longer be fairly compared across models. 
periphyton.plot.1 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                      Latitude
                    + Temperature.C
                    + Salinity.ppt
                    # + Day.length.hours
                    + log10.mean.fetch
                    # + Depth.Categorical
                    + pop.density.2015
                    + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    + log10.Leaf.PercN.imputed
                    + log10.mesograzer.mass.per.g.plant.imputed
                    + log10.crustacean.mass.per.g.plant.imputed
                    + log10.gastropod.mass.per.g.plant.imputed
                    + grazer.richness.site
                    # + squid.survival.24hr.imputed
                    # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.1)
# N = 980
# AIC = 538.5664
sem.coefs(list(periphyton.plot.1), ZEN_2014_master_data_imputed, standardize = "scale")
#                                       response                                 predictor     estimate  std.error p.value
# 9  log10.periphyton.mass.per.g.zostera.imputed         log10.Zostera.longest.leaf.length  0.161501068 0.04790209  0.0008
# 8  log10.periphyton.mass.per.g.zostera.imputed               log10.Zostera.sheath.length  0.116148781 0.04511524  0.0102
# 13 log10.periphyton.mass.per.g.zostera.imputed log10.crustacean.mass.per.g.plant.imputed  0.064155770 0.02980550  0.0316
# 1  log10.periphyton.mass.per.g.zostera.imputed                                  Latitude -0.323611074 0.16710516  0.0603
# 2  log10.periphyton.mass.per.g.zostera.imputed                             Temperature.C -0.197337143 0.12856952  0.1331


# Hierarchical model (plot level): SEAGRASS STRUCTURE
periphyton.plot.2 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                         #   Latitude
                         # + Temperature.C
                         # + Salinity.ppt
                         # # + Day.length.hours
                         # + log10.mean.fetch
                         # # + Depth.Categorical
                         # + pop.density.2015
                         + log10.Zostera.AG.mass.imputed
                         + log10.Zostera.shoots.core.imputed
                         # + log10.Zostera.sheath.length
                         + log10.Zostera.longest.leaf.length
                         # + AllelicRichness
                         # + log10.Leaf.PercN.imputed
                         # + log10.mesograzer.mass.per.g.plant.imputed
                         # + log10.crustacean.mass.per.g.plant.imputed
                         # + log10.gastropod.mass.per.g.plant.imputed
                         # + grazer.richness.site
                         # # + squid.survival.24hr.imputed
                         # + amphipod.survival.24hr # many missing value because only done at some sites
                         , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.2)
# N = 980
# AIC = 477.3763
sem.coefs(list(periphyton.plot.2), ZEN_2014_master_data_imputed, standardize = "scale")
#                                      response                         predictor    estimate  std.error p.value
# 3 log10.periphyton.mass.per.g.zostera.imputed log10.Zostera.longest.leaf.length  0.21786708 0.04081790  0.0000
# 2 log10.periphyton.mass.per.g.zostera.imputed log10.Zostera.shoots.core.imputed  0.03425733 0.02565333  0.1821
# 1 log10.periphyton.mass.per.g.zostera.imputed     log10.Zostera.AG.mass.imputed -0.01059450 0.02204977  0.6310


# Hierarchical model (plot level): PRODUCTIVITY
periphyton.plot.3 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                           #   Latitude
                           # + Temperature.C
                           # + Salinity.ppt
                           # # + Day.length.hours
                           # + log10.mean.fetch
                           # # + Depth.Categorical
                           # + pop.density.2015
                         #   + log10.Zostera.AG.mass.imputed
                         # + log10.Zostera.shoots.core.imputed
                         + log10.Zostera.sheath.length
                         # + log10.Zostera.longest.leaf.length
                         # + AllelicRichness
                         + log10.Leaf.PercN.imputed
                         # + log10.mesograzer.mass.per.g.plant.imputed
                         # + log10.crustacean.mass.per.g.plant.imputed
                         # + log10.gastropod.mass.per.g.plant.imputed
                         # + grazer.richness.site
                         # # + squid.survival.24hr.imputed
                         # + amphipod.survival.24hr # many missing value because only done at some sites
                         , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.3)
# N = 980
# AIC = 471.4786
sem.coefs(list(periphyton.plot.3), ZEN_2014_master_data_imputed, standardize = "scale")
#                                      response                   predictor   estimate  std.error p.value
# 1 log10.periphyton.mass.per.g.zostera.imputed log10.Zostera.sheath.length 0.20260704 0.03817256   0.000
# 2 log10.periphyton.mass.per.g.zostera.imputed    log10.Leaf.PercN.imputed 0.02647563 0.02505656   0.291


# Hierarchical model (plot level): SEAGRASS STRUCTURE + PRODUCTIVITY
periphyton.plot.4 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                           #   Latitude
                           # + Temperature.C
                           # + Salinity.ppt
                           # # + Day.length.hours
                           # + log10.mean.fetch
                           # # + Depth.Categorical
                           # + pop.density.2015
                           + log10.Zostera.AG.mass.imputed
                         + log10.Zostera.shoots.core.imputed
                         + log10.Zostera.sheath.length
                         + log10.Zostera.longest.leaf.length
                         # + AllelicRichness
                         + log10.Leaf.PercN.imputed
                         # + log10.mesograzer.mass.per.g.plant.imputed
                         # + log10.crustacean.mass.per.g.plant.imputed
                         # + log10.gastropod.mass.per.g.plant.imputed
                         # + grazer.richness.site
                         # # + squid.survival.24hr.imputed
                         # + amphipod.survival.24hr # many missing value because only done at some sites
                         , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.4)
# N = 980
# AIC = 476.8926
sem.coefs(list(periphyton.plot.4), ZEN_2014_master_data_imputed, standardize = "scale")
#                                      response                         predictor     estimate  std.error p.value
# 4 log10.periphyton.mass.per.g.zostera.imputed log10.Zostera.longest.leaf.length  0.150831335 0.04799111  0.0017
# 3 log10.periphyton.mass.per.g.zostera.imputed       log10.Zostera.sheath.length  0.126014777 0.04501492  0.0052
# 2 log10.periphyton.mass.per.g.zostera.imputed log10.Zostera.shoots.core.imputed  0.032037090 0.02556849  0.2105


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
periphyton.plot.5 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                                            #   Latitude
                                            + Temperature.C
                                          + Salinity.ppt
                                          + Day.length.hours
                                          + log10.mean.fetch
                                          # + Depth.Categorical
                                          # + pop.density.2015
                                          # + log10.Zostera.AG.mass.imputed
                                          # + log10.Zostera.shoots.core.imputed
                                          # + log10.Zostera.sheath.length
                                          # + log10.Zostera.longest.leaf.length
                                          # + AllelicRichness
                                          # + log10.periphyton.mass.per.g.zostera.imputed
                                          # + log10.Leaf.PercN.imputed
                                          # # + log10.mesograzer.mass.per.g.plant.imputed
                                          # # + log10.crustacean.mass.per.g.plant.imputed
                                          # # + log10.gastropod.mass.per.g.plant.imputed
                                          # # + grazer.richness.site
                                          # # + squid.survival.24hr.imputed
                                          # # + amphipod.survival.24hr # many missing value because only done at some sites
                                          , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.5)
# N = 980
# AIC = 511.0036
sem.coefs(list(periphyton.plot.5), ZEN_2014_master_data_imputed, standardize = "scale")
#                                      response        predictor    estimate std.error p.value
# 3 log10.periphyton.mass.per.g.zostera.imputed Day.length.hours -0.16745299 0.1587770  0.2978
# 4 log10.periphyton.mass.per.g.zostera.imputed log10.mean.fetch -0.14392084 0.1377173  0.3021
# 2 log10.periphyton.mass.per.g.zostera.imputed     Salinity.ppt  0.12799201 0.1428670  0.3755
# 1 log10.periphyton.mass.per.g.zostera.imputed    Temperature.C -0.07348919 0.1392671  0.6006


# Hierarchical model (plot level): ZOSTERA GENETICS
periphyton.plot.6 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                           #   Latitude
                           # + Temperature.C
                           # + Salinity.ppt
                           # + Day.length.hours
                           # + log10.mean.fetch
                           # + Depth.Categorical
                           # + pop.density.2015
                           # + log10.Zostera.AG.mass.imputed
                           # + log10.Zostera.shoots.core.imputed
                           # + log10.Zostera.sheath.length
                           # + log10.Zostera.longest.leaf.length
                         + AllelicRichness
                         # + log10.periphyton.mass.per.g.zostera.imputed
                         # + log10.Leaf.PercN.imputed
                         # # + log10.mesograzer.mass.per.g.plant.imputed
                         # # + log10.crustacean.mass.per.g.plant.imputed
                         # # + log10.gastropod.mass.per.g.plant.imputed
                         # # + grazer.richness.site
                         # # + squid.survival.24hr.imputed
                         # # + amphipod.survival.24hr # many missing value because only done at some sites
                         , random = ~1 | Coast/Site
                         , na.action = na.omit
                         , data = ZEN_2014_master_data_imputed)
summary(periphyton.plot.6)
# N = 980
# AIC = 493.7975
sem.coefs(list(periphyton.plot.6), ZEN_2014_master_data_imputed, standardize = "scale")
#                                      response       predictor  estimate std.error p.value
# 1 log10.periphyton.mass.per.g.zostera.imputed AllelicRichness 0.1738336 0.1247322  0.1704


# SUMMARY OF RESULTS: PLOT level
# Best model (AIC = 471.4786) is productivity with strong positive effect of sheath length 
# (proxy for Zostera growth rate). Next best is seagrass structure + productivity model, 
# which adds a highly significant positive effect of leaf length. These results seem best 
# interpreted as simple bottom-up control: periphyton biomass accumulation responds to the
# same factors that stimulate Zostera growth and leaf elongation. Note that this interpretation 
# is opposite the hypothesis suggested at workshop that periphyton biomass should be greatest on 
# slow-growing leaves.


###################################################################################
# MODELS OF PERIPHYTON BIOMASS: SITE LEVEL (SIMPLE)                               #
###################################################################################

# Linear model (site level): SEAGRASS STRUCTURE
periphyton.lm.1 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      # + log10.Zostera.AG.mass.site 
                      , data = ZEN_2014_site_means_49)
summary(periphyton.lm.1)
#                                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                -4.1074     1.1706  -3.509  0.00102 ** 
#   log10.Zostera.longest.leaf.length.cm.site   1.6922     0.3866   4.377 6.86e-05 ***
#   log10.Zostera.shoots.core.site              0.1117     0.2720   0.411  0.68325    
# 
# Residual standard error: 0.6089 on 46 degrees of freedom
# Multiple R-squared:  0.3534,	Adjusted R-squared:  0.3253 
# F-statistic: 12.57 on 2 and 46 DF,  p-value: 4.417e-05
AIC(periphyton.lm.1) # 95.3359


# Linear model (site level): PRODUCTIVITY
periphyton.lm.2 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                        log10.Zostera.sheath.length.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(periphyton.lm.2)
#                                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      -2.49090    0.35580  -7.001 9.09e-09 ***
#   log10.Zostera.sheath.length.site  1.55954    0.26174   5.958 3.34e-07 ***
#   Leaf.PercN.site                  -0.02562    0.15671  -0.163    0.871    
# 
# Residual standard error: 0.5587 on 46 degrees of freedom
# Multiple R-squared:  0.4555,	Adjusted R-squared:  0.4319 
# F-statistic: 19.24 on 2 and 46 DF,  p-value: 8.46e-07
AIC(periphyton.lm.2) # 86.90933


# Linear model (site level): SEAGRASS STRUCTURE + PRODUCTIVITY
periphyton.lm.3 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      + log10.Zostera.sheath.length.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(periphyton.lm.3)
#                                           Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                               -2.20753    1.26066  -1.751  0.08690 . 
# log10.Zostera.longest.leaf.length.cm.site -0.81855    0.86841  -0.943  0.35104   
# log10.Zostera.shoots.core.site             0.16149    0.25171   0.642  0.52448   
# log10.Zostera.sheath.length.site           2.33199    0.74345   3.137  0.00304 **
# Leaf.PercN.site                           -0.06309    0.16048  -0.393  0.69612   
# 
# Residual standard error: 0.5615 on 44 degrees of freedom
# Multiple R-squared:  0.474,	Adjusted R-squared:  0.4262 
# F-statistic: 9.914 on 4 and 44 DF,  p-value: 8.3e-06
AIC(periphyton.lm.3) # 89.21539


# Linear model (site level): ABIOTIC ENVIRONMENT
periphyton.lm.4 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                        Temperature.C 
                      + Salinity.ppt
                      + Day.length.hours 
                      + log10.mean.fetch
                      , data = ZEN_2014_site_means_49)
summary(periphyton.lm.4)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       0.97226    1.15835   0.839   0.4058  
# Temperature.C    -0.02853    0.02136  -1.335   0.1887  
# Salinity.ppt      0.01865    0.01468   1.270   0.2106  
# Day.length.hours -0.12077    0.05017  -2.407   0.0203 *
#   log10.mean.fetch -0.18314    0.24366  -0.752   0.4563  
# 
# Residual standard error: 0.6817 on 44 degrees of freedom
# Multiple R-squared:  0.2247,	Adjusted R-squared:  0.1542 
# F-statistic: 3.188 on 4 and 44 DF,  p-value: 0.022
AIC(periphyton.lm.4) # 108.2279


# Linear model (site level): ZOSTERA GENETICS
periphyton.lm.5 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                        log10.AllelicRichness 
                      , data = ZEN_2014_site_means_49)
summary(periphyton.lm.5)
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -1.9095     0.4263  -4.480 4.76e-05 ***
#   log10.AllelicRichness   1.4242     0.6338   2.247   0.0294 *  
# 
# Residual standard error: 0.7118 on 47 degrees of freedom
# Multiple R-squared:  0.09702,	Adjusted R-squared:  0.0778 
# F-statistic:  5.05 on 1 and 47 DF,  p-value: 0.02937
AIC(periphyton.lm.4) # 108.2279


# Linear model (site level): GEOGRAPHY
periphyton.lm.6 <- lm(log10.periphyton.mass.per.g.zostera.site ~ 
                         Latitude
                       + Ocean 
                       + Coast
                       , data = ZEN_2014_site_means_49)
summary(periphyton.lm.6)
#                    Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         0.04525    0.63823   0.071    0.944  
# Latitude           -0.02734    0.01161  -2.355    0.023 *
#   OceanPacific        0.40574    0.34156   1.188    0.241  
# CoastEast.Pacific   0.18169    0.29505   0.616    0.541  
# CoastWest.Atlantic -0.11028    0.27402  -0.402    0.689  
# CoastWest.Pacific        NA         NA      NA       NA  
# 
# Residual standard error: 0.6014 on 44 degrees of freedom
# Multiple R-squared:  0.3965,	Adjusted R-squared:  0.3416 
# F-statistic: 7.226 on 4 and 44 DF,  p-value: 0.0001456
AIC(periphyton.lm.6) # 95.95605

# RESULTS OF SITE LEVEL ANALYSIS  
# Strongest model is productivity (AIC = 86.90933), with strong positive effect of 
# Zostera sheath length (proxy for Zostera productivity rate). Reasonably close model 
# is seagrass structure + productivity, which however does not add any significant predictors. 
# As is true of the model predicting Zostera productivity, the close statistical correlation
# between Zostera productivity and periphyton biomass presumably reflects favorable
# growth conditions for both. 


###################################################################################
# MODELS OF CRUSTACEAN BIOMASS: RANDOM FORESTS                                    #
###################################################################################

# Explore graphically

pairs.panels(ZEN_2014_site_means[,c("log10.crustacean.mass.per.g.plant.site", "log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
  "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "log10.mean.fetch", "pop.density.2015",
  "log10.periphyton.mass.per.g.zostera.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 3) 

# Construct random forest: CRUSTACEAN BIOMASS (with geographic variables)
RF.crustacean.1 = randomForest(
  log10.crustacean.mass.per.g.plant ~ 
    Ocean
  + Coast
  + basin
  + Latitude
  + Longitude
  + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass
  + log10.macrophytes.total.AG.mass.core
  + log10.Zostera.shoots.core
  + Zostera.sheath.width
  + Zostera.sheath.length
  + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera  
  #  + log10.mesograzer.abund.per.g.plant
  # + log10.crustacean.abund.per.g.plant
  # + log10.amphipod.abund.per.g.plant
  # + log10.gammarid.abund.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.crustacean.1 # % Var explained: 73.51

# Plot variable importance
varImpPlot(RF.crustacean.1) # Look at % Inc MSE (first plot) not Inc Node Purity
# RESULTS: Longitude clearly separates as strongest predictor


# Construct random forest: CRUSTACEAN BIOMASS (omitting geographic variables)
RF.crustacean.2 = randomForest(
  log10.crustacean.mass.per.g.plant ~ 
    #  Ocean
    # + Coast
    # + Latitude
    # + Longitude
    # + basin 
    + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass
  + log10.macrophytes.total.AG.mass.core
  + log10.Zostera.shoots.core
  + Zostera.sheath.width
  + Zostera.sheath.length
  + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera  
  #  + log10.mesograzer.abund.per.g.plant
  # + log10.crustacean.abund.per.g.plant
  # + log10.amphipod.abund.per.g.plant
  # + log10.gammarid.abund.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.crustacean.2 # % Var explained: 72.54

# Plot variable importance
varImpPlot(RF.crustacean.2) 
# RESULT: Allelic richness, temperature, fetch, allelic richness ...), Go figure ...


###################################################################################
# MODELS OF CRUSTACEAN BIOMASS: PLOT LEVEL (HIERARCHICAL)                         #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable
crustaceans1 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                     Latitude
                   + Temperature.C
                   + Salinity.ppt
                   #+ Day.length.hours
                   + log10.mean.fetch
                   # + Depth.Categorical
                   + pop.density.2015
                   + log10.Zostera.AG.mass.imputed
                   + log10.Zostera.shoots.core.imputed
                   + log10.Zostera.sheath.length
                   + log10.Zostera.longest.leaf.length
                   + AllelicRichness
                   + log10.periphyton.mass.per.g.zostera.imputed
                   + log10.Leaf.PercN.imputed
                   # + log10.mesograzer.mass.per.g.plant.imputed
                   # + log10.crustacean.mass.per.g.plant.imputed
                   # + log10.gastropod.mass.per.g.plant.imputed
                   # + grazer.richness.site
                   # + squid.survival.24hr.imputed
                   # + amphipod.survival.24hr # many missing value because only done at some sites
                   , random = ~1 | Coast/Site
                   , na.action = na.omit
                   , data = ZEN_2014_master_data_imputed)
summary(crustaceans1)
# N = 980
# AIC = 1450.742
sem.coefs(list(crustaceans1), ZEN_2014_master_data_imputed, standardize = "scale")
#                                     response                                   predictor    estimate  std.error p.value
# 9  log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.longest.leaf.length -0.18014958 0.06629213  0.0067
# 11 log10.crustacean.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed  0.10399288 0.04390389  0.0181
# 7  log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.shoots.core.imputed -0.07883481 0.03513792  0.0251
# 2  log10.crustacean.mass.per.g.plant.imputed                               Temperature.C  0.28029351 0.13351867  0.0423
# 3  log10.crustacean.mass.per.g.plant.imputed                                Salinity.ppt -0.23931032 0.14200080  0.0999


# Hierarchical model (plot level): SEAGRASS STRUCTURE
crustaceans2 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                    #   Latitude
                    # + Temperature.C
                    # + Salinity.ppt
                    # #+ Day.length.hours
                    # + log10.mean.fetch
                    # # + Depth.Categorical
                    # + pop.density.2015
                    + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(crustaceans2)
# N = 980
# AIC = 1409.433
sem.coefs(list(crustaceans2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                    response                         predictor    estimate  std.error p.value
# 3 log10.crustacean.mass.per.g.plant.imputed log10.Zostera.longest.leaf.length -0.13731635 0.05521691  0.0131
# 2 log10.crustacean.mass.per.g.plant.imputed log10.Zostera.shoots.core.imputed -0.07164030 0.03458990  0.0386
# 1 log10.crustacean.mass.per.g.plant.imputed     log10.Zostera.AG.mass.imputed  0.04613211 0.03064740  0.1326


# Hierarchical model (plot level): PRODUCTIVITY
crustaceans3 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                    #   Latitude
                    # + Temperature.C
                    # + Salinity.ppt
                    # #+ Day.length.hours
                    # + log10.mean.fetch
                    # # + Depth.Categorical
                    # + pop.density.2015
                    # + log10.Zostera.AG.mass.imputed
                    # + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(crustaceans3)
# N = 980
# AIC = 1410.553
sem.coefs(list(crustaceans3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                    response                                   predictor   estimate  std.error p.value
# 1 log10.crustacean.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed 0.07662081 0.04266336  0.0728
# 2 log10.crustacean.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed 0.03137279 0.03334154  0.3470


# Hierarchical model (plot level): SEAGRASS STRUCTURE + PRODUCTIVITY
crustaceans4 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(crustaceans4)
# N = 980
# AIC = 1412.368

sem.coefs(list(crustaceans4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                    response                                   predictor    estimate  std.error p.value
# 3 log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.longest.leaf.length -0.15832154 0.05597971  0.0048
# 4 log10.crustacean.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed  0.10077984 0.04329399  0.0201
# 2 log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.shoots.core.imputed -0.07505989 0.03453282  0.0300
# 1 log10.crustacean.mass.per.g.plant.imputed               log10.Zostera.AG.mass.imputed  0.04801576 0.03060133  0.1170
# 5 log10.crustacean.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed  0.02851642 0.03322925  0.3910


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
crustaceans5 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      + Temperature.C
                      + Salinity.ppt
                      + Day.length.hours
                      + log10.mean.fetch
                      # + Depth.Categorical
                      # + pop.density.2015
                      # + log10.Zostera.AG.mass.imputed
                    # + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(crustaceans5)
# N = 980
# AIC = 1418.446
sem.coefs(list(crustaceans5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                    response        predictor    estimate std.error p.value
# 1 log10.crustacean.mass.per.g.plant.imputed    Temperature.C  0.27023655 0.1048581  0.0136
# 2 log10.crustacean.mass.per.g.plant.imputed     Salinity.ppt -0.21786890 0.1221374  0.0817
# 3 log10.crustacean.mass.per.g.plant.imputed log10.mean.fetch -0.06455769 0.1179343  0.5870


# Hierarchical model (plot level): ZOSTERA GENETICS
crustaceans6 <- lme(log10.crustacean.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                    # + Salinity.ppt
                    # + Day.length.hours
                    # + log10.mean.fetch
                    # + Depth.Categorical
                    # + pop.density.2015
                    # + log10.Zostera.AG.mass.imputed
                    # + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    # + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(crustaceans6)
# N = 980
# AIC = 1410.625
sem.coefs(list(crustaceans6), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                    response       predictor   estimate std.error p.value
# 1 log10.crustacean.mass.per.g.plant.imputed AllelicRichness 0.05010693 0.1144681  0.6637

x <- list(crustaceans1, crustaceans2, crustaceans3, crustaceans4, crustaceans5, crustaceans6)
sem.model.fits(x) # This outputs an AIC table 


# SUMMARY OF RESULTS: PLOT level
# Best model (AIC = 1409.433) is seagrass structure: crustacean biomass higher in sparser, 
# shorter-leaved eelgrass stands. Close second is productivity model (AIC = 1410.553) with
# marginally greater positive effect of periphyton biomass; and genetics model (AIC = 1410.625),
# which paradoxically is not even close to significant ... Main message I draw from this is 
# that at PLOT level (within a site, controlling for variation among sites), crustacean biomass
# tends to be higher in sparse, short-leaved clumps with periphyton. This presumably reflects
# some combination of attraction of mobile epifauna to epiphytic food and accumulationg on 
# nearest structure in a patchy landscape. Note that results differ at site level (see below).  


###################################################################################
# MODELS OF CRUSTACEAN BIOMASS: SITE LEVEL (SIMPLE)                               #
###################################################################################

# Linear model (site level): SEAGRASS STRUCTURE
crustacean.lm.1 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                log10.Zostera.longest.leaf.length.cm.site 
             + log10.Zostera.shoots.core.site
             # + log10.Zostera.AG.mass.site 
             , data = ZEN_2014_site_means_49)
summary(crustacean.lm.1)
#                                           Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                -4.0007     1.4722  -2.718  0.00924 **
#   log10.Zostera.longest.leaf.length.cm.site   1.5092     0.4862   3.104  0.00326 **
#   log10.Zostera.shoots.core.site              0.6578     0.3421   1.923  0.06072 . 
# 
# Residual standard error: 0.7657 on 46 degrees of freedom
# Multiple R-squared:  0.174,	Adjusted R-squared:  0.138 
# F-statistic: 4.844 on 2 and 46 DF,  p-value: 0.01233
AIC(crustacean.lm.1) # 117.8012


# Linear model (site level): PRODUCTIVITY
crustacean.lm.2 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                        log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(crustacean.lm.2)
#                                          Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                               -0.5551     0.4722  -1.176   0.2458  
# log10.periphyton.mass.per.g.zostera.site   0.3199     0.1467   2.180   0.0344 *
# Leaf.PercN.site                            0.5066     0.2013   2.516   0.0154 *
# 
# Residual standard error: 0.74 on 46 degrees of freedom
# Multiple R-squared:  0.2285,	Adjusted R-squared:  0.195 
# F-statistic: 6.814 on 2 and 46 DF,  p-value: 0.002559
AIC(crustacean.lm.2) # 114.451


# Linear model (site level): SEAGRASS STRUCTURE + PRODUCTIVITY
crustacean.lm.3 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      + log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(crustacean.lm.3)
#                                           Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                -3.5925     1.5868  -2.264   0.0286 *
# log10.Zostera.longest.leaf.length.cm.site   0.9087     0.5523   1.645   0.1070  
# log10.Zostera.shoots.core.site              0.6008     0.3241   1.854   0.0705 .
# log10.periphyton.mass.per.g.zostera.site    0.2176     0.1756   1.239   0.2220  
# Leaf.PercN.site                             0.4691     0.1990   2.357   0.0229 *
# 
# Residual standard error: 0.7235 on 44 degrees of freedom
# Multiple R-squared:  0.2947,	Adjusted R-squared:  0.2305 
# F-statistic: 4.595 on 4 and 44 DF,  p-value: 0.003459
AIC(crustacean.lm.3) # 114.0613


# Linear model (site level): ABIOTIC ENVIRONMENT
crustacean.lm.4 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                        Temperature.C 
                      + Salinity.ppt
                      + Day.length.hours 
                      + log10.mean.fetch
                      , data = ZEN_2014_site_means_49)
summary(crustacean.lm.4)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       2.39042    1.32285   1.807   0.0776 .
# Temperature.C     0.01251    0.02440   0.513   0.6108  
# Salinity.ppt     -0.01079    0.01677  -0.643   0.5233  
# Day.length.hours -0.14197    0.05730  -2.478   0.0171 *
#   log10.mean.fetch  0.11382    0.27826   0.409   0.6845  
# 
# Residual standard error: 0.7785 on 44 degrees of freedom
# Multiple R-squared:  0.1833,	Adjusted R-squared:  0.1091 
# F-statistic: 2.469 on 4 and 44 DF,  p-value: 0.05845
AIC(crustacean.lm.4) # 121.2415


# Linear model (site level): ZOSTERA GENETICS
crustacean.lm.5 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                        log10.AllelicRichness 
                      , data = ZEN_2014_site_means_49)
summary(crustacean.lm.5)
#                       Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            -0.6555     0.4836  -1.355   0.1818  
# log10.AllelicRichness   1.2576     0.7191   1.749   0.0869 .
# 
# Residual standard error: 0.8076 on 47 degrees of freedom
# Multiple R-squared:  0.0611,	Adjusted R-squared:  0.04112 
# F-statistic: 3.058 on 1 and 47 DF,  p-value: 0.08685
AIC(crustacean.lm.5) # 122.0764


# SUMMARY OF RESULTS: SITE level
# Best model (AIC = 114.0613) predicting crustcean biomass is seagrass structure + productivity, 
# showing positive effects of leaf N and marginally positive effect of shoot density. Very close model 
# (AIC = 114.451) is productivity model, which adds positive effect of periphyton biomass.
# These results seem to clearly indicate bottom-up resource control of crustacean biomass 
# by nitrogen availability fueling periphyton growth. 


###################################################################################
# MODELS OF GAMMARID BIOMASS: RANDOM FORESTS                                      #
###################################################################################

# Explore graphically

pairs.panels(ZEN_2014_site_means[,c("log10.gammarid.mass.per.g.plant.site", "log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
  "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "log10.mean.fetch", "pop.density.2015",
  "log10.periphyton.mass.per.g.zostera.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 3) 

# NOTE: even log10 transform of gammarid mass is highly non-normal (right-skewed).
hist(ZEN_2014_site_means$log10.gammarid.mass.per.g.plant.site)

# Construct random forest: GAMMARID BIOMASS (with geographic variables)
RF.gammarid.1 = randomForest(log10.gammarid.mass.per.g.plant ~ 
    Ocean
  + Coast
  + basin
  + Latitude
  + Longitude
  + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass
  + log10.macrophytes.total.AG.mass.core
  + log10.Zostera.shoots.core
  + Zostera.sheath.width
  + Zostera.sheath.length
  + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera  
  #  + log10.mesograzer.abund.per.g.plant
  # + log10.crustacean.abund.per.g.plant
  # + log10.amphipod.abund.per.g.plant
  # + log10.gammarid.abund.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.gammarid.1 # % Var explained: 68.47

# Plot variable importance
varImpPlot(RF.gammarid.1) 
# RESULTS: Latitude and Longitude clearly separate as strongest predictors


# Construct random forest: GAMMARID BIOMASS (omitting geographic variables)
RF.gammarid.2 = randomForest(log10.gammarid.mass.per.g.plant ~ 
    #  Ocean
    # + Coast
    # + Latitude
    # + Longitude
    # + basin 
    + Temperature.C
  + Salinity.ppt
  # + Day.length.hours
  + log10.mean.fetch
  + pop.density.2015
  + Depth.Categorical
  # + Depth.m
  # + Perc.Silt
  # + Perc.Sand
  # + Perc.Gravel
  # + GenotypicRichness
  + AllelicRichness
  # + Inbreeding
  + log10.Zostera.AG.mass
  + log10.macrophytes.total.AG.mass.core
  + log10.Zostera.shoots.core
  + Zostera.sheath.width
  + Zostera.sheath.length
  + Zostera.longest.leaf.length
  # + log10.pct.cover.macroalgae
  # + pct.cover.seagrass
  + log10.periphyton.mass.per.g.zostera  
  #  + log10.mesograzer.abund.per.g.plant
  # + log10.crustacean.abund.per.g.plant
  # + log10.amphipod.abund.per.g.plant
  # + log10.gammarid.abund.per.g.plant
  #   + log10.gastropod.abund.per.g.plant
  #   + log10.grazer.richness.plot
  #   + log10.grazer.richness.site
  + amphipod.survival.24hr
  + log10.Leaf.PercN
  # + log10.GenotypicRichness
  ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.gammarid.2 # % Var explained: 67.46

# Plot variable importance
varImpPlot(RF.gammarid.2) 
# RESULT: Temperature clearly strongest.


###################################################################################
# MODELS OF GAMMARID BIOMASS: PLOT LEVEL (HIERARCHICAL)                         #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable
gammarid.plot.1 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      Latitude
                    + Temperature.C
                    + Salinity.ppt
                    #+ Day.length.hours
                    + log10.mean.fetch
                    # + Depth.Categorical
                    + pop.density.2015
                    + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # + log10.mesograzer.mass.per.g.plant.imputed
                    # + log10.crustacean.mass.per.g.plant.imputed
                    # + log10.gastropod.mass.per.g.plant.imputed
                    # + grazer.richness.site
                    # + squid.survival.24hr.imputed
                    # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_imputed)
summary(gammarid.plot.1)
# N = 980
# AIC = 1295.423
sem.coefs(list(gammarid.plot.1), ZEN_2014_master_data_imputed, standardize = "scale")
#                                   response                                   predictor     estimate  std.error p.value
# 11 log10.gammarid.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed  0.161519778 0.04713476  0.0006
# 1  log10.gammarid.mass.per.g.plant.imputed                                    Latitude -0.543912240 0.16156342  0.0017
# 3  log10.gammarid.mass.per.g.plant.imputed                                Salinity.ppt -0.374418171 0.14987136  0.0168
# 5  log10.gammarid.mass.per.g.plant.imputed                            pop.density.2015 -0.236634940 0.12571699  0.0673
# 7  log10.gammarid.mass.per.g.plant.imputed           log10.Zostera.shoots.core.imputed -0.057375339 0.03747273  0.1261


# Hierarchical model (plot level): SEAGRASS STRUCTURE
gammarid.plot.2 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gammarid.plot.2)
# N = 980
# AIC = 1258.341
sem.coefs(list(gammarid.plot.2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                  response                         predictor     estimate  std.error p.value
# 2 log10.gammarid.mass.per.g.plant.imputed log10.Zostera.shoots.core.imputed -0.043984362 0.03710272  0.2361
# 3 log10.gammarid.mass.per.g.plant.imputed log10.Zostera.longest.leaf.length -0.057392070 0.05895979  0.3306
# 1 log10.gammarid.mass.per.g.plant.imputed     log10.Zostera.AG.mass.imputed  0.003216865 0.03287751  0.9221


# Hierarchical model (plot level): PRODUCTIVITY
gammarid.plot.3 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      # + log10.Zostera.AG.mass.imputed
                      # + log10.Zostera.shoots.core.imputed
                      # + log10.Zostera.sheath.length
                      # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gammarid.plot.3)
# N = 980
# AIC = 1244.814
sem.coefs(list(gammarid.plot.3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                  response                                   predictor   estimate  std.error p.value
# 1 log10.gammarid.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed 0.14631661 0.04534972  0.0013
# 2 log10.gammarid.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed 0.01990135 0.03552004  0.5754


# Hierarchical model (plot level): SEAGRASS STRUCTURE + PRODUCTIVITY
gammarid.plot.4 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gammarid.plot.4)
# N = 980
# AIC = 1255.58

sem.coefs(list(gammarid.plot.4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                  response                                   predictor     estimate  std.error p.value
# 4 log10.gammarid.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed  0.160905493 0.04629921  0.0005
# 3 log10.gammarid.mass.per.g.plant.imputed           log10.Zostera.longest.leaf.length -0.092380932 0.05949203  0.1208
# 2 log10.gammarid.mass.per.g.plant.imputed           log10.Zostera.shoots.core.imputed -0.048315877 0.03690081  0.1907
# 5 log10.gammarid.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed  0.017411895 0.03555180  0.6244
# 1 log10.gammarid.mass.per.g.plant.imputed               log10.Zostera.AG.mass.imputed  0.005597302 0.03273338  0.8643


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
gammarid.plot.5 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      + Temperature.C
                    + Salinity.ppt
                    + Day.length.hours
                    + log10.mean.fetch
                    # + Depth.Categorical
                    # + pop.density.2015
                    # + log10.Zostera.AG.mass.imputed
                    # + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gammarid.plot.5)
# N = 980
# AIC = 1266.779
sem.coefs(list(gammarid.plot.5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                  response        predictor    estimate std.error p.value
# 1 log10.gammarid.mass.per.g.plant.imputed    Temperature.C  0.23563472 0.1354332  0.0894
# 2 log10.gammarid.mass.per.g.plant.imputed     Salinity.ppt -0.21087144 0.1388147  0.1364
# 4 log10.gammarid.mass.per.g.plant.imputed log10.mean.fetch -0.12422710 0.1333891  0.3571
# 3 log10.gammarid.mass.per.g.plant.imputed Day.length.hours -0.06723316 0.1544346  0.6656


# Hierarchical model (plot level): ZOSTERA GENETICS
gammarid.plot.6 <- lme(log10.gammarid.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # + Day.length.hours
                      # + log10.mean.fetch
                      # + Depth.Categorical
                      # + pop.density.2015
                      # + log10.Zostera.AG.mass.imputed
                      # + log10.Zostera.shoots.core.imputed
                      # + log10.Zostera.sheath.length
                      # + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gammarid.plot.6)
# N = 980
# AIC = 1251.041
sem.coefs(list(gammarid.plot.6), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                  response       predictor   estimate std.error p.value
# 1 log10.gammarid.mass.per.g.plant.imputed AllelicRichness 0.06553935 0.1237555  0.5991


# SUMMARY OF RESULTS: PLOT level
# Results in this case seem quite clear: at the level of plots within a site, gammarids 
# respond strongly to availability of periphyton. There is little evidence for significant
# effects of any other predictor. 


###################################################################################
# MODELS OF GAMMARID BIOMASS: SITE LEVEL (SIMPLE)                                 #
###################################################################################

# Linear model (site level): SEAGRASS STRUCTURE
gammarid.lm.1 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      # + log10.Zostera.AG.mass.site 
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.1)
#                                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                -4.5733     1.2079  -3.786 0.000442 ***
# log10.Zostera.longest.leaf.length.cm.site   1.3510     0.3989   3.386 0.001459 ** 
# log10.Zostera.shoots.core.site              0.8241     0.2807   2.936 0.005180 ** 
# 
# Residual standard error: 0.6283 on 46 degrees of freedom
# Multiple R-squared:  0.2211,	Adjusted R-squared:  0.1873 
# F-statistic:  6.53 on 2 and 46 DF,  p-value: 0.003189
AIC(gammarid.lm.1) # 98.4112


# Linear model (site level): PRODUCTIVITY
gammarid.lm.2 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.2)
#                                          Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                               -0.3600     0.4102  -0.878  0.38474   
# log10.periphyton.mass.per.g.zostera.site   0.3449     0.1275   2.706  0.00952 **
# Leaf.PercN.site                            0.2129     0.1749   1.217  0.22980   
# 
# Residual standard error: 0.6429 on 46 degrees of freedom
# Multiple R-squared:  0.1846,	Adjusted R-squared:  0.1491 
# F-statistic: 5.205 on 2 and 46 DF,  p-value: 0.009164
AIC(gammarid.lm.2) # 100.66


# Linear model (site level): SEAGRASS STRUCTURE + PRODUCTIVITY
gammarid.lm.3 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      + log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.3)
#                                           Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                -3.5221     1.3220  -2.664  0.01075 * 
# log10.Zostera.longest.leaf.length.cm.site   0.7431     0.4602   1.615  0.11350   
# log10.Zostera.shoots.core.site              0.7770     0.2700   2.877  0.00616 **
# log10.periphyton.mass.per.g.zostera.site    0.3037     0.1463   2.076  0.04379 * 
# Leaf.PercN.site                             0.1897     0.1658   1.144  0.25867   
# 
# Residual standard error: 0.6028 on 44 degrees of freedom
# Multiple R-squared:  0.3142,	Adjusted R-squared:  0.2519 
# F-statistic:  5.04 on 4 and 44 DF,  p-value: 0.001969
AIC(gammarid.lm.3) # 96.17355


# Linear model (site level): ABIOTIC ENVIRONMENT
gammarid.lm.4 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        Temperature.C 
                      + Salinity.ppt
                      + Day.length.hours 
                      + log10.mean.fetch
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.4)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       1.45224    1.13234   1.283   0.2064  
# Temperature.C     0.01808    0.02088   0.866   0.3914  
# Salinity.ppt     -0.01418    0.01435  -0.988   0.3284  
# Day.length.hours -0.10686    0.04905  -2.179   0.0347 *
#   log10.mean.fetch -0.09659    0.23819  -0.406   0.6871  
# 
# Residual standard error: 0.6664 on 44 degrees of freedom
# Multiple R-squared:  0.1619,	Adjusted R-squared:  0.08571 
# F-statistic: 2.125 on 4 and 44 DF,  p-value: 0.09368
AIC(gammarid.lm.4) # 106.0028


# Linear model (site level): ZOSTERA GENETICS
gammarid.lm.5 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        log10.AllelicRichness 
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.5)
#                       Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            -0.8059     0.4138  -1.948   0.0574 .
# log10.AllelicRichness   0.8310     0.6152   1.351   0.1833  
# 
# Residual standard error: 0.691 on 47 degrees of freedom
# Multiple R-squared:  0.03736,	Adjusted R-squared:  0.01688 
# F-statistic: 1.824 on 1 and 47 DF,  p-value: 0.1833
AIC(gammarid.lm.5) # 106.7912


# Linear model (site level): GEOGRAPHY
gammarid.lm.6 <- lm(log10.gammarid.mass.per.g.plant.site ~ 
                        Latitude
                      + Ocean
                      + Coast
                      , data = ZEN_2014_site_means_49)
summary(gammarid.lm.6)
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)         0.002866   0.667820   0.004    0.997
# Latitude           -0.013424   0.012145  -1.105    0.275
# OceanPacific        0.462935   0.357396   1.295    0.202
# CoastEast.Pacific   0.099777   0.308730   0.323    0.748
# CoastWest.Atlantic  0.462007   0.286723   1.611    0.114
# CoastWest.Pacific         NA         NA      NA       NA
# 
# Residual standard error: 0.6293 on 44 degrees of freedom
# Multiple R-squared:  0.2525,	Adjusted R-squared:  0.1845 
# F-statistic: 3.716 on 4 and 44 DF,  p-value: 0.01087
AIC(gammarid.lm.6) # 100.3971


# SUMMARY OF RESULTS: SITE level
# Clear results for gammarids at site scale, as for plot scale, and they are largely 
# concordant: best model is seagrass structure + productivity (AIC = 96.17355), with 
# significantly positive effects of shoot density and periphyton biomass. Plot-level
# analysis also showed strong positive effect of periphyton. By contrast there is little
# evidence of any effect of geography or abiotic factors. 


###################################################################################
# MODELS OF GASTROPOD BIOMASS: RANDOM FORESTS                                     #
###################################################################################

# Explore graphically

pairs.panels(ZEN_2014_site_means[,c("log10.gastropod.mass.per.g.plant.site", "log10.Zostera.AG.mass.site", "log10.Zostera.shoots.core.site",
  "AllelicRichness", "Latitude", "Temperature.C", "Salinity.ppt", "Leaf.PercN.site", "log10.mean.fetch", "pop.density.2015",
  "log10.periphyton.mass.per.g.zostera.site")],
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 2) 

# Construct random forest: GASTROPOD BIOMASS (with geographic variables)
RF.gastropod.1 = randomForest(log10.gastropod.mass.per.g.plant ~ 
                                Ocean
                              + Coast
                              + basin
                              + Latitude
                              + Longitude
                              + Temperature.C
                              + Salinity.ppt
                              # + Day.length.hours
                              + log10.mean.fetch
                              + pop.density.2015
                              + Depth.Categorical
                              # + Depth.m
                              # + Perc.Silt
                              # + Perc.Sand
                              # + Perc.Gravel
                              # + GenotypicRichness
                              + AllelicRichness
                              # + Inbreeding
                              + log10.Zostera.AG.mass
                              + log10.macrophytes.total.AG.mass.core
                              + log10.Zostera.shoots.core
                              + Zostera.sheath.width
                              + Zostera.sheath.length
                              + Zostera.longest.leaf.length
                              # + log10.pct.cover.macroalgae
                              # + pct.cover.seagrass
                              + log10.periphyton.mass.per.g.zostera  
                              #  + log10.mesograzer.abund.per.g.plant
                              # + log10.crustacean.abund.per.g.plant
                              # + log10.amphipod.abund.per.g.plant
                              # + log10.gammarid.abund.per.g.plant
                              #   + log10.gastropod.abund.per.g.plant
                              #   + log10.grazer.richness.plot
                              #   + log10.grazer.richness.site
                              + amphipod.survival.24hr
                              + log10.Leaf.PercN
                              # + log10.GenotypicRichness
                              ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.gastropod.1 # % Var explained: 66.62

# Plot variable importance
varImpPlot(RF.gastropod.1) # Temperature, latitude ...


# Construct random forest: GASTROPOD BIOMASS (omitting geographic variables)
RF.gastropod.2 = randomForest(log10.gastropod.mass.per.g.plant ~ 
                              #   Ocean
                              # + Coast
                              # + basin
                              # + Latitude
                              # + Longitude
                              + Temperature.C
                              + Salinity.ppt
                              # + Day.length.hours
                              + log10.mean.fetch
                              + pop.density.2015
                              + Depth.Categorical
                              # + Depth.m
                              # + Perc.Silt
                              # + Perc.Sand
                              # + Perc.Gravel
                              # + GenotypicRichness
                              + AllelicRichness
                              # + Inbreeding
                              + log10.Zostera.AG.mass
                              + log10.macrophytes.total.AG.mass.core
                              + log10.Zostera.shoots.core
                              + Zostera.sheath.width
                              + Zostera.sheath.length
                              + Zostera.longest.leaf.length
                              # + log10.pct.cover.macroalgae
                              # + pct.cover.seagrass
                              + log10.periphyton.mass.per.g.zostera  
                              #  + log10.mesograzer.abund.per.g.plant
                              # + log10.crustacean.abund.per.g.plant
                              # + log10.amphipod.abund.per.g.plant
                              # + log10.gammarid.abund.per.g.plant
                              #   + log10.gastropod.abund.per.g.plant
                              #   + log10.grazer.richness.plot
                              #   + log10.grazer.richness.site
                              + amphipod.survival.24hr
                              + log10.Leaf.PercN
                              # + log10.GenotypicRichness
                              ,
  data = ZEN_2014_master_data,
  na.action = na.roughfix,
  ntrees = 500,
  importance = T
)

# Examine summary output
RF.gastropod.2 # % Var explained: 65.83

# Plot variable importance
varImpPlot(RF.gastropod.2) 
# RESULT: Temperature is #1. 


###################################################################################
# MODELS OF GASTROPOD BIOMASS: PLOT LEVEL (HIERARCHICAL)                          #
###################################################################################

# Hierarchical model (plot level): FULL
# NOTE: Using coast as random variable
gastropods.plot.1 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      Latitude
                    + Temperature.C
                    + Salinity.ppt
                    #+ Day.length.hours
                    + log10.mean.fetch
                    # + Depth.Categorical
                    + pop.density.2015
                    + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # + log10.mesograzer.mass.per.g.plant.imputed
                    # + log10.crustacean.mass.per.g.plant.imputed
                    # + log10.gastropod.mass.per.g.plant.imputed
                    # + grazer.richness.site
                    # + squid.survival.24hr.imputed
                    # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_imputed)
summary(gastropods.plot.1)
# N = 980
# AIC = 1688.144
sem.coefs(list(gastropods.plot.1), ZEN_2014_master_data_imputed, standardize = "scale")
#                                    response                                   predictor    estimate  std.error p.value
# 2  log10.gastropod.mass.per.g.plant.imputed                               Temperature.C  0.30156061 0.14241543  0.0406
# 12 log10.gastropod.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed  0.06217693 0.03823934  0.1043


# Hierarchical model (plot level): SEAGRASS STRUCTURE
gastropods.plot.2 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gastropods.plot.2)
# N = 980
# AIC = 1642.289
sem.coefs(list(gastropods.plot.2), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                   response                         predictor    estimate  std.error p.value
# 1 log10.gastropod.mass.per.g.plant.imputed     log10.Zostera.AG.mass.imputed -0.04400171 0.03406952  0.1968
# 3 log10.gastropod.mass.per.g.plant.imputed log10.Zostera.longest.leaf.length -0.07479943 0.05988726  0.2120
# 2 log10.gastropod.mass.per.g.plant.imputed log10.Zostera.shoots.core.imputed -0.02608778 0.03833231  0.4963


# Hierarchical model (plot level): PRODUCTIVITY
gastropods.plot.3 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      # + log10.Zostera.AG.mass.imputed
                      # + log10.Zostera.shoots.core.imputed
                      # + log10.Zostera.sheath.length
                      # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gastropods.plot.3)
# N = 980
# AIC = 1637.9
sem.coefs(list(gastropods.plot.3), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                   response                                   predictor    estimate  std.error p.value
# 2 log10.gastropod.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed  0.06436014 0.03698237  0.0821
# 1 log10.gastropod.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed -0.05961671 0.04715044  0.2064


# Hierarchical model (plot level): SEAGRASS STRUCTURE + PRODUCTIVITY
gastropods.plot.4 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # #+ Day.length.hours
                      # + log10.mean.fetch
                      # # + Depth.Categorical
                      # + pop.density.2015
                      + log10.Zostera.AG.mass.imputed
                    + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    + log10.periphyton.mass.per.g.zostera.imputed
                    + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gastropods.plot.4)
# N = 980
# AIC = 1647.42
sem.coefs(list(gastropods.plot.4), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                   response                                   predictor    estimate  std.error p.value
# 5 log10.gastropod.mass.per.g.plant.imputed                    log10.Leaf.PercN.imputed  0.06171054 0.03696748  0.0954
# 1 log10.gastropod.mass.per.g.plant.imputed               log10.Zostera.AG.mass.imputed -0.04330713 0.03406708  0.2040


# Hierarchical model (plot level): ABIOTIC ENVIRONMENT
gastropods.plot.5 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      + Temperature.C
                    + Salinity.ppt
                    + Day.length.hours
                    + log10.mean.fetch
                    # + Depth.Categorical
                    # + pop.density.2015
                    # + log10.Zostera.AG.mass.imputed
                    # + log10.Zostera.shoots.core.imputed
                    # + log10.Zostera.sheath.length
                    # + log10.Zostera.longest.leaf.length
                    # + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gastropods.plot.5)
# N = 980
# AIC = 1650.255
sem.coefs(list(gastropods.plot.5), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                   response        predictor    estimate std.error p.value
# 1 log10.gastropod.mass.per.g.plant.imputed    Temperature.C  0.25931220 0.1355219  0.0627
# 3 log10.gastropod.mass.per.g.plant.imputed Day.length.hours  0.26150025 0.1470710  0.0828
# 2 log10.gastropod.mass.per.g.plant.imputed     Salinity.ppt -0.16816496 0.1413537  0.2410


# Hierarchical model (plot level): ZOSTERA GENETICS
gastropods.plot.6 <- lme(log10.gastropod.mass.per.g.plant.imputed ~ 
                      #   Latitude
                      # + Temperature.C
                      # + Salinity.ppt
                      # + Day.length.hours
                      # + log10.mean.fetch
                      # + Depth.Categorical
                      # + pop.density.2015
                      # + log10.Zostera.AG.mass.imputed
                      # + log10.Zostera.shoots.core.imputed
                      # + log10.Zostera.sheath.length
                      # + log10.Zostera.longest.leaf.length
                    + AllelicRichness
                    # + log10.periphyton.mass.per.g.zostera.imputed
                    # + log10.Leaf.PercN.imputed
                    # # + log10.mesograzer.mass.per.g.plant.imputed
                    # # + log10.crustacean.mass.per.g.plant.imputed
                    # # + log10.gastropod.mass.per.g.plant.imputed
                    # # + grazer.richness.site
                    # # + squid.survival.24hr.imputed
                    # # + amphipod.survival.24hr # many missing value because only done at some sites
                    , random = ~1 | Coast/Site
                    , na.action = na.omit
                    , data = ZEN_2014_master_data_49_imputed)
summary(gastropods.plot.6)
# N = 980
# AIC = 1636.848
sem.coefs(list(gastropods.plot.6), ZEN_2014_master_data_49_imputed, standardize = "scale")
#                                   response       predictor   estimate std.error p.value
# 1 log10.gastropod.mass.per.g.plant.imputed AllelicRichness -0.1794264 0.1235667  0.1536

# SUMMARY OF RESULTS: PLOT level
# I conclude that gastropod biomass is essentially unpredictable at the plot scale. The 'best' 
# model is allelic richness (AIC = 1636.848), whicH however, is not even close to significant. 
# It appears to have lowest AIC score merely because it has only a single predictor so fewer
# parameters than the other worthless models. The unpredictablity of gastropod biomass at plot 
# scale contrasts with the situation for crustaceans, which may reflect important differences in 
# ecology between these sedentary vs mobile grazers. 


###################################################################################
# MODELS OF GASTROPOD BIOMASS: SITE LEVEL (SIMPLE)                                #
###################################################################################

# Linear model (site level): SEAGRASS STRUCTURE
gastropod.lm.1 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      # + log10.Zostera.AG.mass.site 
                      , data = ZEN_2014_site_means_49)
summary(gastropod.lm.1)
#                                           Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                1.85428    1.48115   1.252   0.2169  
# log10.Zostera.longest.leaf.length.cm.site -1.01218    0.48920  -2.069   0.0442 *
# log10.Zostera.shoots.core.site             0.04419    0.34422   0.128   0.8984  
# 
# Residual standard error: 0.7704 on 46 degrees of freedom
# Multiple R-squared:  0.1266,	Adjusted R-squared:  0.0886 
# F-statistic: 3.333 on 2 and 46 DF,  p-value: 0.04449
AIC(gastropod.lm.1) # 118.3982


# Linear model (site level): PRODUCTIVITY
gastropod.lm.2 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                        log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(gastropod.lm.2)
#                                          Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                               -0.7874     0.4881  -1.613   0.1135  
# log10.periphyton.mass.per.g.zostera.site  -0.3742     0.1517  -2.467   0.0174 *
# Leaf.PercN.site                            0.3322     0.2081   1.596   0.1174  
# 
# Residual standard error: 0.765 on 46 degrees of freedom
# Multiple R-squared:  0.1388,	Adjusted R-squared:  0.1014 
# F-statistic: 3.708 on 2 and 46 DF,  p-value: 0.03215
AIC(gastropod.lm.2) # 117.7059


# Linear model (site level): SEAGRASS STRUCTURE + PRODUCTIVITY
gastropod.lm.3 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                        log10.Zostera.longest.leaf.length.cm.site 
                      + log10.Zostera.shoots.core.site
                      + log10.periphyton.mass.per.g.zostera.site 
                      + Leaf.PercN.site
                      , data = ZEN_2014_site_means_49)
summary(gastropod.lm.3)
#                                           Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                0.72074    1.65060   0.437   0.6645  
# log10.Zostera.longest.leaf.length.cm.site -0.90085    0.57454  -1.568   0.1241  
# log10.Zostera.shoots.core.site             0.03725    0.33713   0.110   0.9125  
# log10.periphyton.mass.per.g.zostera.site  -0.17872    0.18268  -0.978   0.3332  
# Leaf.PercN.site                            0.38602    0.20700   1.865   0.0689 .
# 
# Residual standard error: 0.7526 on 44 degrees of freedom
# Multiple R-squared:  0.2027,	Adjusted R-squared:  0.1303 
# F-statistic: 2.797 on 4 and 44 DF,  p-value: 0.03736
AIC(gastropod.lm.3) # 117.9269


# Linear model (site level): ABIOTIC ENVIRONMENT
gastropod.lm.4 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                        Temperature.C 
                      + Salinity.ppt
                      + Day.length.hours 
                      + log10.mean.fetch
                      , data = ZEN_2014_site_means_49)
summary(gastropod.lm.4)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      -2.04830    1.28638  -1.592   0.1185  
# Temperature.C     0.04518    0.02373   1.904   0.0635 .
# Salinity.ppt     -0.01916    0.01630  -1.175   0.2463  
# Day.length.hours  0.12393    0.05572   2.224   0.0313 *
#   log10.mean.fetch  0.13589    0.27059   0.502   0.6180  
# 
# Residual standard error: 0.757 on 44 degrees of freedom
# Multiple R-squared:  0.1933,	Adjusted R-squared:   0.12 
# F-statistic: 2.636 on 4 and 44 DF,  p-value: 0.04653
AIC(gastropod.lm.4) # 118.502


# Linear model (site level): ZOSTERA GENETICS
gastropod.lm.5 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                        log10.AllelicRichness 
                      , data = ZEN_2014_site_means_49)
summary(gastropod.lm.5)
#                       Estimate Std. Error t value Pr(>|t|)  
# (Intercept)             1.1817     0.4682   2.524   0.0150 *
#   log10.AllelicRichness  -1.4157     0.6961  -2.034   0.0477 *
# 
# Residual standard error: 0.7818 on 47 degrees of freedom
# Multiple R-squared:  0.08088,	Adjusted R-squared:  0.06132 
# F-statistic: 4.136 on 1 and 47 DF,  p-value: 0.04766
AIC(gastropod.lm.5) # 118.8969


# Linear model (site level): GEOGRAPHY
gastropod.lm.6 <- lm(log10.gastropod.mass.per.g.plant.site ~ 
                     Latitude
                     + Ocean
                     + Coast
                     , data = ZEN_2014_site_means_49)
summary(gastropod.lm.6)
#                    Estimate Std. Error t value Pr(>|t|)
# (Intercept)         0.01512    0.82426   0.018    0.985
# Latitude            0.01143    0.01499   0.762    0.450
# OceanPacific       -0.51004    0.44112  -1.156    0.254
# CoastEast.Pacific   0.27808    0.38105   0.730    0.469
# CoastWest.Atlantic -0.56003    0.35389  -1.583    0.121
# CoastWest.Pacific        NA         NA      NA       NA
# 
# Residual standard error: 0.7767 on 44 degrees of freedom
# Multiple R-squared:  0.1508,	Adjusted R-squared:  0.07355 
# F-statistic: 1.953 on 4 and 44 DF,  p-value: 0.1186
AIC(gastropod.lm.6) # 121.0225


# SUMMARY OF RESULTS: SITE level
# Little explanatory power. 'Best' model is productivity (AIC = 117.7059) with periphyton
# being significant predictor; however the effect of periphyton is negative, indicating 
# either top-down control by gastropods or that gastropods dominate at sites where periphyton 
# is low for other reasons. Very close second model is seagrass structure + productivity, 
# which returns a marginally non-significant positive effect of leaf % N (and, curiously, no 
# effect of periphyton). Surprisingly, the geography model showed nothing. Go figure. 


###################################################################################
# SEM: CRUSTACEAN BIOMASS                                                         #
###################################################################################

# Following models include geography, i.e., use only site (not Ocean or Coast) as 
# random nesting variables.  Note: Model can have either Ocean or Coast but not both,
# or it hangs up. 

# NOTE: the workaround below to incorporate simple lm for site-level varaibles doesn't 
# seem to work. The model appears ot be still using N = 1000 points (inctead of correct N = 50)
# as evidenced by highly significant P values of nearly all predictors for leaf %N. 

# Full model: Bottom-up

crust.list.1 <- list(
  
  crust.1 <- lme(log10.crustacean.mass.per.g.plant.imputed ~
    Ocean
    + Latitude
    # + Coast
    + Temperature.C
    + Salinity.ppt
    + log10.mean.fetch
    + log10.Leaf.PercN.imputed
    + AllelicRichness
    + log10.Zostera.proxy.production.rate
    + log10.Zostera.shoots.core.imputed
    + log10.Zostera.longest.leaf.length
    + log10.Zostera.AG.mass.imputed
    + log10.periphyton.mass.per.g.zostera.imputed
    , random = ~1 | Site
    , na.action = na.omit
    , data = ZEN_2014_master_data_49_imputed),
  
  gast.1 <- lme(log10.gastropod.mass.per.g.plant.imputed ~
                  Ocean
                + Latitude
                # + Coast
                + Temperature.C
                + Salinity.ppt
                + log10.mean.fetch
                + log10.Leaf.PercN.imputed
                + AllelicRichness
                + log10.Zostera.proxy.production.rate
                + log10.Zostera.shoots.core.imputed
                + log10.Zostera.longest.leaf.length
                + log10.Zostera.AG.mass.imputed
                + log10.periphyton.mass.per.g.zostera.imputed
                , random = ~1 | Site
                , na.action = na.omit
                , data = ZEN_2014_master_data_49_imputed),
  
  peri.1 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                  Ocean
                + Latitude
                # + Coast
                + Temperature.C
                + Salinity.ppt
                + log10.mean.fetch
                + log10.Leaf.PercN.imputed
                + AllelicRichness
                + log10.Zostera.shoots.core.imputed
                + log10.Zostera.longest.leaf.length
                + log10.Zostera.AG.mass.imputed
                + log10.Zostera.proxy.production.rate
                # + log10.mesograzer.mass.per.g.plant.imputed
                # + log10.periphyton.mass.per.g.zostera.imputed
                # + log10.crustacean.mass.per.g.plant.imputed
                # + log10.gastropod.mass.per.g.plant.imputed
                # + log10.grazer.richness.site # NOTE: left out because only relevant fo this resposne variable, and missing paths with others kill the chi square fit of SEM
                , random = ~1 | Site
                , na.action = na.omit
                , data = ZEN_2014_master_data_49_imputed),

  ZAG.1 <- lme(log10.Zostera.AG.mass.imputed ~ 
                 Ocean
               + Latitude
               # + Coast
               + Temperature.C
               + Salinity.ppt
               + log10.mean.fetch
               + log10.Leaf.PercN.imputed
               + AllelicRichness
               + log10.Zostera.shoots.core.imputed
               + log10.Zostera.longest.leaf.length
               + log10.Zostera.proxy.production.rate
               + pop.density.2015
               # + log10.Zostera.AG.mass.imputed
               # + log10.mesograzer.mass.per.g.plant.imputed
               # + log10.periphyton.mass.per.g.zostera.imputed
               # + log10.crustacean.mass.per.g.plant.imputed
               # + log10.gastropod.mass.per.g.plant.imputed
               # + log10.grazer.richness.site
               , random = ~1 | Site
               , na.action = na.omit
               , data = ZEN_2014_master_data_49_imputed),

  Zprodn.1 <- lme(log10.Zostera.proxy.production.rate ~ 
                    Ocean
                  + Latitude
                  # + Coast
                  + Temperature.C
               + Salinity.ppt
               + log10.mean.fetch
               + log10.Leaf.PercN.imputed
               + AllelicRichness
               + log10.Zostera.shoots.core.imputed
               + log10.Zostera.longest.leaf.length
               + pop.density.2015
               # + log10.Zostera.AG.mass.imputed
               # + log10.mesograzer.mass.per.g.plant.imputed
               # + log10.periphyton.mass.per.g.zostera.imputed
               # + log10.crustacean.mass.per.g.plant.imputed
               # + log10.gastropod.mass.per.g.plant.imputed
               # + log10.grazer.richness.site
               , random = ~1 | Site
               , na.action = na.omit
               , data = ZEN_2014_master_data_49_imputed),
  
  Zshoot.1 <- lme(log10.Zostera.shoots.core.imputed ~ 
                    Ocean
                  + Latitude
                  # + Coast
                  + Temperature.C
               + Salinity.ppt
               + log10.mean.fetch
               + log10.Leaf.PercN.imputed
               + AllelicRichness
               # + log10.Zostera.shoots.core.imputed
               # + log10.Zostera.longest.leaf.length
               + pop.density.2015
               # + log10.Zostera.AG.mass.imputed
               # + log10.mesograzer.mass.per.g.plant.imputed
               # + log10.periphyton.mass.per.g.zostera.imputed
               # + log10.crustacean.mass.per.g.plant.imputed
               # + log10.gastropod.mass.per.g.plant.imputed
               # + log10.grazer.richness.site
               , random = ~1 | Site
               , na.action = na.omit
               , data = ZEN_2014_master_data_49_imputed),
  
  leaf.1 <- lme(log10.Zostera.longest.leaf.length ~ # should this be lm since most predictors are site level?
                  Ocean
                + Latitude
                # + Coast
                + Temperature.C
                + Salinity.ppt
                + log10.mean.fetch
                + log10.Leaf.PercN.imputed
                + AllelicRichness
                # + log10.Zostera.shoots.core.imputed
                # + log10.Zostera.longest.leaf.length
                + pop.density.2015
                # + log10.Zostera.AG.mass.imputed
                # + log10.mesograzer.mass.per.g.plant.imputed
                # + log10.periphyton.mass.per.g.zostera.imputed
                # + log10.crustacean.mass.per.g.plant.imputed
                # + log10.gastropod.mass.per.g.plant.imputed
                # + log10.grazer.richness.site
                , random = ~1 | Site
                , na.action = na.omit
                , data = ZEN_2014_master_data_49_imputed),
  
  leafN.1 <- lm(log10.Leaf.PercN ~                 
                  # Ocean
                  + Latitude
                  # + Coast
                  + Temperature.C
                  + Salinity.ppt
                  + log10.mean.fetch
                  + AllelicRichness
                  + pop.density.2015
                , data = ddply(ZEN_2014_master_data_49_imputed, "Site", summarize, 
                               log10.Leaf.PercN = mean(log10.Leaf.PercN),
                        # Ocean = mean(Ocean),
                        Latitude = mean(Latitude),
                        Temperature.C = mean(Temperature.C),
                        Salinity.ppt = mean(Salinity.ppt),
                        log10.mean.fetch = mean(log10.mean.fetch),
                        AllelicRichness = mean(AllelicRichness),
                        pop.density.2015 = mean(pop.density.2015)))
)  
  

# Run goodness-of-fit tests
sem.fit(
  # List of structured equations
  modelList = crust.list.1, 
  data = ZEN_2014_master_data_49_imputed,
  # Additional variables who have no directed paths in the SEM but should be included in d-sep tests
  # add.vars = c("log10.Leaf.PercN.imputed"), 
  # Define the grouping variable
  grouping.vars = c("Site"),
  # Variables that should have correlated errors
  corr.errors = c("log10.Zostera.shoots.core.imputed ~~ log10.Zostera.longest.leaf.length", "log10.crustacean.mass.per.g.plant.imputed ~~ log10.gastropod.mass.per.g.plant.imputed")
  # , conditional = T # This displays all conditioning variables. Comment it out to suppress. 
)
#   fisher.c df p.value
# 1     1.13  6    0.98
# 
#      AIC    AICc  K   n
# 1 187.13 207.203 93 965

sem.coefs(crust.list.1, ZEN_2014_master_data_49_imputed, corr.errors = c("log10.Zostera.proxy.production.rate ~~ log10.Zostera.AG.mass.imputed"), standardize = "scale")
#                                       response                                   predictor     estimate  std.error p.value
# 9    log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.longest.leaf.length -0.174779743 0.05631161  0.0020
# 1    log10.crustacean.mass.per.g.plant.imputed                                OceanPacific  0.672269739 0.27072358  0.0171
# 11   log10.crustacean.mass.per.g.plant.imputed log10.periphyton.mass.per.g.zostera.imputed  0.102696950 0.04386693  0.0194
# 8    log10.crustacean.mass.per.g.plant.imputed           log10.Zostera.shoots.core.imputed -0.072596901 0.03468466  0.0366
# 4    log10.crustacean.mass.per.g.plant.imputed                                Salinity.ppt -0.278106270 0.14771815  0.0667
# 2    log10.crustacean.mass.per.g.plant.imputed                                    Latitude -0.287751145 0.15789127  0.0755
# 20 log10.periphyton.mass.per.g.zostera.imputed           log10.Zostera.longest.leaf.length  0.149803356 0.04832241  0.0020
# 22 log10.periphyton.mass.per.g.zostera.imputed         log10.Zostera.proxy.production.rate  0.135225434 0.05549985  0.0150
# 19 log10.periphyton.mass.per.g.zostera.imputed           log10.Zostera.shoots.core.imputed -0.106103431 0.06415681  0.0985
# 31               log10.Zostera.AG.mass.imputed           log10.Zostera.longest.leaf.length  0.370082921 0.06680652  0.0000
# 30               log10.Zostera.AG.mass.imputed           log10.Zostera.shoots.core.imputed  0.340066196 0.08933974  0.0002
# 24               log10.Zostera.AG.mass.imputed                                    Latitude -0.326378122 0.11323058  0.0063
# 23               log10.Zostera.AG.mass.imputed                                OceanPacific  0.344761776 0.18791724  0.0738
# 41         log10.Zostera.proxy.production.rate           log10.Zostera.shoots.core.imputed  1.049049619 0.01426202  0.0000
# 42         log10.Zostera.proxy.production.rate           log10.Zostera.longest.leaf.length  0.473723203 0.02296917  0.0000
# 35         log10.Zostera.proxy.production.rate                                    Latitude -0.272106188 0.06291065  0.0001
# 34         log10.Zostera.proxy.production.rate                                OceanPacific  0.245472839 0.10296976  0.0218
# 43         log10.Zostera.proxy.production.rate                            pop.density.2015 -0.107799315 0.04737086  0.0282
# 44           log10.Zostera.shoots.core.imputed                                OceanPacific -0.908303189 0.28375201  0.0026
# 47           log10.Zostera.shoots.core.imputed                                Salinity.ppt -0.341224576 0.16265110  0.0421
# 51           log10.Zostera.shoots.core.imputed                            pop.density.2015 -0.247805021 0.13657444  0.0769
# 52           log10.Zostera.longest.leaf.length                                OceanPacific  1.188650705 0.25189399  0.0000
# 60      ~~ log10.Zostera.proxy.production.rate            ~~ log10.Zostera.AG.mass.imputed  0.048746658         NA  0.0636


###################################################################################
# SEM: CRUSTACEAN BIOMASS                                                         #
###################################################################################

# Following models include geography, i.e., use only site (not Ocean or Coast) as 
# random nesting variables.  Note: Model can have either Ocean or Coast but not both,
# or it hangs up. 

# Full model: Bottom-up

gast.list.1 <- list(
  
  gast.1 <- lme(log10.gastropod.mass.per.g.plant.imputed ~
                   Ocean
                 + Latitude
                 # + Coast
                 + Temperature.C
                 + Salinity.ppt
                 + log10.mean.fetch
                 + log10.Leaf.PercN.imputed
                 + AllelicRichness
                 # + log10.Zostera.proxy.production.rate
                 + log10.Zostera.shoots.core.imputed
                 + log10.Zostera.longest.leaf.length
                 + log10.Zostera.AG.mass.imputed
                 + log10.periphyton.mass.per.g.zostera.imputed
                 , random = ~1 | Site
                 , na.action = na.omit
                 , data = ZEN_2014_master_data_49_imputed),
  
  peri.1 <- lme(log10.periphyton.mass.per.g.zostera.imputed ~ 
                  Ocean
                + Latitude
                # + Coast
                + Temperature.C
                + Salinity.ppt
                + log10.mean.fetch
                + log10.Leaf.PercN.imputed
                + AllelicRichness
                + log10.Zostera.shoots.core.imputed
                + log10.Zostera.longest.leaf.length
                + log10.Zostera.AG.mass.imputed
                + log10.Zostera.proxy.production.rate
                # + log10.mesograzer.mass.per.g.plant.imputed
                # + log10.periphyton.mass.per.g.zostera.imputed
                # + log10.crustacean.mass.per.g.plant.imputed
                # + log10.gastropod.mass.per.g.plant.imputed
                # + log10.grazer.richness.site # NOTE: left out because only relevant fo this resposne variable, and missing paths with others kill the chi square fit of SEM
                , random = ~1 | Site
                , na.action = na.omit
                , data = ZEN_2014_master_data_49_imputed),
  
  ZAG.1 <- lme(log10.Zostera.AG.mass.imputed ~ 
                 Ocean
               + Latitude
               # + Coast
               + Temperature.C
               + Salinity.ppt
               + log10.mean.fetch
               + log10.Leaf.PercN.imputed
               + AllelicRichness
               + log10.Zostera.shoots.core.imputed
               + log10.Zostera.longest.leaf.length
               + log10.Zostera.proxy.production.rate
               + pop.density.2015
               # + log10.Zostera.AG.mass.imputed
               # + log10.mesograzer.mass.per.g.plant.imputed
               # + log10.periphyton.mass.per.g.zostera.imputed
               # + log10.crustacean.mass.per.g.plant.imputed
               # + log10.gastropod.mass.per.g.plant.imputed
               # + log10.grazer.richness.site
               , random = ~1 | Site
               , na.action = na.omit
               , data = ZEN_2014_master_data_49_imputed),
  
  Zprodn.1 <- lme(log10.Zostera.proxy.production.rate ~ 
                    Ocean
                  + Latitude
                  # + Coast
                  + Temperature.C
                  + Salinity.ppt
                  + log10.mean.fetch
                  + log10.Leaf.PercN.imputed
                  + AllelicRichness
                  + log10.Zostera.shoots.core.imputed
                  + log10.Zostera.longest.leaf.length
                  + pop.density.2015
                  # + log10.Zostera.AG.mass.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + log10.grazer.richness.site
                  , random = ~1 | Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed),
  
  Zshoot.1 <- lme(log10.Zostera.shoots.core.imputed ~ 
                    Ocean
                  + Latitude
                  # + Coast
                  + Temperature.C
                  + Salinity.ppt
                  + log10.mean.fetch
                  + log10.Leaf.PercN.imputed
                  + AllelicRichness
                  # + log10.Zostera.shoots.core.imputed
                  # + log10.Zostera.longest.leaf.length
                  + pop.density.2015
                  # + log10.Zostera.AG.mass.imputed
                  # + log10.mesograzer.mass.per.g.plant.imputed
                  # + log10.periphyton.mass.per.g.zostera.imputed
                  # + log10.crustacean.mass.per.g.plant.imputed
                  # + log10.gastropod.mass.per.g.plant.imputed
                  # + log10.grazer.richness.site
                  , random = ~1 | Site
                  , na.action = na.omit
                  , data = ZEN_2014_master_data_49_imputed),
  
  leaf.1 <- lme(log10.Zostera.longest.leaf.length ~ # should this be lm since most predictors are site level?
                  Ocean
                + Latitude
                # + Coast
                + Temperature.C
                + Salinity.ppt
                + log10.mean.fetch
                + log10.Leaf.PercN.imputed
                + AllelicRichness
                # + log10.Zostera.shoots.core.imputed
                # + log10.Zostera.longest.leaf.length
                + pop.density.2015
                # + log10.Zostera.AG.mass.imputed
                # + log10.mesograzer.mass.per.g.plant.imputed
                # + log10.periphyton.mass.per.g.zostera.imputed
                # + log10.crustacean.mass.per.g.plant.imputed
                # + log10.gastropod.mass.per.g.plant.imputed
                # + log10.grazer.richness.site
                , random = ~1 | Site
                , na.action = na.omit
                , data = ZEN_2014_master_data_49_imputed)
)


# Run goodness-of-fit tests
sem.fit(
  # List of structured equations
  modelList = crust.list.1, 
  data = ZEN_2014_master_data_49_imputed,
  # Additional variables who have no directed paths in the SEM but should be included in d-sep tests
  # add.vars = c("log10.Leaf.PercN.imputed"), 
  # Define the grouping variable
  grouping.vars = c("Site"),
  # Variables that should have correlated errors
  corr.errors = c("log10.Zostera.shoots.core.imputed ~~ log10.Zostera.longest.leaf.length")
  # , conditional = T # This displays all conditioning variables. Comment it out to suppress. 
)
# fisher.c df p.value
# 
sem.coefs(gast.list.1, ZEN_2014_master_data_49_imputed, corr.errors = c("log10.Zostera.proxy.production.rate ~~ log10.Zostera.AG.mass.imputed"), standardize = "scale")



###################################################################################
# FOSSIL CODE                                                                     #
###################################################################################

# Dredge with linear models
## FOR EXPLORATORY PURPOSES ONLY!!! NEVER EVER DO THIS FOR A REAL 
## STUDY 
# fit model with all parameters 
test <- lm(log10.crustacean.mass.per.g.plant.site ~ Ocean + Coast + Basin + Latitude 
           + Longitude + Temperature.C + Salinity.ppt + glaciation + mean.fetch + pop.density.2015 
           + AllelicRichness + Zostera.AG.mass.site + macrophytes.total.AG.mass.core.site
           , data = ZEN_2014_site_means) 

test.nogeo <- lm(log10.crustacean.mass.per.g.plant.site ~ Temperature.C + Salinity.ppt + glaciation 
                 + mean.fetch + pop.density.2015 + AllelicRichness + Zostera.AG.mass.site + macrophytes.total.AG.mass.core.site
                 , data = ZEN_2014_site_means) 

# the dredge function fits all combinations of the variables in the model fit above 
options(na.action = "na.fail") 
results <- dredge(test) 
# subset(results, delta <5)
subset(results, delta == 0) 

crust.0 <- lm(log10.crustacean.mass.per.g.plant.site ~ 
                Longitude + Temperature.C + Salinity.ppt
              , data = ZEN_2014_site_means)
summary(crust.0)










