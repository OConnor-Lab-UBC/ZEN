###################################################################################
#                                                                                ##
# ZEN 2014 data analysis (JED) V1                                                ##
# Data are current as of 17 May 2016                                             ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2016-05-17                                                        ##
#                                                                                ##
###################################################################################

# TO DO:
#   on main data spreadsheet: harmonize variable names (column heads) in data tabs
#   with associated variable names on metadata tab. The  remove the fossil 'tab' 
#   headers in metadata tab. 

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA                                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SUMMARIZE MESOGRAZER RICHNESS AND OUTPUT SPECIES X PLOT MATRIX                  #
# ASSEMBLE DATA FRAMES INTO A MASTER DATA SET                                     #
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

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Read in summary SUBSITE data: 
SubsiteSummary <- read.csv(file = "ZEN_2014_Site&PlotData_2016_05_17_Released_subsite_summary.csv", header = TRUE)
names(SubsiteSummary)
str(SubsiteSummary)
# NOTE: Variables of interest from this file are primarily site-level metadata, plus 
# genetic data (which are site-level metrics). The many site-level means are recalucklated below.  

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
# CREATE MESOGRAZER RICHNESS SPECIES X PLOT MATRIX                                #
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
                            "Longitude", "Sampling.Time", "Month", "Date.Collected", "Temperature.C", "Salinity.ppt",
                            "Day.length.hours", "Depth.Categorical", "Depth.m", "Perc.Silt", "Perc.Sand", "Perc.Gravel", 
                            "GenotypicRichness", "AllelicRichness", "Inbreeding")]

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

hist(ZEN_2014_master_data$Zostera.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$total.macrophytes.aboveground.mean.mass, col = "cyan")    
hist(ZEN_2014_master_data$periphyton.mass.per.g.zostera, col = "cyan")    
hist(ZEN_2014_master_data$crustacean.abund.per.g.plant, col = "cyan")    

hist(ZEN_2014_master_data$gastropod.abund.per.g.plant, col = "cyan")    
hist(ZEN_2014_master_data$grazer.richness.plot, col = "cyan")    
hist(ZEN_2014_master_data$amphipod.survival, col = "cyan")    


###################################################################################
# LOG TRANSFORMS                                                                  #
###################################################################################

names(ZEN_2014_master_data)

# NOTE: For many variables I add a constant roughly equal to the smallest value recorded 
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
write.csv(ZEN_2014_master_data, "ZEN_2014_master_data_2016-05-17.csv", row.names = F)


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

# Add site-level environmental variables back in
ZEN_2014_site_means$Ocean <- ZEN_2014_master_data$Ocean[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Coast <- ZEN_2014_master_data$Coast[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Latitude <- ZEN_2014_master_data$Latitude[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Longitude <- ZEN_2014_master_data$Longitude[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Temperature.C <- ZEN_2014_master_data$Temperature.C[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Salinity.ppt <- ZEN_2014_master_data$Salinity.ppt[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Day.length.hours <- ZEN_2014_master_data$Day.length.hours[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Depth.Categorical <- ZEN_2014_master_data$Depth.Categorical[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Depth.m <- ZEN_2014_master_data$Depth.m[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Perc.Silt <- ZEN_2014_master_data$Perc.Silt[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Perc.Sand <- ZEN_2014_master_data$Perc.Sand[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Perc.Gravel <- ZEN_2014_master_data$Perc.Gravel[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$GenotypicRichness <- ZEN_2014_master_data$GenotypicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$AllelicRichness <- ZEN_2014_master_data$AllelicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$Inbreeding <- ZEN_2014_master_data$Inbreeding[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.GenotypicRichness <- ZEN_2014_master_data$log10.GenotypicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.AllelicRichness <- ZEN_2014_master_data$log10.AllelicRichness[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$grazer.richness.site <- ZEN_2014_master_data$grazer.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.grazer.richness.site <- ZEN_2014_master_data$log10.grazer.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]

ZEN_2014_site_means$grazer.estimated.richness.site <- ZEN_2014_master_data$grazer.estimated.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$grazer.estimated.Shannon.diversity.site <- ZEN_2014_master_data$grazer.estimated.Shannon.diversity.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$grazer.estimated.Simpson.diversity.site <- ZEN_2014_master_data$grazer.estimated.Simpson.diversity.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.grazer.estimated.richness.site <- ZEN_2014_master_data$log10.grazer.estimated.richness.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.grazer.estimated.Shannon.diversity.site <- ZEN_2014_master_data$log10.grazer.estimated.Shannon.diversity.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]
ZEN_2014_site_means$log10.grazer.estimated.Simpson.diversity.site <- ZEN_2014_master_data$log10.grazer.estimated.Simpson.diversity.site[match(ZEN_2014_site_means$Site, ZEN_2014_master_data$Site)]

names(ZEN_2014_site_means)

# Export the data frame
write.csv(ZEN_2014_site_means, "ZEN_2014_site_means_2016-05-18.csv", row.names = F)

# Create subset of site means for only Atlantic

# Create separate data sets for Atlantic and Pacific sites (subsite summary)
ZEN_2014_site_means.Atl <- droplevels(subset(ZEN_2014_site_means, Ocean == "Atlantic"))
ZEN_2014_site_means.Pac <- droplevels(subset(ZEN_2014_site_means, Ocean == "Pacific"))


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






