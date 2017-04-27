## validate soilgrids100 gg and pscs predictions using independent validation points
## code by Colby Brungard (cbrung@nmsu.edu)

# Set working directory and load necessary packages
setwd("D:/SoilGrids100")
library(plyr)
library(dplyr)
library(caret)
library(sp)
library(raster)
library(rgdal)
library(reshape)

load("./Validation Data/indepVal.Rdata")

# Validation of greatgroup predictions

# Some observations contain only partial data for greatgroup and particlesizeclass. 
# Split into seperate validation datasets and subset to remove missing observations
gg.val <- val[, -5]
gg.val2 <- na.omit(gg.val) # 2063 validation observations

# how many observations by area?  
ddply(gg.val2, .(area), .fun = nrow)

# area  V1
# 1    BeaverCountyUT 326
# 2     BeaverCreekWY  51
# 3  BoundaryWatersMN 173
# 4             CRCUT 361
# 5       FortBlissNM 103
# 6              Iowa  69
# 7      JuabCountyUT 151
# 8             Maine 193
# 9   MillardCountyUT 604
# 10    SnakeValleyUT  32

# Check for classes in validation data that are not in the soilgrids100 legend.
# Define function to enforce consistent capitilization
proper = function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

# Read in legend, select the column of names, convert to dataframe and use anti_join function to compare levels not in each dataset
gg.levs <- read.csv("./Predictions/TAXgg_M_100m.tif.csv")
gg.levs$Class <- proper(gg.levs$Class)
sg.gg <- data.frame(gg.levs$Class) # must be a dataframe for anti_join to work
names(sg.gg) <- 'greatgroup'
vl.gg <- data.frame(gg.val2$greatgroup) # must be a dataframe for anti_join to work
names(vl.gg) <- 'greatgroup'
unique(anti_join(vl.gg, sg.gg)) # list the factor levels in the validation data that do not match the gg legend 

# Rename classes with no counterpart in the soilgrids legend to match those in the soilgrids legends. Mismatches were mostly misspelled words.   
gg.val2$greatgroup <- revalue(gg.val2$greatgroup, c(
  'Aquicambid' = 'Aquicambids',
  'Calciargid' = 'Calciargids', 
  'Endoaquent' = 'Endoaquents', 
  'Halargids' = 'Haplargids',
  'Haplargid' = 'Haplargids',
  'Haplocalcid' = 'Haplocalcids',
  'Haplocambid' = 'Haplocambids',
  'Hapudalfs' = 'Hapludalfs',
  'Natrargid' = 'Natrargids', 
  'Torrirothents' = 'Torriorthents'))

# Convert to spatial points dataframe and Write to shapefile for travis
coordinates(gg.val2) <- ~ longitude+latitude
proj4string(gg.val2) <- '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs' 

# writeOGR(gg.val2, "./Validation Data/shapefiles", "greatgroups", driver="ESRI Shapefile", overwrite = TRUE)

# Load soilgrids predictions
gg <- raster("./Predictions/TAXgg_M_100m.tif")

# Extract predictions at each validation point then rename the extracted values to enable joining
ggE <- extract(gg, gg.val2, sp = TRUE)
names(ggE@data)[4] <- 'Value'

# There are some pedons < 250m apart. I'm not going to worry about this.

# Convert extracted (i.e., predicted) numbers to class names in the legend
gg.pred <- join(ggE@data, gg.levs, by = 'Value', type = 'left')
gg.pred$Class <- as.factor(gg.pred$Class)

# Calculate confusion matricies
# Inorder to calculate meaningful confusion matricies I had to match the existing factor levels between predicted and observed as not all classes were observed in the validation datasets then drop all factor levels that were not used. 
# Drop all empty levels, get levles common between datasets, then force each factor to have the same levels
gg.notdrop <- unique(c(levels(droplevels(gg.pred$greatgroup)), levels(droplevels(gg.pred$Class))))
gg.pred$greatgroup <- factor(gg.pred$greatgroup, levels = gg.notdrop)
gg.pred$Class <- factor(gg.pred$Class, levels = gg.notdrop)

# All validation points 
gg.cm.all <- confusionMatrix(gg.pred$Class, gg.pred$greatgroup)
print(gg.cm.all$overall)
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 3.869223e-01   2.922213e-01   3.656862e-01   4.084850e-01   3.092429e-01   6.341386e-14            NaN

# By individual validation region. 
cm.gg <- function(x){ 
  notdrop <- unique(c(levels(droplevels(x$greatgroup)), levels(droplevels(x$Class))))
  x$greatgroup <- factor(x$greatgroup, levels = notdrop)
  x$Class <- factor(x$Class, levels = notdrop)
  CM <- confusionMatrix(x$Class, x$greatgroup)
  return(CM)
}

gg.cm.area <- dlply(gg.pred, .(area), .fun = cm.gg)

# Make nice plots
gg.overall <- ldply(.data = gg.cm.area, function(x) x$overall)

# Rename each area to be more meaningful
gg.overall$area <- mapvalues(gg.overall$area, c('BeaverCountyUT','BeaverCreekWY','BoundaryWatersMN','CRCUT','FortBlissNM','JuabCountyUT','MillardCountyUT','SnakeValleyUT', 'Iowa'), c('Western UT (Beaver County)','Northcentral WY','Northern MN','Southeastern UT','Southcentral NM','Western UT (Juab County)','Western UT (Millard County)','Western UT (Snake Valley)', 'Eastern Iowa'))
 
# Great Group Plots
ggplot(gg.overall, aes(x=area, y = Accuracy)) + geom_point() + 
  xlab('') + coord_flip() + theme_bw() + scale_x_discrete(limits=gg.overall[order(gg.overall$Accuracy), ]$area)
ggsave('./Validation Data/ggOverall.png', height = 7, width = 7, units = 'in')

ggplot(gg.overall, aes(x=area, y = Kappa)) + geom_point() + 
  xlab('') + coord_flip() + theme_bw() + scale_x_discrete(limits=gg.overall[order(gg.overall$Kappa), ]$area)
ggsave('./Validation Data/ggKappa.png', height = 7, width = 7, units = 'in')


# ------------------------------------------------------------------------------------------------
# Validation of particle size class predictions

# Omit missing particle size class observations
psc.val <- val[, -4]
psc.val2 <- na.omit(psc.val) 
dim(psc.val2) # 2095 validation observations

# Travis changed several of the Psaments to be lithic psamments. Rather than reading this new data in it is easier to just modify the existing dataset. 
levels(psc.val2$particlesizeclass) <- c(levels(psc.val2$particlesizeclass), 'Lithic-Psamments') #Add a factor level
psc.val2[psc.val2$id == 'C293', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C313', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D019', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D024', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D025', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D035', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D196', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'D199', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C040', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C049', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C205', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C221', ]$particlesizeclass <- 'Lithic-Psamments'
psc.val2[psc.val2$id == 'C249', ]$particlesizeclass <- 'Lithic-Psamments'

# how many observations by area?  
ddply(psc.val2, .(area), .fun = nrow)

# area  V1
# 1    BeaverCountyUT 326
# 2     BeaverCreekWY  51
# 3  BoundaryWatersMN 161
# 4             CRCUT 356
# 5       FortBlissNM 103
# 6              Iowa  69
# 7      JuabCountyUT 204
# 8             Maine 187
# 9   MillardCountyUT 606
# 10    SnakeValleyUT  32

# Rename classes with no counterpart in the soilgrids legend to match those in the soilgrids legends. For greatgroup mismatches were mostly misspelled words.   
psc.levs <- read.table("./Predictions/pscs_lookup_20161227_ttab.txt", sep = '\t', header = TRUE)
sg.psc <- data.frame(psc.levs[,2]) 
names(sg.psc) <- 'particlesizeclass'
vl.psc <- data.frame(toupper(psc.val2$particlesizeclass)) # make upper case to match legend 
names(vl.psc) <- 'particlesizeclass'
unique(anti_join(vl.psc, sg.psc)) # list the factor levels in the validation data that do not match the gg legend

# particlesizeclass mismatches were because soilgrids lumped X over sandy or sandy-skeletal. Use toupper to make sure captilization matches
psc.val2$particlesizeclass <- revalue(toupper(psc.val2$particlesizeclass), c(
  'FINE-SILTY OVER SANDY-SKELETAL' = 'FINE-SILTY OVER SANDY OR SANDY-SKELETAL',
  'LOAMY OVER SANDY-SKELETAL' = 'LOAMY OVER SANDY OR SANDY-SKELETAL',
  'LOAMY-SKELETAL OVER SANDY' = 'LOAMY-SKELETAL OVER SANDY OR SANDY-SKELETAL',
  'COARSE-LOAMY OVER SANDY-SKELETAL' = 'COARSE-LOAMY OVER SANDY OR SANDY-SKELETAL',
  'FINE-LOAMY OVER SANDY' = 'FINE-LOAMY OVER SANDY OR SANDY-SKELETAL',
  'FINE-SILTY OVER SANDY' = 'FINE-SILTY OVER SANDY OR SANDY-SKELETAL'))

# Convert to spatial points dataframe and Write to shapefile for travis
coordinates(psc.val2) <- ~ longitude+latitude
proj4string(psc.val2) <- '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs' 

# writeOGR(psc.val2, "./Validation Data/shapefiles", "particlesizeclass", driver="ESRI Shapefile", overwrite = TRUE)

# Load soilgrids predictions. This is Travis' hardened class map. https://dl.dropboxusercontent.com/u/99205734/pscs_hardened_20161227.zip
psc <- raster("./Predictions/pscs_c1_20161227comp.tif") 

# Extract predictions at each validation point then rename the extracted values to enable joining
pscE <- extract(psc, psc.val2, sp = TRUE) 
names(pscE@data)[4] <- 'Value'

# There are some pedons < 250m apart. I'm not going to worry about this.

# Convert extracted (i.e., predicted) numbers to class names in the legend
psc.pred <- join(pscE@data, psc.levs, by = 'Value', type = 'left')
psc.pred$particlesizeclass <- as.factor(psc.pred$particlesizeclass)

# Calculate confusion matricies
# Inorder to calculate meaningful confusion matricies I had to match the existing factor levels between predicted and observed, 
# then drop all factor levels that were not used. 
# Drop all empty levels, get levles common between datasets, then force each factor to have the same levels.
# Note: particlesizeclass is observed, PSCS is predicted
psc.notdrop <- unique(c(levels(droplevels(psc.pred$particlesizeclass)), levels(droplevels(psc.pred$PSCS))))
psc.pred$particlesizeclass <- factor(psc.pred$particlesizeclass, levels = psc.notdrop)
psc.pred$PSCS <- factor(psc.pred$PSCS, levels = psc.notdrop)

# All validation points 
psc.cm.all <- confusionMatrix(psc.pred$PSCS, psc.pred$particlesizeclass)
print(psc.cm.all$overall)
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 4.893514e-01   3.442235e-01   4.675822e-01   5.111509e-01   2.938045e-01   1.092028e-77            NaN 

# By individual validation region. 
cm.psc <- function(x){ 
  notdrop <- unique(c(levels(droplevels(x$particlesizeclass)), levels(droplevels(x$PSCS))))
  x$particlesizeclass <- factor(x$particlesizeclass, levels = notdrop)
  x$PSCS <- factor(x$PSCS, levels = notdrop)
  CM <- confusionMatrix(x$PSCS, x$particlesizeclass)
  return(CM)
}

psc.cm.area <- dlply(psc.pred, .(area), .fun = cm.psc)

# Make nice plots
psc.overall <- ldply(.data = psc.cm.area, function(x) x$overall)

# Rename each area to be more meaningful
psc.overall$area <- mapvalues(gg.overall$area, c('BeaverCountyUT','BeaverCreekWY','BoundaryWatersMN','CRCUT','FortBlissNM','JuabCountyUT','MillardCountyUT','SnakeValleyUT', 'Iowa'), c('Western UT (Beaver County)','Northcentral WY','Northern MN','Southeastern UT','Southcentral NM','Western UT (Juab County)','Western UT (Millard County)','Western UT (Snake Valley)', 'Eastern Iowa'))

# Particle Size Class Plots
ggplot(psc.overall, aes(x=area, y = Accuracy)) + geom_point() + 
  xlab('') + coord_flip() + theme_bw() + scale_x_discrete(limits=psc.overall[order(psc.overall$Accuracy), ]$area)
ggsave('./Validation Data/pscOverall.png', height = 7, width = 7, units = 'in')

ggplot(psc.overall, aes(x=area, y = Kappa)) + geom_point() + 
  xlab('') + coord_flip() + theme_bw() + scale_x_discrete(limits=psc.overall[order(psc.overall$Kappa), ]$area)
ggsave('./Validation Data/pscKappa.png', height = 7, width = 7, units = 'in')

# end
