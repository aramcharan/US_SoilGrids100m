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
library(geosphere) 

validationData <- read.csv("./Validation Data/Travis/Valid_SoilGrids_201704_GG_pscsmodorg_cwb.csv") 

# Rename each area to be more meaningful
validationData$area <- mapvalues(validationData$area, c('BeaverCountyUT','BeaverCreekWY','BoundaryWatersMN','CRCUT','FortBlissNM','JuabCountyUT','MillardCountyUT', 'Iowa'), c('Western UT (Beaver County)','Northcentral WY','Northern MN','Southeastern UT','Southcentral NM','Western UT (Juab County)','Western UT (Millard County)', 'Eastern Iowa'))

#Reorder the classes for consistency 
validationData$area <- factor(validationData$area, levels = c('Eastern Iowa',  'Maine','Northcentral WY', 'Northern MN','Southcentral NM','Southeastern UT','Western UT (Beaver County)','Western UT (Juab County)','Western UT (Millard County)'))


# 1. Validation of greatgroup predictions
# Subset data for only greatgroup observations and remove missing values
gg.val1 <- validationData[,c(1:4,6)]
gg.val <- na.omit(gg.val1) # 1999 greatgroup validation observations

# 1.1 how many observations by area?  
gg.nobsPerArea <- ddply(gg.val, .(area), .fun = nrow)
gg.nobsPerArea

#                          area  V1
# 1                Eastern Iowa  69
# 2                       Maine 183
# 3             Northcentral WY  51
# 4                 Northern MN  97
# 5             Southcentral NM 103
# 6             Southeastern UT 361
# 7  Western UT (Beaver County) 326
# 8    Western UT (Juab County) 204
# 9 Western UT (Millard County) 605

# 1.2 Check for classes in validation data that are not in the soilgrids100 legend.
# Define function to enforce consistent capitilization
proper = function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

# validation observations. Convert to dataframe and use anti_join function to compare levels not in each dataset
gg.val$greatgroup <- as.factor(proper(gg.val$greatgroup))
vl.gg <- as.data.frame(gg.val$greatgroup)
names(vl.gg) <- 'greatgroup'

# Read in legend, enforce capitilization, select the column of names, and match column name to validation dataset column name 
gg.levs <- read.csv("./Predictions/TAXgg_M_100m.tif.csv")
gg.levs$Class <- proper(gg.levs$Class)
sg.gg <- data.frame(gg.levs$Class) # must be a dataframe for anti_join to work
names(sg.gg) <- 'greatgroup'

# list the factor levels in the validation data that do not match the gg legend
unique(anti_join(vl.gg, sg.gg)) #all factor levels match (use revalue to change factor levels if needed for missspelled words, or change in data directly). 


# 1.3 Number of greatgroup classes per area
ggPerArea <- gg.val %>% 
              group_by(area) %>%
               summarise(n_gg = length(unique(greatgroup)))
ggPerArea

#                        area  n_gg
#                Eastern Iowa     4
#                       Maine    12
#             Northcentral WY     5
#                 Northern MN    13
#             Southcentral NM     5
#             Southeastern UT     9
#  Western UT (Beaver County)    11
#    Western UT (Juab County)     8
# Western UT (Millard County)    12


# 1.4 Calculate observation density obs/km2 per area.
# Code modified from: https://chitchatr.wordpress.com/2015/01/23/calculating-the-area-of-a-convex-hull/
# and https://www.rdocumentation.org/packages/geosphere/versions/1.5-5/topics/areaPolygon

# Split gg.val into a list by area
gglist <- split(gg.val[,c(2,3)], f = gg.val$area)

# Calculate convex hull around each set of points and get area in km2
ggarealist <- list()
for(i in seq(gglist)){
 test <- chull(gglist[[i]]) # convex hull
 test2 <- gglist[[i]][test,] #Extract the coordinates of the convex hull
 ggarealist[[i]] <- areaPolygon(test2[,c(2,1)])/1e6 #areaPolygon returns m2 so divide by 1e6 to convert to km2
 }

# Calculate point density ( Number of points/area km2)
ggpointdens <- data.frame()
for(i in seq(gglist)){
 ggpointdens[i,1] <- ggPerArea[i,2]/ggarealist[[i]]
 }
names(ggpointdens) <- 'PointDensity'
print(cbind(ggPerArea$area,ggpointdens))

#              ggPerArea$area PointDensity
#                Eastern Iowa 2.5308259330
#                       Maine 0.0005307507
#             Northcentral WY 0.0205269708
#                 Northern MN 0.0022483874
#             Southcentral NM 0.0262662896
#             Southeastern UT 0.0064692880
#  Western UT (Beaver County) 0.0341146861
#    Western UT (Juab County) 0.0099128066
# Western UT (Millard County) 0.0239207803


# 1.5 Convert to spatial points dataframe and Write to shapefile for travis
coordinates(gg.val) <- ~ longitude+latitude
proj4string(gg.val) <- '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs'

# writeOGR(gg.val2, "./Validation Data/shapefiles", "greatgroups", driver="ESRI Shapefile", overwrite = TRUE)


# 1.6 Load soilgrids predictions
gg <- raster("./Predictions/TAXgg_M_100m.tif")

# Extract predictions at each validation point then rename the extracted valuesvl.gg$greatgroup <- proper(vl.gg$greatgroup) to enable joining
ggE <- extract(gg, gg.val, sp = TRUE)
names(ggE@data)[4] <- 'Value'

# There are some pedons < 250m apart. I'm not going to worry about this.

# Convert extracted (i.e., predicted) numbers to class names in the legend
gg.pred <- join(ggE@data, gg.levs, by = 'Value', type = 'left')
gg.pred$Class <- as.factor(gg.pred$Class)


# 1.7 Calculate confusion matricies
# Inorder to calculate meaningful confusion matricies I had to match the existing factor levels between predicted and observed as not all classes were observed in the validation datasets then drop all factor levels that were not used. 
# Drop all empty levels, get levles common between datasets, then force each factor to have the same levels
gg.notdrop <- unique(c(levels(droplevels(gg.pred$greatgroup)), levels(droplevels(gg.pred$Class))))
gg.pred$greatgroup <- factor(gg.pred$greatgroup, levels = gg.notdrop)
gg.pred$Class <- factor(gg.pred$Class, levels = gg.notdrop)

# All validation points 
gg.cm.all <- confusionMatrix(gg.pred$Class, gg.pred$greatgroup)
print(gg.cm.all$overall)
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 3.677873e-01   2.623293e-01   3.465738e-01   3.893905e-01   3.266433e-01   5.744117e-05            NaN

# By individual validation region. 
cm.gg <- function(x){ 
  notdrop <- unique(c(levels(droplevels(x$greatgroup)), levels(droplevels(x$Class))))
  x$greatgroup <- factor(x$greatgroup, levels = notdrop)
  x$Class <- factor(x$Class, levels = notdrop)
  CM <- confusionMatrix(x$Class, x$greatgroup)
  return(CM)
}

gg.cm.area <- dlply(gg.pred, .(area), .fun = cm.gg)
save(gg.cm.area, file = "./Validation Data/SoilGrids GreatGroup Confusion Matricies.RData")

# 1.9 Make nice plots. Reverse factor levels to make a nicer plot. 
gg.overall <- ldply(.data = gg.cm.area, function(x) x$overall)
levels(gg.overall$area) <- rev(levels(gg.overall$area))

# Add Number of great group classes per area
gg.overall$ggPerArea <- ggPerArea[,2]

# Add number of observation number per area
gg.overall$obsPerArea <- gg.nobsPerArea[,2]

# add point density in points/km2
gg.overall$pointDens <- ggpointdens[,1]

# Great Group Plots (can't use theme_bw() as it doesn't center the title.)
ggplot(gg.overall, aes(x=Accuracy, y = rev(area))) + geom_point(aes(size = ggPerArea, color = obsPerArea)) + 
  scale_color_gradient(low = 'blue', high = 'red', name = 'Number of \nObservations') + 
  scale_size_area(name = 'Number of\nClasses', breaks = c(4,6,8,10,12,14))  + 
  xlim(0,1) + 
  ylab('') + 
  ggtitle("Great Groups") + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = 'white', colour = 'black', linetype = 'solid'),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey90"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))


# ggsave('./Validation Data/ggOverall 042017.png', height = 4, width = 7, units = 'in')

# 1.8 Distribution of observations in each greatgroup by area 

gglist2 <- split(gg.val, f = gg.val$area)
barchart(droplevels(gglist2[[1]]$greatgroup), main = names(gglist2)[[1]], col = 'gray')
barchart(droplevels(gglist2[[2]]$greatgroup), main = names(gglist2)[[2]], col = 'gray')
barchart(droplevels(gglist2[[3]]$greatgroup), main = names(gglist2)[[3]], col = 'gray')
barchart(droplevels(gglist2[[4]]$greatgroup), main = names(gglist2)[[4]], col = 'gray')
barchart(droplevels(gglist2[[5]]$greatgroup), main = names(gglist2)[[5]], col = 'gray')
barchart(droplevels(gglist2[[6]]$greatgroup), main = names(gglist2)[[6]], col = 'gray')
barchart(droplevels(gglist2[[7]]$greatgroup), main = names(gglist2)[[7]], col = 'gray')
barchart(droplevels(gglist2[[8]]$greatgroup), main = names(gglist2)[[8]], col = 'gray')
barchart(droplevels(gglist2[[9]]$greatgroup), main = names(gglist2)[[9]], col = 'gray')

# 1.9 Thoughts
# These graphs show that most of the validation datasets are dominated by one or two greatgroup classes (one or two classes dominate at all taxonomic levels in my experience). 
# There doesn't seem to be a direct relationship with overall validation accuracy (compare MN-high accuracy (plot 3) with Beaver County UT-low accuracy (plot1)). This suggests to me that the accuracy of class predictions is mostly dominated by how well the model predicts the majority class. 
# Since predictive accuracy doesn't seem to be related to the number of classes, the number of observations, or observation density, and the models used are highly flexible and usally find relationships if they exist, I conclude that model accuracy is mostly dependent upon the strength of the covariate-soil relationships. In some landscapes, the covariates better represent soil distribution than others. This suggests that need for better covariates.  

# ------------------------------------------------------------------------------------------------
# 2. Validation of particle size class predictions

# Omit missing particle size class observations
psc.val1 <- validationData[,c(1:3,5:6)]
psc.val <- na.omit(psc.val1) # 2012 validation observations

# 2.1 how many observations by area?  
psc.nobsPerArea <- ddply(psc.val, .(area), .fun = nrow)
psc.nobsPerArea

 #                      area  V1
 #               Eastern Iowa  69
 #                      Maine 183
 #            Northcentral WY  51
 #                Northern MN 108
 #            Southcentral NM 103
 #            Southeastern UT 356
 # Western UT (Beaver County) 329
 #   Western UT (Juab County) 208
 # Western UT (Millard County) 605

# 2.2 Check for classes in validation data that are not in the soilgrids100 legend.

# validation observations. Convert to dataframe and use anti_join function to compare levels not in each dataset
psc.val$particlesizeclass <- as.factor(proper(psc.val$particlesizeclass))
vl.psc <- as.data.frame(psc.val$particlesizeclass)
names(vl.psc) <- 'particlesizeclass'

# Read in legend, enforce capitilization, select the column of names, and match column name to validation dataset column name 
psc.levs <- read.csv("./Predictions/PSCS_M_100m.tif.csv")
psc.levs$Class <- proper(psc.levs$Class)
sg.psc <- data.frame(psc.levs$Class) # must be a dataframe for anti_join to work
names(sg.psc) <- 'particlesizeclass'

# list the factor levels in the validation data that do not match the psc legend.
unique(anti_join(vl.psc, sg.psc)) 

psc.val$particlesizeclass <-revalue(psc.val$particlesizeclass, 
       c("Loamy over sandy-skeletal" = "Loamy over sandy or sandy-skeletal",
           "Fine-loamy over sandy" = "Fine-loamy over sandy or sandy-skeletal",
       "Loamy-skeletal over sandy" = "Loamy-skeletal over sandy or sandy-skeletal",
"Coarse-loamy over sandy-skeletal" = "Coarse-loamy over sandy or sandy-skeletal"))


# 2.3 Number of greatgroup classes per area
pscPerArea <- psc.val %>% 
  group_by(area) %>%
  summarise(n_psc = length(unique(particlesizeclass)))
pscPerArea

#                          area n_psc
# 1                Eastern Iowa     3
# 2                       Maine    14
# 3             Northcentral WY     7
# 4                 Northern MN    12
# 5             Southcentral NM     8
# 6             Southeastern UT     7
# 7  Western UT (Beaver County)    13
# 8    Western UT (Juab County)    16
# 9 Western UT (Millard County)    10


# 2.4 Calculate observation density obs/km2 per area.

# Split psc.val into a list by area
psclist <- split(psc.val[,c(2,3)], f = psc.val$area)

# Calculate convex hull around each set of points and get area in km2
pscarealist <- list()
for(i in seq(psclist)){
  test <- chull(psclist[[i]]) # convex hull
  test2 <- psclist[[i]][test,] #Extract the coordinates of the convex hull
  pscarealist[[i]] <- areaPolygon(test2[,c(2,1)])/1e6 #areaPolygon returns m2 so divide by 1e6 to convert to km2
}

# Calculate point density ( Number of points/area km2)
pscpointdens <- data.frame()
for(i in seq(psclist)){
  pscpointdens[i,1] <- pscPerArea[i,2]/pscarealist[[i]]
}
names(pscpointdens) <- 'PointDensity'
print(cbind(pscPerArea$area,pscpointdens))

#              pscPerArea$area PointDensity
#                 Eastern Iowa 1.8981194497
#                        Maine 0.0006192091
#              Northcentral WY 0.0287377591
#                  Northern MN 0.0020684061
#              Southcentral NM 0.0420260633
#              Southeastern UT 0.0050966034
#   Western UT (Beaver County) 0.0403173563
#     Western UT (Juab County) 0.0196421218
#  Western UT (Millard County) 0.0199339836


# 2.5 Convert to spatial points dataframe and Write to shapefile
coordinates(psc.val) <- ~ longitude+latitude
proj4string(psc.val) <- '+proj=longlat +ellps=GRS80 +datum=WGS84 +no_defs'

# writeOGR(gg.val2, "./Validation Data/shapefiles", "greatgroups", driver="ESRI Shapefile", overwrite = TRUE)


# 2.6 Load soilgrids predictions
psc <- raster("./Predictions/PSCS_M_100m.tif")

# Extract predictions at each validation point then rename the extracted valuesvl.gg$greatgroup <- proper(vl.gg$greatgroup) to enable joining
pscE <- extract(psc, psc.val, sp = TRUE)
names(pscE@data)[4] <- 'Value'

# There are some pedons < 250m apart. I'm not going to worry about this.

# Convert extracted (i.e., predicted) numbers to class names in the legend
psc.pred <- join(pscE@data, psc.levs, by = 'Value', type = 'left')
psc.pred$Class <- as.factor(psc.pred$Class)


# 2.6 Calculate confusion matricies
# Inorder to calculate meaningful confusion matricies I had to match the existing factor levels between predicted and observed as not all classes were observed in the validation datasets then drop all factor levels that were not used. 
# Drop all empty levels, get levles common between datasets, then force each factor to have the same levels
psc.notdrop <- unique(c(levels(droplevels(psc.pred$particlesizeclass)), levels(droplevels(psc.pred$Class))))
psc.pred$particlesizeclass <- factor(psc.pred$particlesizeclass, levels = psc.notdrop)
psc.pred$Class <- factor(psc.pred$Class, levels = psc.notdrop)

# All validation points 
psc.cm.all <- confusionMatrix(psc.pred$Class, psc.pred$particlesizeclass)
print(psc.cm.all$overall)
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 4.146707e-01   2.830211e-01   3.929893e-01   4.366021e-01   2.460080e-01   1.208664e-61            NaN

# By individual validation region. 
cm.psc <- function(x){ 
  notdrop <- unique(c(levels(droplevels(x$particlesizeclass)), levels(droplevels(x$Class))))
  x$particlesizeclass <- factor(x$particlesizeclass, levels = notdrop)
  x$Class <- factor(x$Class, levels = notdrop)
  CM <- confusionMatrix(x$Class, x$particlesizeclass)
  return(CM)
}

psc.cm.area <- dlply(psc.pred, .(area), .fun = cm.psc)
save(psc.cm.area, file = "./Validation Data/SoilGrids Particle Size Class Confusion Matricies.RData")


# 2.7 Make nice plots
psc.overall <- ldply(.data = psc.cm.area, function(x) x$overall)
levels(psc.overall$area) <- rev(levels(psc.overall$area))

# Add Number of great group classes per area
psc.overall$pscPerArea <- pscPerArea[,2]

# Add number of observation number per area
psc.overall$obsPerArea <- psc.nobsPerArea[,2]

# add point density in points/km2
psc.overall$pointDens <- pscpointdens[,1]

# Particle Size Class Plots
ggplot(psc.overall, aes(x=Accuracy, y = rev(area))) + geom_point(aes(size = pscPerArea, color = obsPerArea)) + 
  scale_color_gradient(low = 'blue', high = 'red', name = 'Number of \nObservations') + 
  scale_size_area(name = 'Number of\nClasses', breaks = c(4,6,8,10,12,14))  + 
  xlim(0,1) + 
  ylab('') + 
  ggtitle("Modified Particle Size Class") + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = 'white', colour = 'black', linetype = 'solid'),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey90"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))

# ggsave('./Validation Data/pscOverall 042017.png', height = 4, width = 7, units = 'in')

# 2.8 Distribution of observations in each greatgroup by area 
psclist2 <- split(psc.val, f = psc.val$area)
barchart(droplevels(psclist2[[1]]$particlesizeclass), main = names(psclist2)[[1]], col = 'gray')
barchart(droplevels(psclist2[[2]]$particlesizeclass), main = names(psclist2)[[2]], col = 'gray')
barchart(droplevels(psclist2[[3]]$particlesizeclass), main = names(psclist2)[[3]], col = 'gray')
barchart(droplevels(psclist2[[4]]$particlesizeclass), main = names(psclist2)[[4]], col = 'gray')
barchart(droplevels(psclist2[[5]]$particlesizeclass), main = names(psclist2)[[5]], col = 'gray')
barchart(droplevels(psclist2[[6]]$particlesizeclass), main = names(psclist2)[[6]], col = 'gray')
barchart(droplevels(psclist2[[7]]$particlesizeclass), main = names(psclist2)[[7]], col = 'gray')
barchart(droplevels(psclist2[[8]]$particlesizeclass), main = names(psclist2)[[8]], col = 'gray')
barchart(droplevels(psclist2[[9]]$particlesizeclass), main = names(psclist2)[[9]], col = 'gray')

# 2.9 Thoughts
#These plots show a similar lack of any trend between observation number, number of classes, point density, data imbalance, and accuracy. I must conclude that differences in accuracy between validation datasets is based on the strength of the relationships between covariates and soil classes. 

# No.... actually its about the quality of the predictions. Predictions in some areas are better than others. Why? The model is the same, so it must be either that 1) covariate-class relationships are better in some areas than others, or 2) there was more training data in some areas than other areas.  The above plots show that it's not about differences in validation datasets. I could answer #2 by calculating the training points within each area, of course how do I deal with similar areas. 
# Predictiosn are better because the model is better. The model would be different because covariate-class relationships are stronger. Tthis could be because the covariate better represent the distribution of classes in some areas or that there was more training training data in similar areas so the model is stronger. 

# the lack of a clear trend between the number of obsevations, the number of classes, observation density, and overall accuracy indicates that differences in validation accuracy between areas is not a result of the independent validation data, but of actual differences in prediction accuracy between areas. Differences in prediction quality between areas is lkely a result of differences in the strength of covariate-class relationships which is driving by the covariate-class relationships in each area and the availability of training data in areas similar to each validation area.  

#The Maine greatgroup predictions were so good because the model accuractly predicted the dominant class (Dystrudepts? - supplementary materials). 

# end
