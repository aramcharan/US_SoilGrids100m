## Plots USA48 SoilGrids for model fitting results / predictions
## by: tom.hengl@isric.org and amanda.ramcharan@gmail.com

setwd("/data/NASIS")
library(caret)
library(rgdal)
library(ranger)
library(raster)
library(scales)
library(lattice)
source("/data/models/saveRDS_functions.R")

rf.lst = paste0("mRF",c("_soiltype","_textype",".soc",".ph_h2o",".sand",".clay",".bd",".n_tot"),".rds")
rf.names = c("Soil great group","Particle size (PSCS)", "Soil organic carbon", "pH in H2O", "Sand fraction", "Clay fraction", "Bulk density", "Total (organic) N")

pdf(file = "Fig_RF_importance_plots_ALL_VARS.pdf", width = 7, height = 10)
par(mfrow=c(3,3), mar=c(2.5,7.5,2.5,0.5), oma=c(1,1,1,1))
for(i in 1:length(rf.lst)){
  mrfX = readRDS.gz(paste0(rf.lst[i]))
  xl = as.list(ranger::importance(mrfX))
  xl = t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]]))
  plot(x=rev(xl)/max(xl)*100, y=1:25, pch = 19, col="blue", xlab="Importance (%)", xlim=c(0,105), ylim=c(0,26), yaxp=c(0,25,25), xaxs="i", yaxs="i", cex=0.8, yaxt="n", ylab="", main=rf.names[i], cex.main=1.1)
  abline(h=1:25, lty=2, col="grey")
  axis(2, at=1:25, labels=substr(rev(attr(xl, "dimnames")[[1]]), 1, 12), las=2)
}
dev.off()

r <- raster("/data/USA48/elev48i0100a.tif")
xy.gg = readRDS.gz("ovA_TAX_gg.pnts.rds")
xy.gg = xy.gg[!is.na(xy.gg$X),c("SOURCEDB","X","Y","soiltype")]
coordinates(xy.gg) = ~ X + Y
proj4string(xy.gg) = proj4string(r)
xy.peds = readRDS.gz("ovA_NCSS_peds.rds")
xy.peds = xy.peds[!duplicated(xy.peds$LOC_ID)&!is.na(xy.peds$soc),c("soc","clay","LONWGS84","LATWGS84")]
xy.peds = xy.peds[!is.na(xy.peds$LONWGS84),]
coordinates(xy.peds) = ~ LONWGS84+LATWGS84
proj4string(xy.peds) = "+proj=longlat +datum=WGS84"
xy.peds = spTransform(xy.peds, CRS(proj4string(r)))

require(maptools)
require(maps)
country.m <- map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
country = as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
proj4string(country) = "+proj=longlat +datum=WGS84"
country <- spTransform(country, CRS(proj4string(r)))

png(file = "Fig_USA48_distribution_training_points.png", res = 150, width = 1650, height = 2200)
par(mfrow=c(2,1), mar=c(1,0,1,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", lwd=2, cex.main=1.2, main=paste0("NASIS points (N=", formatC(length(xy.gg), format="d", big.mark=","),")"))
points(xy.gg, pch=21, bg=alpha("red", 0.2), cex=.5, col=alpha("black", 0.3))
plot(country, col="darkgrey", lwd=2, cex.main=1.2, main=paste0("NCSS DB (N=", formatC(length(xy.peds), format="d", big.mark=","),")"))
points(xy.peds, pch=21, bg=alpha("blue", 0.2), cex=.5, col=alpha("black", 0.3))
dev.off()

