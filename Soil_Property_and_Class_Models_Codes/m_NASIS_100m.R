## Fit models for USDA GreatGroups, texture classes, and several soil properties using the NASIS points (ca 350,000 with classification and 45,000 points with soil properties)
## Code by Tom.Hengl@isric.org / points prepared by Amanda Rachmaran (a.m.ramcharan@gmail.com) and Travis Nauman (tnauman@usgs.gov)

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "ROCR", "randomForest", "R.utils", "plyr", "parallel", "psych", "mda", "dismo", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "mxnet", "doParallel", "caret", "plotKML")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/NASIS")
load(".RData")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
library(caret)
library(ranger)
library(xgboost)
library(nnet)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(raster)
library(plotKML)
library(GSIF)
#library(randomForestSRC)
#options(rf.cores=detectCores(), mc.cores=detectCores())
library(parallel)
library(doParallel)
#library(mxnet)
source("/data/models/saveRDS_functions.R")
source("/data/models/wrapper.predict_cs.R")
source("/data/NASIS/functions_NASIS.R")

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "gdal_translate"
  gdalwarp =  "gdalwarp"
  gdalbuildvrt = "gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
}

## USDA soil classes ----
levsGG <- read.csv("TAXOUSDA_GreatGroups.csv")
## Grid definition
r <- raster("/data/USA48/elev48i0100a.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
NODATA_value = -32768
te = paste(extent(r)[c(1,3,2,4)], collapse = " ")

## check if all layers are ready:
des <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_Covs100m.csv")
#s = raster::stack(paste0('/data/GEOG/NASIS/covs100m/', des$WORLDGRIDS_CODE, ".tif"))
#str(s[1])

## legends ----
#LNDCOV6.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_LandCover.csv")
NLCD116.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_LandCover_2011.csv")
PMTGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_pmaterial.csv")
DRNGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_drainage.csv")
PVEGKT6.leg <- read.csv("/data/GEOG/NASIS/covs100m/potential_vegetation_legend.csv")
COUNTY6.leg <- read.csv("/data/GEOG/NASIS/covs100m/Counties_legend.csv")
#GEOUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/geology_legend.csv")
GESUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/surfacegeology_legend.csv")

## Texture families ----
#system('7za x NASIS_L48_pscs_edt.zip')
#NASISpscs.pnts = readOGR("pscsmod_L48_edits.shp", "pscsmod_L48_edits")
system('7za x nasispts_pscsmodorgL48.zip')
NASISpscs.pnts = readOGR("nasispts_pscsmodorgL48.shp", "nasispts_pscsmodorgL48")
#with 308,291 features
#It has 5 fields
summary(NASISpscs.pnts$xwgs84)
NASISpscs.pnts$LOC_ID = paste("ID", NASISpscs.pnts$xwgs84, NASISpscs.pnts$ywgs84, sep="_")

## Training points (taxonomy, texture classes):
TAX_gg.pnts = readRDS("TAX_gg.pnts.rds")
#str(TAX_gg.pnts@data)
# 339,977 obs. of  7 variables:

## Soil properties NRCS NCSS ----
#system("7za e Geo_Peds07272016.csv.gz")
#NCSS_peds = read.csv("Geo_Peds07272016.csv")
NCSS_peds = readRDS("Geo_Peds12.20.2016.rds")
## 213,531 horizons
summary(NCSS_peds$long)
NCSS_peds$LOC_ID <- paste("ID", NCSS_peds$long, NCSS_peds$lat, sep="_")
t.props = c("clay", "sand", "soc", "bd", "n_tot", "ph_h2o", "ec_12pre", "zn_ppm","mg_pct", "k_pct", "p_ppm")  
## clay, sand, soil organic carbon, total nitrogen, bulk density and electrical conductivity
## https://github.com/ISRICWorldSoil/US48SoilGrids/issues/26
NCSS_peds = NCSS_peds[,c("pedon_key","LOC_ID","long","lat","hzn_top","hzn_bot",t.props[1:7])]
str(NCSS_peds)
summary(NCSS_peds$ec_12pre)

## RaCA points ----
#RaCA.pnts = readOGR("Export_RaCAsites.shp", "Export_RaCAsites")
## https://github.com/aramcharan/SoilGrids100m/blob/master/RaCA_data_prep.R
RaCA.pnts = read.csv("RaCA_12.6.2016.csv")
summary(RaCA.pnts$soc)
summary(RaCA.pnts$Measure_BD)
summary(RaCA.pnts$n_tot_ncs)
RaCA.pnts = RaCA.pnts[,c("upedonid","Lat","Long","soc","TOP","BOT","Measure_BD","n_tot_ncs")]
RaCA.pnts = plyr::rename(RaCA.pnts, replace=c("upedonid"="pedonid", "Lat"="lat","Long"="long","TOP"="hzn_top", "BOT"="hzn_bot", "Measure_BD"="bd","n_tot_ncs"="n_tot"))
summary(RaCA.pnts$long)
## 113920 NA
RaCA.pnts = RaCA.pnts[!is.na(RaCA.pnts$long),]
RaCA.pnts$LOC_ID <- paste("ID", RaCA.pnts$long, RaCA.pnts$lat, sep="_")

## Geochemical points ----
## https://mrdata.usgs.gov/ds-801/
## 28,749
geochem = readRDS("/data/Geochem/geochem_mrdata.rds")
length(!duplicated(geochem$LOC_ID))
geochem = plyr::rename(geochem, c("latitude"="lat", "longitude"="long", "depth_upper"="hzn_top", "depth_lower"="hzn_bot"))
geochem$pedonid = paste0("MR", geochem$site_id, sep="_")

NCSS_pedsM = dplyr::bind_rows(list(NCSS_peds, RaCA.pnts, geochem[,c("LOC_ID","pedonid","lat","long","hzn_top","hzn_bot","zn_ppm","mg_pct","k_pct","p_ppm")]))
str(NCSS_pedsM)
## 273495 obs. of  18 variables
tsum = lapply(t.props, function(i){summary(NCSS_pedsM[,i])})
names(tsum) = t.props
tsum
write.csv(do.call(rbind, tsum), "tvars_quantiles.csv")

## Unique locations ----
pnts = data.frame(LOC_ID=unique(c(TAX_gg.pnts$LOC_ID, NASISpscs.pnts$LOC_ID, NCSS_pedsM$LOC_ID)), LONWGS84=NA, LATWGS84=NA)
pnts$LONWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][2]}))
pnts$LATWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][3]}))
sel.pnts = pnts$LONWGS84< -59 & pnts$LONWGS84> -127 & pnts$LATWGS84< 52 & pnts$LATWGS84> 24 & !is.na(pnts$LONWGS84)
pnts = pnts[sel.pnts,]
str(pnts)
## 319397 obs. of  3 variables
coordinates(pnts) <- ~ LONWGS84 + LATWGS84
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")

## Tiling system ----
obj <- GDALinfo("/data/GEOG/NASIS/covs100m/COUNTY6.tif")
tile.lst <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=TRUE)
proj4string(tile.lst)
tile.tbl <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## 6300 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
writeOGR(tile.pol, "tiles50km.shp", "tiles50km", "ESRI Shapefile")
saveRDS(tile.pol, "tiles50km.rds")
plot(spTransform(pnts, CRS(proj4string(tile.lst))), pch="+")
lines(tile.lst, col="red")
## Overlay tiles and admin units (fully parallelized):
system(paste('gdal_translate /data/GEOG/NASIS/covs100m/COUNTY6.tif /data/USA48/COUNTY6.sdat -ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
system(paste0(saga_cmd, ' -c=48 shapes_grid 2 -GRIDS=\"/data/USA48/COUNTY6.sgrd\" -POLYGONS=\"tiles50km.shp\" -PARALLELIZED=1 -RESULT=\"ov_COUNTY6_tiles50km.shp\"'))
ov_COUNTY6 = readOGR("/data/NASIS/ov_COUNTY6_tiles50km.shp", "ov_COUNTY6_tiles50km")
summary(sel.t <- !ov_COUNTY6@data[,"COUNTY6..MA"]==-32767)
## 3392 tiles with values
ov_COUNTY6 = ov_COUNTY6[sel.t,]
str(ov_COUNTY6@data)

## Prepare covariates as tiles ----
t.sel = as.character(ov_COUNTY6$ID)
new.dirs <- paste0("/data/tt/NASIS/covs100m/T", t.sel)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
x <- lapply(gsub("covs100m", "predicted100m", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## run in parallel:
covs.lst = as.character(des$WORLDGRIDS_CODE[-which(des$WORLDGRIDS_CODE=="COUNTY6")])

## Clean up:
#rds.lst = list.files(path="/data/tt/NASIS/covs100m", pattern=glob2rx("*.rds"), recursive = TRUE, full.names = TRUE)
#unlink(rds.lst)
#tif.lst = list.files(path="/data/tt/NASIS/covs100m", pattern=glob2rx("*_predicted_*.tif"), recursive = TRUE, full.names = TRUE)
#unlink(tif.lst)
## Test:
#make_RDS_tiles(i=5208, tile.tbl=tile.tbl, covs.lst=covs.lst, NLCD116.leg=NLCD116.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)
#make_RDS_tiles(i=2018, tile.tbl=tile.tbl, covs.lst=covs.lst, NLCD116.leg=NLCD116.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)
#make_RDS_tiles(i=5609, tile.tbl=tile.tbl, covs.lst=covs.lst, NLCD116.leg=NLCD116.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)
#make_RDS_tiles(i=5525, tile.tbl=tile.tbl, covs.lst=covs.lst, NLCD116.leg=NLCD116.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)

## TAKES ca 8 HOURS TO PREPARE
sfInit(parallel=TRUE, cpus=46) ## 25
sfExport("make_RDS_tiles", "tile.tbl", "covs.lst", "NLCD116.leg", "PMTGSS7.leg", "DRNGSS7.leg", "PVEGKT6.leg", "COUNTY6.leg", "GESUSG6.leg")
sfLibrary(plyr)
sfLibrary(rgdal)
sfLibrary(ranger)
sfLibrary(sp)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( make_RDS_tiles(i, tile.tbl=tile.tbl, covs.lst=covs.lst, NLCD116.leg=NLCD116.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg) ) } )
sfStop()
save.image()

## Clean up empty dirs
pr.dirs <- basename(dirname(list.files(path="/data/tt/NASIS/covs100m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
## 3392
pr.dirs.c <- list.dirs("/data/tt/NASIS/covs100m")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
#x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})

## Gap-filled Parent material and drainage maps:
unlink(paste0('/data/GEOG/NASIS/covs100m/', c("PMTGSS7","DRNGSS7"), '_f.tif'))
for(j in c("PMTGSS7","DRNGSS7")){
  if(!file.exists(paste0('/data/GEOG/NASIS/covs100m/', j, '_f.tif'))){
    pm.tif = list.files(path="/data/tt/NASIS/covs100m", pattern=paste0(j, "_predicted"), recursive = TRUE, full.names = TRUE)
    cat(pm.tif, sep="\n", file="my_liste.txt")
    system(paste0(gdalbuildvrt, ' -input_file_list my_liste.txt ', j, '_f.vrt'))
    system(paste0(gdalwarp, ' ', j, '_f.vrt /data/GEOG/NASIS/covs100m/', j, '_f.tif -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000 -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
  }
}
save.image()

## SP OVERLAY ----  
## (takes ca 15 mins)
ov <- extract.tiled(x=pnts, tile.pol=tile.pol, path="/data/tt/NASIS/covs100m", ID="ID", cpus=48)
## filter out missing values
#ov$PMTGSS7 = as.factor(ifelse(ov$PMTGSS7=="", "NULL", ifelse(ov$PMTGSS7=="NA", NA, paste(ov$PMTGSS7))))
ov$PMTGSS7 = as.factor(ov$PMTGSS7)
#summary(ov$PMTGSS7)
levels(ov$PMTGSS7)
ov$DRNGSS7 = as.factor(ov$DRNGSS7)
summary(ov$DRNGSS7)
#str(ov)
## "PMTGSS7" has probably too many levels (for a RF model)
## http://stackoverflow.com/questions/14921805/convert-a-factor-to-indicator-variables
PMTGSS7.mat = data.frame(model.matrix(~PMTGSS7-1, ov))
xsum.PMTGSS7 = sapply(PMTGSS7.mat, sum, na.rm=TRUE)
## 10 categories have <10 observations
sel.rm.PMTGSS7 = names(PMTGSS7.mat)[which(xsum.PMTGSS7<6)]
## Land cover:
summary(ov$NLCD116)
NLCD116.mat = data.frame(model.matrix(~NLCD116-1, ov))
xsum.NLCD116 = sapply(NLCD116.mat, sum, na.rm=TRUE)
sel.rm.NLCD116 = names(NLCD116.mat)[which(xsum.NLCD116<6)]
## Potential vegetation:
PVEGKT6.mat = data.frame(model.matrix(~PVEGKT6-1, ov))
xsum.PVEGKT6 = sapply(PVEGKT6.mat, sum, na.rm=TRUE)
sel.rm.PVEGKT6 = names(PVEGKT6.mat)[which(xsum.PVEGKT6<6)]
## Drainage class:
summary(ov$DRNGSS7)
DRNGSS7.mat = data.frame(model.matrix(~DRNGSS7-1, ov))
xsum.DRNGSS7 = sapply(DRNGSS7.mat, sum, na.rm=TRUE)
sel.rm.DRNGSS7 = names(DRNGSS7.mat)[which(xsum.DRNGSS7<6)]

## Merge all indicators with other covariates ----
ov.fs = plyr::join_all(list(ov, 
                            cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PMTGSS7.mat))]), PMTGSS7.mat[,names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)]]), 
                            cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(NLCD116.mat))]), NLCD116.mat[,names(NLCD116.mat)[which(xsum.NLCD116>5)]]), 
                            cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PVEGKT6.mat))]), PVEGKT6.mat[,names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)]]), 
                            cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(DRNGSS7.mat))]), DRNGSS7.mat[,names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)]])))
## remove layers with too many missing values?
c.check = sapply(ov.fs, function(i){sum(is.na(i))})
names(ov.fs)[c.check>4000]
#ov.fs$WETGSS7 = ifelse(is.na(ov.fs$WETGSS7), 0, ov.fs$WETGSS7)
#summary(ov.fs$WETGSS7)

## join sampled values and covariates ----
ovA_NCSS_peds = plyr::join(NCSS_pedsM, ov.fs, by="LOC_ID", type="left")
## Fix some non-sensical values:
ovA_NCSS_peds$soc = ifelse(ovA_NCSS_peds$soc<0, 0, ovA_NCSS_peds$soc)
ovA_NCSS_peds$ph_h2o = ifelse(ovA_NCSS_peds$ph_h2o<2, NA, ovA_NCSS_peds$ph_h2o)
ovA_TAX_gg.pnts = plyr::join(as.data.frame(TAX_gg.pnts), ov.fs, by="LOC_ID", type="left")
ovA_NASISpscs.pnts = plyr::join(as.data.frame(NASISpscs.pnts), ov.fs, by="LOC_ID", type="left")
#summary(ovA_NASISpscs.pnts$pscsmodorg)
xs = summary(ovA_NASISpscs.pnts$pscsmodorg, maxsum=length(levels(ovA_NASISpscs.pnts$pscsmodorg)))
sel.levs = attr(xs, "names")[xs > 5]
ovA_NASISpscs.pnts$textype <- ovA_NASISpscs.pnts$pscsmodorg
ovA_NASISpscs.pnts$textype[which(!ovA_NASISpscs.pnts$pscsmodorg %in% sel.levs)] <- NA
ovA_NASISpscs.pnts$textype <- droplevels(ovA_NASISpscs.pnts$textype)
## Regression matrix per variable:
saveRDS.gz(ovA_NCSS_peds, file="ovA_NCSS_peds.rds")
saveRDS.gz(ovA_TAX_gg.pnts, file="ovA_TAX_gg.pnts.rds")
saveRDS.gz(ovA_NASISpscs.pnts, file="ovA_NASISpscs.pnts.rds")
#write.csv(ovA_TAX_gg.pnts[1:200,], "ovA_TAX_gg_pnts_example.csv")
save.image()

## covariates - USE ONLY INDICATORS (NO FACTOR-TYPE predictors)
pr.lst <- des$WORLDGRIDS_CODE
## tidy up list of covariates for soil-class mapping
prC.lst = c(paste(pr.lst[-c(grep(pattern = "GESUSG6",pr.lst), grep(pattern = "PMTGSS7",pr.lst), grep(pattern = "PVEGKT6",pr.lst), grep(pattern = "NLCD116",pr.lst), grep(pattern = "DRNGSS7",pr.lst), grep(pattern = "COUNTY6",pr.lst), grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst), grep(pattern = glob2rx("N??PRI5"),pr.lst), grep(pattern = glob2rx("V??PRI5"),pr.lst))]), names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)], names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)], names(NLCD116.mat)[which(xsum.NLCD116>5)], names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)])
str(prC.lst)
## 225 covs
prT.lst = c(prC.lst, paste(pr.lst[c(grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst))]))
str(prT.lst)
## 253
## '"DRNGSS7NULL"'? 
## https://docs.google.com/spreadsheets/d/1y0nMRrdPahnI5Et3AsV8GZV469ug34II5ZPJs4gbhKk/edit#gid=384507662

## Clean up:
#rds.lst = list.files(path="./covs100m", pattern=glob2rx("cT*.rds"), recursive = TRUE, full.names = TRUE)
#unlink(rds.lst)

mean.val = lapply(ov.fs[,prT.lst], function(x){quantile(x, probs=.5, na.rm=TRUE)})
names(mean.val) = prT.lst
#make_newdata(i=2019, in.path="/data/tt/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prT.lst, mean.val=mean.val)
#make_newdata(i=5208, in.path="/data/tt/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prT.lst, mean.val=mean.val)
## Make prediction locations in parallel (TAKES 1.5 HRS):
sfInit(parallel=TRUE, cpus=48)
sfExport("t.sel", "make_newdata", "PMTGSS7.leg", "DRNGSS7.leg", "prT.lst", "mean.val")
sfLibrary(plyr)
sfLibrary(sp)
x <- sfClusterApplyLB(as.numeric(t.sel), fun=function(i){ try( make_newdata(i, in.path="/data/tt/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prT.lst, mean.val=mean.val) ) } ) 
sfStop()
save.image()

## FIT MODELS FOR FACTORS ----

t.fvars = c("soiltype", "textype")
t.names = list("USDA Great Group", "SPCS classes")
names(t.names) = t.fvars

unlink(paste0("t.mRF_", t.fvars,".rds"))
unlink(paste0("mRF_", t.fvars,".rds"))

## fit models for factor variables
## TAKES ca 40 mins per variable and >190GB RAM
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
Nsub = 2e4
for(j in t.fvars){
  if(!file.exists(paste0("mRF_", j,".rds"))){
    formulaString = as.formula(paste(j, ' ~ ', paste(prC.lst, collapse="+")))
    if(j=="soiltype"){ x = ovA_TAX_gg.pnts[,all.vars(formulaString)]  }
    if(j=="textype"){ x = ovA_NASISpscs.pnts[,all.vars(formulaString)]  }
    x = x[complete.cases(x),] ## 325,215 obs for 'soiltype'
    x[,j] = droplevels(x[,j])
    rf.tuneGrid <- expand.grid(mtry = seq(5,round(sqrt(ncol(x))*2),by=5))
    if(!file.exists(paste0("t.mRF_", j,".rds"))){
      ## estimate Mtry using cross-validation (subsample original RM)
      xs = x[sample.int(nrow(x), Nsub),]
      xs[,j] = droplevels(xs[,j])
      cl <- makeCluster(48)
      registerDoParallel(cl)
      t.mrfX <- caret::train(formulaString, data=xs, method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid, num.trees=85)
      stopCluster(cl)    
      saveRDS.gz(t.mrfX, paste0("t.mRF_", j,".rds"))
    } else {
      t.mrfX = readRDS.gz(paste0("t.mRF_", j,".rds"))
    }
    ## How Many Trees in a Random Forest? 65-120 is more than enough (http://dx.doi.org/10.1007/978-3-642-31537-4_13)
    gc(); gc()
    mRF = ranger::ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, probability=TRUE, num.trees=85, mtry=t.mrfX$bestTune$mtry) ## optimize mtry based on caret
    saveRDS.gz(mRF, file=paste0("mRF_", j,".rds"))
    cat("Results of model fitting 'randomForest':\n", file=paste0("NASIS_resultsFit_", j, ".txt"))
    cat("\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    cat(paste("Variable:", t.names[[j]]), file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    cat("\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    sink(file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE, type="output")
    cat("\n Random forest model:", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    print(mRF)
    cat("\n Variable importance:\n", file=paste0("NASIS_resultsFit_", j, ".txt"), append=TRUE)
    xl <- as.list(ranger::importance(mRF))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:150]])))
    sink()
  }
}
closeAllConnections()
gc(); gc()
#object.size(mRF)/1e6

## PREDICTIONS - CLASSES ----

# del.lst <- list.files(path="/data/tt/NASIS/predicted100m/", pattern=glob2rx("^TAXgg*.rds$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst) 
# del.lst <- list.files(path="/data/tt/NASIS/predicted100m/", pattern=glob2rx("^PSCS*.tif$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)
# del.lst <- list.files(path="/data/tt/NASIS/predicted100m/", pattern=glob2rx("^TAXgg*.tif$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst) 
## Delete all pridictions:
#del.lst <- list.files(path="/data/tt/NASIS/predicted100m/", pattern=glob2rx("*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst) 

## Predict classes in parallel:
cl.lst = c("soiltype", "textype")
varn.lst = c("TAXgg", "PSCS")
#out.path = "/data/tt/NASIS/predicted100m/"
## Test predictions:
#mRF = readRDS.gz(paste0("mRF_", cl.lst[1], ".rds")
#varn = varn.lst[1]
#make_newdata(i=1923, in.path="/data/tt/NASIS/covs100m", PMTGSS7.leg, DRNGSS7.leg, independent.variable.names=prC.lst)
#predict_factor_tile(i=1922, mRF, varn, levs, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/")

for(k in 1:length(cl.lst)){
  mRF = readRDS.gz(paste0("mRF_", cl.lst[k], ".rds"))
  varn = varn.lst[k]
  ## 291 class for soiltype
  ## 65 for textype
  levs = attr(mRF$predictions, "dimnames")[[2]]
  ## run predictions in parallel: ca 30 HRS = 3600 tiles * 1 min / 2 tiles / 60 mins
  #n.cores = round(200/(object.size(mRF)/1e9))
  detach("package:snowfall", unload=TRUE)
  detach("package:snow", unload=TRUE)
  if(varn=="TAXgg"){ 
    cl <- parallel::makeCluster(12, type="FORK")
    x = parallel::parLapply(cl, as.numeric(t.sel), fun=function(i){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/T", i, "/", varn, "_T", i, ".tif"))){ try( predict_factor_tile(i, mRF, varn, levs, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m", mc.cores=5) ) } } )
    parallel::stopCluster(cl)
  }
  if(varn=="PSCS"){ 
    cl <- parallel::makeCluster(24, type="FORK")
    x = parallel::parLapply(cl, as.numeric(t.sel), fun=function(i){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/T", i, "/", varn, "_T", i, ".tif"))){ try( predict_factor_tile(i, mRF, varn, levs, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/") ) } } )
    parallel::stopCluster(cl)
  } 
  gc(); gc()
  ## mosaic tiles:
  library(snowfall)
  unlink(list.files(path="/data/GEOG/NASIS/predicted100m", pattern=paste0(varn, "_"), full.names=TRUE))
  levs = list.files(path="/data/tt/NASIS/predicted100m/T5445", pattern=glob2rx(paste0("^",varn,"_*_T*.tif$")))
  levs = sapply(basename(levs), function(x){strsplit(x, "_")[[1]][2]})
  sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
  sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_100m", "varn", "te")
  out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_100m(x, in.path="/data/tt/NASIS/predicted100m/", varn=varn, te=te) )})
  sfStop()
  if(varn=="TAXgg"){ write.csv(data.frame(Value=1:length(mRF$forest$levels), Class=mRF$forest$levels), file="/data/GEOG/NASIS/predicted100m/TAXgg_M_100m.tif.csv") }
  if(varn=="PSCS"){ write.csv(data.frame(Value=1:length(mRF$forest$levels), Class=mRF$forest$levels), file="/data/GEOG/NASIS/predicted100m/PSCS_M_100m.tif.csv") }
}

## Dominant class:
mosaic_tiles_100m(j="M", in.path="/data/tt/NASIS/predicted100m/", varn="TAXgg", te=te, r="near", ot="Int16", dstnodata=-32768)
mosaic_tiles_100m(j="M", in.path="/data/tt/NASIS/predicted100m/", varn="PSCS", te=te, r="near")

## Clean up emty dirs ----
pr.dirs <- basename(dirname(list.files(path="/data/tt/NASIS/predicted100m/", pattern=glob2rx("PSCS_T*.tif$"), recursive = TRUE, full.names = TRUE)))
## 3392 - 3379 = 13 empty dirs
pr.dirs.c = list.dirs("/data/tt/NASIS/predicted100m/")[-1]
pr.dirsT <- list.dirs("/data/tt/NASIS/covs100m")[-1]
selD <- which(!basename(pr.dirs.c) %in% basename(pr.dirsT))
#x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})
save.image()
## Clean up new.data objects that did not produce any predictions:
pr.crds <- list.files(path="/data/tt/NASIS/covs100m", pattern=glob2rx("^cT*.rds$"), recursive = TRUE, full.names = TRUE)
selD <- pr.crds[which(!basename(dirname(pr.crds)) %in% basename(pr.dirs))]
unlink(selD)

## FIT MODELS FOR NUMERIC VARIABLES ----

#str(ovA_NCSS_peds)
ovA_NCSS_peds$DEPTH = ovA_NCSS_peds$hzn_top + (ovA_NCSS_peds$hzn_bot - ovA_NCSS_peds$hzn_top)/2
summary(ovA_NCSS_peds$DEPTH)
## 7 standard depths
breaks = c(0, rowMeans(data.frame(c(0,5,15,30,60,100), c(5,15,30,60,100,200))), 200, 4500)
ovA_NCSS_peds$DEPTH_c = cut(ovA_NCSS_peds$DEPTH, breaks, labels = paste0("sl", 1:8))
summary(ovA_NCSS_peds$DEPTH_c)
## 254 NA's
hist(ovA_NCSS_peds$DEPTH) ## some profiles are up to 35m deep?
## Estimate SG variables per depth:
sg.props = c("CLYPPT", "SNDPPT", "ORCDRC", "BLDFIE", rep("", length(t.props)-4))
#View(data.frame(t.props,sg.props))
for(j in 1:length(t.props)){
  if(!(t.props[j]=="n_tot"|t.props[j]=="ec_12pre"|t.props[j]=="zn_ppm"|t.props[j]=="mg_pct"|t.props[j]=="k_pct"|t.props[j]=="p_ppm"|t.props[j]=="ph_h2o")){
    ovA_NCSS_peds[,paste0("sg_", t.props[j])] = ifelse(ovA_NCSS_peds$DEPTH_c=="sl1", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl1_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl2", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl2_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl3", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl3_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl4", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl4_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl5", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl5_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl6", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl6_100m")], ifelse(ovA_NCSS_peds$DEPTH_c=="sl7", ovA_NCSS_peds[,paste0(sg.props[j], "_M_sl7_100m")], NA)))))))
  }
}
#xyplot(clay~sg_clay, ovA_NCSS_peds)
NCSS_peds.pnt = ovA_NCSS_peds[,c("LOC_ID","long","lat","hzn_top","hzn_bot","DEPTH",t.props)]
NCSS_peds.pnt = NCSS_peds.pnt[!is.na(NCSS_peds.pnt$long),]
coordinates(NCSS_peds.pnt) <- ~ long + lat
proj4string(NCSS_peds.pnt) <- CRS("+proj=longlat +datum=WGS84")
#unlink("NCSS_peds.shp")
#writeOGR(NCSS_peds.pnt, "NCSS_peds.shp", "NCSS_peds", "ESRI Shapefile")
unlink("NCSS_peds.gpkg")
writeOGR(NCSS_peds.pnt[,c("LOC_ID","hzn_top","hzn_bot","DEPTH",t.props)], "NCSS_peds.gpkg", "NCSS_peds", "GPKG")

## fit models in a loop
## Generic settings for caret:
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:4, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample=1)
rf.tuneGrid <- expand.grid(mtry = seq(5,50,by=5))
Nsub = 1e4

#unlink(paste0("t.mRF.", t.props,".rds"))
#unlink(paste0("mRF.", t.props,".rds"))
#unlink(paste0("mGB.", t.props,".rds"))
save.image()

library(parallel)
cl <- parallel::makeCluster(48)
doParallel::registerDoParallel(cl)
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file="SPROPS_resultsFit.txt")
for(j in t.props){
  out.rf <- paste0("mRF.", j,".rds")
  if(!file.exists(out.rf)){
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    cat(paste("Variable:", j), file="SPROPS_resultsFit.txt", append=TRUE)
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    if(j=="n_tot"|j =="ec_12pre"|j =="ph_h2o"|j =="zn_ppm"|j =="mg_pct"|j =="k_pct"|j =="p_ppm"){
      formulaString = as.formula(paste(j, ' ~ DEPTH + ', paste(prC.lst, collapse="+")))
    } else {
      formulaString = as.formula(paste0(j, ' ~ DEPTH + ', paste(prC.lst, collapse="+"), "+ sg_", j))
    }
    x = ovA_NCSS_peds[,all.vars(formulaString)]
    x = x[complete.cases(x),]
    ## Caret training settings (reduce number of combinations to speed up):
    ## optimize mtry parameter:
    if(!file.exists(gsub("mRF","t.mRF",out.rf))){
      t.mrfX <- caret::train(formulaString, data=x[sample.int(nrow(x), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
      saveRDS.gz(t.mrfX, file=gsub("mRF","t.mRF",out.rf))
    } else {
      t.mrfX <- readRDS.gz(gsub("mRF","t.mRF",out.rf))
    }
    ## fit RF model using 'ranger' (fully parallelized)
    ## reduce number of trees so the output objects do not get TOO LARGE i.e. >5GB
    mrfX <- ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
    saveRDS.gz(mrfX, file=out.rf)
    ## Top covariates:
    sink(file="SPROPS_resultsFit.txt", append=TRUE, type="output")
    print(mrfX)
    cat("\n Variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
    xl <- as.list(ranger::importance(mrfX))
    print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:150]])))
    ## XGBoost
    if(!file.exists(paste0("mGB.",j,".rds"))){
      mGBX <- caret::train(formulaString, data=x, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid)
      saveRDS.gz(mGBX, file=paste0("mGB.",j,".rds"))
      ## save also binary model for prediction purposes:
      xgb.save(mGBX$finalModel, paste0("Xgb.",j))
    } else {
      mGBX <- readRDS.gz(paste0("mGB.",j,".rds"))
    }
    importance_matrix <- xgb.importance(mGBX$coefnames, model = mGBX$finalModel)
    cat("\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(mGBX)
    cat("\n XGBoost variable importance:\n", file="SPROPS_resultsFit.txt", append=TRUE)
    print(importance_matrix[1:150,])
    cat("--------------------------------------\n", file="SPROPS_resultsFit.txt", append=TRUE)
    sink()
  }
}
rm(mrfX); rm(mGBX)
stopCluster(cl); closeAllConnections()

## PREDICTIONS - NUMERIC VARIABLES ----

z.min <- as.list(c(0,0,0,50,0,20,0,0,0,0,0))
names(z.min) = t.props
z.max <- as.list(c(100,100,710,3500,580,120,1600,15000,25000,1000,35000))
names(z.max) = t.props
type.lst <- c("Byte","Byte","Int16","Int16","Int16","Byte","Byte","Int16","Int16","Int16","Int16")
names(type.lst) = t.props
mvFlag.lst <- c(255, 255, -32768, -32768, -32768, 255, 255, -32768, -32768, -32768, -32768)
names(mvFlag.lst) = t.props
sg.var <- c("CLYPPT", "SNDPPT", "ORCDRC", "BLDFIE", rep("", length(t.props)-4))
names(sg.var) = t.props

## Run per property
## TAKES ABOUT 10hrs for RF models and about 2hrs for XGB models
for(j in t.props[c(1:7,9:10)]){
  if(!file.exists(paste0("/data/GEOG/NASIS/predicted100m/", j, "_M_sl1_100m.tif"))){
    ## it appears that 'parallel' package conflicts with snowfall package so needs to be turned off:
    detach("package:snowfall", unload=TRUE)
    detach("package:snow", unload=TRUE)
    if(j=="ph_h2o"|j=="soc"|j=="ec_12pre"|j=="p_ppm"|j=="zn_ppm"){ multiplier = 10 }
    if(j %in% c("n_tot","n_tot","k_pct")){ multiplier = 100 }
    if(j %in% c("bd")){ multiplier = 1000 }
    if(j %in% c("mg_pct")){ multiplier = 10000 }
    if(j %in% c("clay","sand")){ multiplier = 1 }
    ## Random forest predictions:
    gm = readRDS.gz(paste0("mRF.", j,".rds"))
    gm1.w = 1/gm$prediction.error
    cpus = unclass(round((256-80)/(3.5*(object.size(gm)/1e9))))
    cl <- makeCluster(ifelse(cpus>35, 35, cpus), type="FORK")
    if(j=="n_tot"|j=="zn_ppm"|j=="mg_pct"|j=="ph_h2o"|j=="ec_12pre"|j=="k_pct"|j=="p_ppm"){
      x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/", x, "/", j,"_", x, "_rf.rds"))){ try( split_predict_n(x, gm, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/", split_no=NULL, varn=j, method="ranger", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/tt/NASIS/covs100m/", x, "/c", x,".rds"), SG.col=NULL ) , silent = TRUE) } } )
    } else {
      x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/", x, "/", j,"_", x, "_rf.rds"))){ try( split_predict_n(x, gm, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/", split_no=NULL, varn=j, method="ranger", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/tt/NASIS/covs100m/", x, "/c", x,".rds"), SG.col=paste0(sg.var[j], "_M_sl", 1:7, "_100m") ) , silent = TRUE) } } )  
    }
    stopCluster(cl)
    gc(); gc()
    ## XGBoost:
    gm = readRDS.gz(paste0("mGB.", j,".rds"))
    gm2.w = 1/(min(gm$results$RMSE, na.rm=TRUE)^2)
    cpus = unclass(round((256-50)/(3.5*(object.size(gm)/1e9))))
    cl <- makeCluster(ifelse(cpus>35, 35, cpus), type="FORK")
    if(j=="n_tot"|j=="zn_ppm"|j=="mg_pct"|j=="ph_h2o"|j=="ec_12pre"|j=="k_pct"|j=="p_ppm"){
      x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/", x, "/", j,"_", x, "_xgb.rds"))){ try( split_predict_n(x, gm, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/", split_no=NULL, varn=j, method="xgboost", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/tt/NASIS/covs100m/", x, "/c", x,".rds"), SG.col=NULL ) , silent = TRUE) } } )
    } else {
      x = parLapply(cl, paste0("T", t.sel), fun=function(x){ if(!file.exists(paste0("/data/tt/NASIS/predicted100m/", x, "/", j,"_", x, "_xgb.rds"))){ try( split_predict_n(x, gm, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/", split_no=NULL, varn=j, method="xgboost", DEPTH.col="DEPTH", multiplier=multiplier, rds.file=paste0("/data/tt/NASIS/covs100m/", x, "/c", x,".rds"), SG.col=paste0(sg.var[j], "_M_sl", 1:7, "_100m") ) , silent = TRUE) } } )  
    }
    stopCluster(cl)
    rm(gm)
    gc(); gc()
    ## sum up predictions:
    if(is.nan(gm1.w)|is.nan(gm2.w)){ gm1.w = 0.55; gm2.w = 0.45 } 
    library(snowfall)
    sfInit(parallel=TRUE, cpus=45)
    sfExport("t.sel", "sum_predict_ensemble", "j", "z.min", "z.max", "gm1.w", "gm2.w", "type.lst", "mvFlag.lst")
    sfLibrary(rgdal)
    sfLibrary(plyr)
    x <- sfClusterApplyLB(paste0("T", t.sel), fun=function(x){ try( if(length(list.files(path = paste0("/data/tt/NASIS/predicted100m/", x, "/"), glob2rx(paste0("^", j, "_M_sl*_", x, ".tif$"))))==0){ try( sum_predict_ensemble(x, in.path="/data/tt/NASIS/covs100m", out.path="/data/tt/NASIS/predicted100m/", varn=j, num_splits=NULL, zmin=z.min[[j]], zmax=z.max[[j]], gm1.w=gm1.w, gm2.w=gm2.w, type=type.lst[[j]], mvFlag=mvFlag.lst[[j]], rds.file=paste0("/data/tt/NASIS/covs100m/", x, "/c", x,".rds")) , silent = TRUE) } )  } )
    sfStop()
    ## mosaic tiles:
    levs = paste0("M_sl", 1:7)
    sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
    sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_100m", "j", "te", "type.lst", "mvFlag.lst")
    out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_100m(x, in.path="/data/tt/NASIS/predicted100m/", varn=j, te=te, ot=type.lst[[j]], dstnodata=mvFlag.lst[[j]]) )})
    sfStop()
    closeAllConnections()
    gc(); gc()
  }
}

## corrupt or missing tiles:
sfInit(parallel=TRUE, cpus=length(t.props))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("missing.tiles", "t.props", "t.sel")
missing.lst <- sfLapply(paste0("T", t.sel), missing.tiles, pr.dirs=pr.dirs)
sfStop()
names(missing.lst) = t.props
str(missing.lst)

## Add metadata to GeoTiffs ----
metasd <- read.csv('/data/GEOG/NASIS/predicted100m/META_GEOTIFF_USA48.csv', stringsAsFactors = FALSE)
sel.metasd = names(metasd)[-sapply(c("FileName","VARIABLE_NAME"), function(x){grep(x, names(metasd))})]
## covariates:
filenames = list.files("/data/GEOG/NASIS/predicted100m/", pattern=glob2rx("*.tif$"))
#write.csv(data.frame(filenames), "/data/NASIS/tif_list_USA48.csv")
for(j in filenames){
  out.tif <- paste0("/data/GEOG/", j)
  metadata = metasd[which(metasd$FileName == j), sel.metasd]
  m = paste('-mo ', '\"', names(metadata), "=", as.vector(metadata), '\"', sep="", collapse = " ")
  command = paste0('gdal_edit.py ', m,' ', out.tif)
  system (command, intern=TRUE)
}
