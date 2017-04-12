## Prepare covariates for USA48 (100 m resolution)
## tom.hengl@isric.org

setwd("/data/USA48")
load(".RData")
library(R.utils)
library(snowfall)
library(raster)
library(RSAGA)
library(rgdal)
library(utils)
library(tools)
library(ranger)
#library(randomForestSRC)
library(parallel)
#options(rf.cores=detectCores(), mc.cores=detectCores())

source("/data/models/saveRDS_functions.R")

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

gz.lst <- list.files(pattern=glob2rx("*.gz$"))
sapply(gz.lst, gunzip)
sapply(gsub("tar.gz", "tar", gz.lst), untar)
## 100 m DEM from http://nationalmap.gov/small_scale/
system(paste(gdal_translate, 'elev48i0100a.tif DEM100m.sdat -ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
system(paste(gdal_translate, 'elev48i0100a.tif /data/GEOG/NASIS/covs100m/DEMNED6.tif -ot \"Int16\" -a_nodata \"-32768\" -co \"COMPRESS=DEFLATE\"'))
## 49810 x 31390 pixels!
system(paste(gdalwarp, 'elev48i0100a.tif DEM250m.sdat -ot \"Int16\" -of \"SAGA\" -tr 250 250'))
## Land cover classes:
system(paste0(gdalwarp, ' ldco48i0100a.tif /data/GEOG/NASIS/covs100m/LNDCOV6.tif -ot \"Byte\" -multi -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -r \"near\" -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
## Land cover class from 2011 (https://www.mrlc.gov/nlcd11_data.php)
system('7za x nlcd_2011_landcover_2011_edition_2014_10_10.zip')
system(paste0(gdalwarp, ' /data/USA48/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img /data/GEOG/NASIS/covs100m/NLCD116.tif -ot \"Byte\" -multi -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -r \"near\" -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
## Bare ground based on landsat
unlink("/data/GEOG/NASIS/covs100m/BARL10.tif")
system(paste0(gdalwarp, ' /data/Landsat/bare2010_100m.tif /data/GEOG/NASIS/covs100m/BARL10.tif -multi -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
system(paste0(gdalwarp, ' /data/Landsat/treecover_100m.tif /data/GEOG/NASIS/covs100m/TREL10.tif -multi -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))


## TWI costs too much time to compute at 100m hence we use the 250m:
#system(paste0(gdalwarp, ' /data/MDEM/MTWI_NA_250m.sdat DEM100m_twi.sdat -multi -of \"SAGA\" -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -tr 100 100 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))

## Derive some standard DEM variables of interest for soil mapping:
saga_DEM_derivatives <- function(INPUT, MASK=NULL, sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM","MRN","TPI"), cpus=46){
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    if(!file.exists(gsub(".sgrd", "_slope.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
    }
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    if(!file.exists(gsub(".sgrd", "_vbf.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
    }
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    if(!file.exists(gsub(".sgrd", "_vdepth.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
    }
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    if(!file.exists(gsub(".sgrd", "_openp.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
    }
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    if(!file.exists(gsub(".sgrd", "_devmean.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=11' ) ) )
    }
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    if(!file.exists(gsub(".sgrd", "_twi.sgrd", INPUT))){
      try( suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_hydrology 15 -DEM=\"', INPUT, '\" -SLOPE_MIN=0.04 -SLOPE_OFF=0.3 -AREA_MOD=\"', gsub(".sgrd", "_catchm.sgrd", INPUT), '\" -SLOPE_TYPE=0 -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) ) ## gsub("100", "250", INPUT)
    }
  }
  ## MRN:
  if(any(sel %in% "MRN")){
    if(!file.exists(gsub(".sgrd", "_mrn.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_hydrology 23 -DEM=\"', INPUT, '\" -AREA=\"tmp.sgrd\" -MRN=\"', gsub(".sgrd", "_mrn.sgrd", INPUT), '\" -ZMAX=\"tmp.sgrd\"' ) ) )
    }
  }
  ## TPI:
  if(any(sel %in% "TPI")){
    if(!file.exists(gsub(".sgrd", "_tpi.sgrd", INPUT))){
      suppressWarnings( system(paste0(saga_cmd, ' -c=', cpus,' ta_morphometry 18 -DEM=\"', INPUT, '\" -STANDARD=1 -TPI=\"', gsub(".sgrd", "_tpi.sgrd", INPUT), '\" -RADIUS_MIN=0 -RADIUS_MAX=2000 -DW_WEIGHTING=3 -DW_BANDWIDTH=75' ) ) )
    }
  }
}

saga_DEM_derivatives(INPUT="DEM100m.sgrd")

## reformat and compress images (use only integers for large geoTiffs):
out.lst = c("SLPNED6","CRVNED6","VBFNED6","DVMNED6","VDPNED6","NEGNED6","POSNED6","TWINED6","MRNNED6","TPINED6")
sf.lst = c(100, 10000, 100, 100, 10, 1000, 1000, 100, 100, 100)
v.lst = c("slope","cprof","vbf","devmean","vdepth","openp","openn","twi","mrn","tpi")

for(i in 1:length(v.lst)){
  if(!file.exists(paste0(out.lst[i],'.tif'))){
    tmp <- set.file.extension(tempfile(), ".sgrd")
    ## http://saga-gis.org/saga_module_doc/2.2.7/grid_calculus_1.html
    system(paste0(saga_cmd, ' -c=46 grid_calculus 1 -GRIDS DEM100m_', v.lst[i],'.sgrd -FORMULA g1*', sf.lst[i], ' -INTERPOLATION 0 -TYPE 7 -RESULT ', tmp))
    system(paste0(gdal_translate, ' ', set.file.extension(tmp, ".sdat"), ' /data/GEOG/NASIS/covs100m/', out.lst[i], '.tif -ot \"Int16\" -a_nodata \"-32768\" -co \"COMPRESS=DEFLATE\"'))
    unlink(set.file.extension(tmp, ".sdat"))
  }
}

## Parent material and drainage class from gSSURGO:
unlink("/data/GEOG/NASIS/covs100m/DRNGSS7.tif")
system('7za x fy2016_drncl.zip')
system(paste0(gdalwarp, ' fy2016_drncl.tif /data/GEOG/NASIS/covs100m/DRNGSS7.tif -t_srs \"', proj4string(r), '\"  -multi -r \"near\" -tr 100 100 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "), ' -dstnodata 255 -co \"COMPRESS=DEFLATE\"'))
fy2016_drncl = read.dbf("fy2016_drncl.dbf")
write.csv(fy2016_drncl[[1]], "fy2016_drainage.csv")

unlink("/data/GEOG/NASIS/covs100m/PMTGSS7.tif")
system('7za x fy2016_pm.zip')
system(paste0(gdalwarp, ' fy2016_pm.tif /data/GEOG/NASIS/covs100m/PMTGSS7.tif -t_srs \"', proj4string(r), '\"  -multi -r \"near\" -tr 100 100 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "), ' -dstnodata 255 -co \"COMPRESS=DEFLATE\"'))
fy2016_pm = read.dbf("fy2016_pm.dbf")
write.csv(fy2016_pm[[1]], "fy2016_pmaterial.csv")

## Potential vegetation:
potveg = read.dbf("US_potential_vegetation.dbf")
potveg$dbf$TYPE_INT = as.integer(potveg$dbf$FORM_LAB)
write.dbf(potveg, "US_potential_vegetation.dbf")
potveg.sum = data.frame(NAMES=levels(potveg$dbf$FORM_LAB), Value=1:length(levels(potveg$dbf$FORM_LAB)))
write.csv(potveg.sum, "potential_vegetation_legend.csv")
system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"US_potential_vegetation.shp\" -FIELD \"TYPE_INT\" -GRID \"potveg_100m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("PVEGKT6.tif")
system(paste(gdal_translate, 'potveg_100m.sdat /data/GEOG/NASIS/covs100m/PVEGKT6.tif -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\"'))

## Counties USA:
counties = read.dbf("cb_2015_us_county_500k.dbf")
counties$dbf$NAME_INT = as.integer(counties$dbf$NAME)
counties$dbf$GEOID_INT = as.integer(counties$dbf$GEOID)
write.dbf(counties, "cb_2015_us_county_500k.dbf")
counties.sum = data.frame(GEOID=levels(counties$dbf$GEOID), Value=1:length(levels(counties$dbf$GEOID)))
counties.sum$NAMES = join(counties.sum, counties$dbf)$NAMES
write.csv(counties.sum, "Counties_legend.csv")
system(paste('ogr2ogr -t_srs \"', proj4string(r), '\" cb_2015_us_county_500k_a.shp cb_2015_us_county_500k.shp'))
system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"cb_2015_us_county_500k_a.shp\" -FIELD \"GEOID_INT\" -GRID \"counties_100m.sgrd\" -GRID_TYPE 1 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("COUNTY6.tif")
system(paste(gdal_translate, 'counties_100m.sdat /data/GEOG/NASIS/covs100m/COUNTY6.tif -ot \"Int16\" -co \"COMPRESS=DEFLATE\"'))

## EarthEnv MODIS cloud images:
cld.lst = list.files(path = "/data/EarthEnv", pattern="MODCF", full.names = TRUE)
cld.lst
raster(cld.lst[1])
cldo.lst = c("MANMCF5.tif", paste0("C0", 1:9, "MCF5.tif"), paste0("C", 10:12, "MCF5.tif"), "CSCMCF5.tif")
sfInit(parallel=TRUE, cpus=length(cld.lst))
sfExport("gdalwarp", "cld.lst", "cldo.lst", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(cld.lst), function(i){ system(paste0(gdalwarp, ' ', cld.lst[i], ' ', cldo.lst[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))} )
sfStop()
## Warning 1: Input file /data/EarthEnv/MODCF_monthlymean_11.tif has a color table, which will likely lead to bad results when using a resampling method other than nearest neighbour or mode. Converting the dataset prior to 24/32 bit is advised.

system("7za e Bioclimatic_Layers.zip")
clm.lst = c("AnnDiurnRge.tif", "PptColdQ.tif", "PptDrst.tif", "PptDryQ.tif", "PptSeas.tif", "PptWarmQ.tif", "PptWetst.tif", "TColdQ.tif", "TDryQ.tif", "TWarmQ.tif", "TWetQ.tif")
clmo.lst = c("ADIUCL5.tif", "PCLQCL5.tif", "PDRSCL5.tif", "PDRYCL5.tif", "PSESCL5.tif", "PWRMCL5.tif", "PWETCL5.tif", "TCLDCL5.tif", "TDRYCL5.tif", "TWRMCL5.tif", "TWETCL5.tif")
raster(clm.lst[1])
sfInit(parallel=TRUE, cpus=length(clm.lst))
sfExport("gdalwarp", "clm.lst", "clmo.lst", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(clm.lst), function(i){ if(!file.exists(clmo.lst[i])){ system(paste0(gdalwarp, ' ', clm.lst[i], ' /data/GEOG/NASIS/covs100m/', clmo.lst[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))} } )
sfStop()

## Chelsa precipitation images:
chl.lst = c(paste0("/data/Climate/CHELSA_prec_",1:12,"_1979-2013.tif"), paste0("/data/Climate/CHELSA_prec_1979-2013_land.tif"))
chl.out = c(paste0("P0",1:9,"CHL5.tif"), paste0("P",10:12,"CHL5.tif"), "PRTCHL5.tif")
sfInit(parallel=TRUE, cpus=length(chl.lst))
sfExport("gdalwarp", "chl.lst", "chl.out", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(chl.lst), function(i){ if(!file.exists(chl.out[i])){ system(paste0(gdalwarp, ' ', chl.lst[i], ' /data/GEOG/NASIS/covs100m/', chl.out[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -ot \"Int16\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()


## MODIS EVI images at 250m:
m.lst = c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
evi.lst <- c(paste0("/data/MOD13Q1/M_liste", m.lst, ".vrt"), paste0("/data/MOD13Q1/SD_liste", m.lst, ".vrt"))
evi.out = c(paste0("EX", 1:6, "MOD5.tif"), paste0("ES", 1:6, "MOD5.tif"))
sfInit(parallel=TRUE, cpus=length(evi.lst))
sfExport("gdalwarp", "evi.lst", "evi.out", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(evi.lst), function(i){ if(!file.exists(evi.out[i])){ system(paste0(gdalwarp, ' ', evi.lst[i], ' ', evi.out[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -ot \"Int16\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()

## MODIS LST images:
m1.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
lst.lst = c(paste0("/data/MOD11A2/M_liste", m1.lst, "_D.vrt"), paste0("/data/MOD11A2/M_liste", m1.lst, "_N.vrt"))
lst.out = c(paste0("T0", 1:9, "MOD3.tif"), paste0("T", 10:12, "MOD3.tif"), paste0("N0", 1:9, "MOD3.tif"), paste0("N", 10:12, "MOD3.tif"))
sfInit(parallel=TRUE, cpus=length(lst.lst))
sfExport("gdalwarp", "lst.lst", "lst.out", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(lst.lst), function(i){ if(!file.exists(lst.out[i])){ system(paste0(gdalwarp, ' ', lst.lst[i], ' ', lst.out[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -ot \"Int16\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()

## PRISM data (http://www.prism.oregonstate.edu/normals/)
system("7za e PRISM_vpdmin_30yr_normal_800mM2_all_bil.zip")
system("7za e PRISM_tmin_30yr_normal_800mM2_all_bil.zip")
system("7za e PRISM_tmax_30yr_normal_800mM2_all_bil.zip")
system("7za e PRISM_ppt_30yr_normal_800mM2_all_bil.zip")
pri.lst = list.files(pattern=glob2rx("*.bil$"))
pri.out = paste0("/data/GEOG/NASIS/covs100m/", c(paste0("P0",1:9,"PRI5.tif"), paste0("P",10:12,"PRI5.tif"), "PRTPRI5.tif", paste0("T0",1:9,"PRI5.tif"), paste0("T",10:12,"PRI5.tif"), "TMTPRI5.tif", paste0("N0",1:9,"PRI5.tif"), paste0("N",10:12,"PRI5.tif"), "TNMPRI5.tif", paste0("V0",1:9,"PRI5.tif"), paste0("V",10:12,"PRI5.tif"), "VPDPRI5.tif"))
write.csv(data.frame(pri.out), "prism_codes.csv")
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalwarp", "pri.lst", "pri.out", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(pri.lst), function(i){ if(!file.exists(pri.out[i])){ system(paste0(gdalwarp, ' ', pri.lst[i], ' /data/GEOG/NASIS/covs100m/', pri.out[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()

## Geology of USA (https://mrdata.usgs.gov/geology/state/)
system("7za e geol_poly.zip")
geol = read.dbf("geol_poly.dbf")
#summary(geol$dbf$LITH62)
geol$dbf$LITH62_INT = as.integer(geol$dbf$LITH62)
write.dbf(geol, "geol_poly.dbf")
geol.sum = data.frame(NAMES=levels(geol$dbf$LITH62), Value=1:length(levels(geol$dbf$LITH62)))
write.csv(geol.sum, "geology_legend.csv")
system(paste('ogr2ogr -t_srs \"', proj4string(r), '\" geol_poly_a.shp geol_poly.shp'))
system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"geol_poly_a.shp\" -FIELD \"LITH62_INT\" -GRID \"geol_100m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("GEOUSG6.tif")
system(paste(gdal_translate, 'geol_100m.sdat /data/GEOG/NASIS/covs100m/GEOUSG6.tif -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\"'))

## Surface geology (http://pubs.usgs.gov/ds/425/):
download.file("http://pubs.usgs.gov/ds/425/USGS_DS_425_SHAPES.zip", "USGS_DS_425_SHAPES.zip")
system("7za e USGS_DS_425_SHAPES.zip")
system("ogrinfo Surficial_materials.shp")
surfgeo = read.dbf("Surficial_materials.dbf")
summary(surfgeo$dbf$UNIT_NAME)
surfgeo$dbf$UNIT_INT = as.integer(surfgeo$dbf$UNIT_NAME)
write.dbf(surfgeo, "Surficial_materials.dbf")
surfgeo.sum = data.frame(NAMES=levels(surfgeo$dbf$UNIT_NAME), Value=1:length(levels(surfgeo$dbf$UNIT_NAME)))
write.csv(surfgeo.sum, "surfacegeology_legend.csv")
unlink("Surficial_materials_a.shp")
system(paste('ogr2ogr -t_srs \"', proj4string(r), '\" Surficial_materials_a.shp Surficial_materials.shp'))
system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"Surficial_materials_a.shp\" -FIELD \"UNIT_INT\" -GRID \"surfgeo_100m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("GESUSG6.tif")
system(paste(gdal_translate, 'surfgeo_100m.sdat /data/GEOG/NASIS/covs100m/GESUSG6.tif -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\"'))

## SoilGrids250m - clay, sand, soil organic carbon, bulk density
sg.lst = sapply(c("BLDFIE","CLYPPT","SNDPPT","ORCDRC"), function(x){paste0("/data/GEOG/",x,"_M_sl", 1:7, "_250m_ll.tif")})
sg.out = gsub("250m_ll", "100m", basename(sg.lst))
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalwarp", "sg.lst", "sg.out", "cellsize", "r")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(sg.lst), function(i){ if(!file.exists(sg.out[i])){ system(paste0(gdalwarp, ' ', sg.lst[i], ' /data/GEOG/NASIS/covs100m/', sg.out[i], ' -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -ot \"Int16\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()


## Filter out missing pixels in parent material and drainage class maps:
pm = raster("/data/GEOG/NASIS/covs100m/PMTGSS7.tif")
rnd.pnts = sampleRandom(pm, size=12e4, na.rm=TRUE, sp=TRUE)
sel.tifs = paste0("/data/GEOG/NASIS/covs100m/", c("PMTGSS7.tif","DRNGSS7.tif","DEMNED6.tif","GESUSG6.tif","SLPNED6.tif","CRVNED6.tif","VBFNED6.tif","MRNNED6.tif","VDPNED6.tif","DVMNED6.tif","TPINED6.tif","NEGNED6.tif","POSNED6.tif","VDCNED6.tif","BARL10.tif"))
sfInit(parallel=TRUE, cpus=length(sel.tifs))
sfExport("rnd.pnts", "sel.tifs")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(sel.tifs, function(i){try( raster::extract(raster(i), rnd.pnts) )}) 
snowfall::sfStop()
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(sel.tifs))
## Convert to factors
PMTGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_pmaterial.csv")
DRNGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_drainage.csv")
#GEOUSG6.leg <- read.csv("geology_legend.csv")
GESUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/surfacegeology_legend.csv")
ov.lst$PMTGSS7.f = plyr::join(data.frame(VALUE=ov.lst$PMTGSS7), PMTGSS7.leg, type="left")$pmaterial_class_f
summary(ov.lst$PMTGSS7.f)
ov.lst$DRNGSS7.f = plyr::join(data.frame(Value=ov.lst$DRNGSS7), DRNGSS7.leg, type="left")$aggregate_class
summary(as.factor(ov.lst$DRNGSS7.f))
ov.lst$GESUSG6.f = plyr::join(data.frame(X=ov.lst$GESUSG6), GESUSG6.leg, type="left")$NAMES
str(ov.lst)
summary(ov.lst$GESUSG6.f)
GESUG6.df = as.data.frame(model.matrix( ~ GESUSG6.f - 1, data=ov.lst))
names(GESUG6.df) = make.names(names(GESUG6.df))
GESUG6.df <- GESUG6.df[match(rownames(ov.lst),rownames(GESUG6.df)),]
ov.lstF = cbind(ov.lst, GESUG6.df)
str(ov.lstF)

fm.pm <- as.formula(paste('PMTGSS7.f ~ DEMNED6 + SLPNED6 + CRVNED6 + VBFNED6 + MRNNED6 + VDPNED6 + DVMNED6 + TPINED6 + NEGNED6 + POSNED6 + VDCNED6 + BARL10 +', paste(names(GESUG6.df), collapse = "+"))) 
fm.dr <- as.formula(paste('DRNGSS7.f ~ DEMNED6 + SLPNED6 + CRVNED6 + VBFNED6 + MRNNED6 + VDPNED6 + DVMNED6 + TPINED6 + NEGNED6 + POSNED6 + VDCNED6 + BARL10 +', paste(names(GESUG6.df), collapse = "+"))) 
## Takes about 2-3 minutes to fit:
rf.pm <- ranger(fm.pm, data=ov.lstF[complete.cases(ov.lstF[,all.vars(fm.pm)]),], write.forest = TRUE, num.trees = 85)
rf.dr <- ranger(fm.dr, data=ov.lstF[complete.cases(ov.lstF[,all.vars(fm.dr)]),], write.forest = TRUE, num.trees = 85)
rf.pm ## 55% accuracy  
summary(rf.pm$predictions)
rf.dr ## 81% accuracy
## TH: This can be used to filter out all missing values
saveRDS.gz(rf.pm, file="RF_parent_material.rds")
saveRDS.gz(rf.dr, file="RF_drainage_class.rds")

## Test it:
m = readRDS("/data/tt/NASIS/covs100m/T4306/T4306.rds")
spplot(m["PMTGSS7"])
m$PMTGSS7.f = predict(rf.pm, m@data)$predictions
spplot(m["PMTGSS7.f"])
m$DRNGSS7.f = predict(rf.dr, m@data)$predictions
spplot(m["DRNGSS7.f"])

## Gamma radiometrics (https://mrdata.usgs.gov/radiometric/)
system('wget https://mrdata.usgs.gov/radiometric/data/NArad_K.zip')
system('wget https://mrdata.usgs.gov/radiometric/data/NArad_Th.zip')
system('wget https://mrdata.usgs.gov/radiometric/data/NArad_U.zip')
system('wget https://mrdata.usgs.gov/radiometric/data/NArad_exp.zip')
system('7za x NAMrad_K.zip')

gamma.lst = paste0("/data/USA48/", c("NAMrad_K.flt","NAMrad_Th.flt","NAMrad_U.flt","NAMrad_exp.flt"))
GDALinfo(gamma.lst[1])
tm.prj = "+proj=tmerc +lat_0=0 +lon_0=-100 +k=0.926 +ellps=clrk66 +datum=NAD27"
gamma.out = paste0("/data/GEOG/NASIS/covs100m/", c("GK0MRD3.tif","GTHMRD3.tif","GU0MRD3.tif","GEXMRD3.tif"))
sfInit(parallel=TRUE, cpus=4)
sfExport("gdalwarp", "gamma.lst", "gamma.out", "cellsize", "r", "tm.prj")
sfLibrary(raster)
out <- sfClusterApplyLB(1:length(gamma.lst), function(i){ if(!file.exists(gamma.out[i])){ system(paste0(gdalwarp, ' ', gamma.lst[i], ' ', gamma.out[i], ' -s_srs \"', tm.prj, '\" -t_srs \"', proj4string(r), '\" -r \"cubicspline\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -ot \"Int16\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))}} )
sfStop()

save.image()
