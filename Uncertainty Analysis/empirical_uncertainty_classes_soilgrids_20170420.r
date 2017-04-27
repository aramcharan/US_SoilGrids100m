## Script converts prediction probabilities (e.g. from a random forest)
## to validation probabilities by relating validation match rates
## to the prediction probabilities. THe assumption is that the
## prediction probabilites reflect some kind of certainty in a prediction
## but must be validated and calibrated by an independent sample
## Credits: Travis Nauman, 4/20/2017

## Libraries
require(ggplot2)

## Import table with validation and prediction values along with prediction probabilities
setwd("C:/models/SoilGrids100m/forGithub")
## Table of independent field observations of soil classes
## must have attributes with field observed soil class (e.g. $pscsmodorg) and the map predicted soil class (e.g. $pscs_cn1)
## as well as the prediction probabilities for the mapped classes ($class_pprob)
pts = read.delim("uncert_pedons_forR_ttab_20170413.txt", header=TRUE)

## Pull out any NA or empty values: gg
ptsgg = subset(pts, pts$gg_c1n != "NA")
ptsgg = subset(pts, pts$gg_c1n != "")
ptsgg = subset(pts, pts$gg != "NA")
ptsgg = subset(pts, pts$gg != "")

## Pull out any NA or empty values: mPSC
pts = subset(pts, pts$pscs_c1n != "NA")
pts = subset(pts, pts$pscs_c1n != "")
pts = subset(pts, pts$pscsmodorg != "NA")
pts = subset(pts, pts$pscsmodorg != "")

## Create nominal categories to split distibution into strata to calculate val probs
## Number of splits will depend on the size of the validation set. Good to have at least 10
## instances in each strata. Much more is better
nbins = 12
nbinsgg = 10

## Use ggplot to cut sample into bins of ~equal sample size
pts$pbins = as.factor(as.numeric(cut_number(pts$pscs_pprob,nbins)))
ptsgg$pbins = as.factor(as.numeric(cut_number(ptsgg$gg_pprob,nbinsgg)))

## Tally matches of map predictions to field observations in each probability bin
pts$pscsmatch = ifelse(as.character(pts$pscsmodorg) == as.character(pts$pscs_c1n), 1, 0)
ptsgg$ggmatch = ifelse(as.character(ptsgg$gg) == as.character(ptsgg$gg_c1n), 1, 0)

## Create Index for rebinning
bindex = seq.int(1,nbins)
bindexgg = seq.int(1,nbinsgg)

## Create new DF for summarized val probs:mPSC
vprobdf = data.frame(bindex)
vprobdf$pprobs_mPSC = "na"
vprobdf$vprobs_mPSC = "na"

## Create new DF for summarized val probs:gg
vprobdfgg = data.frame(bindexgg)
vprobdfgg$pprobs_gg = "na"
vprobdfgg$vprobs_gg = "na"

## Loop through to subset and calculate new probs for mPSC
for(i in bindex){
  bindf = subset(pts, pts$pbins == i)
  vprob_mPSC = sum(bindf$pscsmatch)/length(bindf$pscsmatch)
  pprob_mPSC = median(bindf$pscs_pprob)
  vprobdf$pprobs_mPSC[vprobdf$bindex == i] = as.numeric(pprob_mPSC)
  vprobdf$vprobs_mPSC[vprobdf$bindex == i] = as.numeric(vprob_mPSC)
  rm(bindf)
  "done with"
  print(i)
}

## Loop through to subset and calculate new probs for gg
for(i in bindexgg){
  bindf = subset(ptsgg, ptsgg$pbins == i)
  vprob_gg = sum(bindf$ggmatch)/length(bindf$ggmatch)
  pprob_gg = median(bindf$gg_pprob)
  vprobdfgg$pprobs_gg[vprobdfgg$bindex == i] = as.numeric(pprob_gg)
  vprobdfgg$vprobs_gg[vprobdfgg$bindex == i] = as.numeric(vprob_gg)
  rm(bindf)
  "done with"
  print(i)
}

## Clean up mPSC
vprobdf$vprobs_mPSC = as.numeric(vprobdf$vprobs_mPSC) * 100
vprobdf$pprobs_mPSC = as.numeric(vprobdf$pprobs_mPSC)
## Clean up gg 
vprobdfgg$vprobs_gg = as.numeric(vprobdfgg$vprobs_gg) * 100
vprobdfgg$pprobs_gg = as.numeric(vprobdfgg$pprobs_gg)

## Plot relationships
par(mfrow=c(1,2))
plot(vprobdfgg$vprobs_gg~vprobdfgg$pprobs_gg, ylim = c(0,100), xlim = c(0,100), main = "GG", xlab = "Prediction probabilities", ylab = "Validation Probabilities")

## Model fits and raster rendering
library(raster)
## max prob rasters must be named pprob_#class.tif
rasters=stack(list.files(getwd(),pattern=".tif$",full.names=FALSE))
x <- seq(1,100)
x1 <- seq(1,100)
y1 <- seq(1,100)

#Linear-log regression: GG
fitloggg = lm(vprobs_gg~log(pprobs_gg), data = vprobdfgg)
yfitloggg = predict(fitloggg, data.frame(pprobs_gg = x))
lines(x, yfitloggg, col = 'blue')
lines(x1,y1)
## To predict with GG prediction probability map (pprobs_gg.tif)
#predict(rasters, fitloggg, type="response",index=2,na.rm=TRUE,progress="window",dataType = 'INT1U',overwrite=TRUE,filename="valprobs_GG.tif")

## Linear-log regression: mPSC
plot(vprobdf$vprobs_mPSC~vprobdf$pprobs_mPSC, ylim = c(0,100), xlim = c(0,100), main = "mPSC", xlab = "Prediction probabilities", ylab = "Validation Probabilities")
fitlog = lm(vprobs_mPSC~log(pprobs_mPSC), data = vprobdf)
yfitlog = predict(fitlog, data.frame(pprobs_mPSC = x))
lines(x, yfitlog, col = 'blue')
lines(x1,y1)
## To predict with mPSC prediction probability map (pprob_mPSC.tif)
#predict(rasters, fitlog, type="response",index=2,na.rm=TRUE,progress="window",dataType = 'INT1U',overwrite=TRUE,filename="valprobs_mPSC.tif")
