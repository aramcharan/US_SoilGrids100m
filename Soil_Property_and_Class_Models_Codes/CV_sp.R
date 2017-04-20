# Cross-validation of soil properties adapted from https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/cv/cv_functions.R
#

list.of.packages <- c("nnet", "plyr", "ROCR", "randomForest", "plyr", "parallel", "psych", "mda", "h2o", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


## --------------------------------------------------------------
## Build cross validation function for soil properties:
## --------------------------------------------------------------

cv_numeric <- function(formulaString, rmatrix, nfold, idcol, cpus, method="ranger", Log=FALSE){
  varn = all.vars(formulaString)[1]
  sel <- dismo::kfold(rmatrix, k=nfold)
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting method ", method," ...", sep=""))
  if(nfold > nrow(rmatrix)){
    stop("'nfold' argument must not exceed total number of points")
  }
  if(missing(cpus)){
    if(method=="randomForest"){
      cpus = nfold
    } else {
      cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE)
    }
  }
  if(method=="h2o"){
    out <- list()
    for(j in 1:nfold){
      out[[j]] <- predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, method=method, cpus=1)
    }
  }
  if(method=="caret"){
    out <- list()
    for(j in 1:nfold){
      out[[j]] <- predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, method=method, cpus=cpus)
    }
  }
  if(method=="ranger"){
    snowfall::sfInit(parallel=TRUE, cpus=cpus)
    snowfall::sfExport("predict_parallelP","idcol","formulaString","rmatrix","sel","varn","method")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="ranger", character.only=TRUE)
    out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, method=method, cpus=1)})
    snowfall::sfStop()
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE)
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  #R.squ ared = 1-sum((out$Observed - out$Predicted)^2, na.rm=TRUE)/(var(out$Observed, na.rm=TRUE)*sum(!is.na(out$Observed)))
  R.squared = 1-var(out$Observed - out$Predicted, na.rm=TRUE)/var(out$Observed, na.rm=TRUE)
  if(Log==TRUE){
    ## If the variable is log-normal then logR.squared is probably more correct
    logRMSE = sqrt(mean((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE))
    #logR.squared = 1-sum((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE)/(var(log1p(out$Observed), na.rm=TRUE)*sum(!is.na(out$Observed)))
    logR.squared = 1-var(log1p(out$Observed) - log1p(out$Predicted), na.rm=TRUE)/var(log1p(out$Observed), na.rm=TRUE)
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, R.squared=R.squared, logRMSE=logRMSE, logR.squared=logR.squared))
  } else {
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, R.squared=R.squared))
  }
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}


################################################
################################################
#set working directory with regression matrix
setwd("")

library(sp)
library(randomForest)
library(nnet)
library(plotKML)
library(GSIF)
library(plyr)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(hexbin)
library(gridExtra)
library(lattice)
library(grDevices)
library(h2o)
library(scales)
library(ranger)
library(xgboost)
library(caret)
library(doParallel)
library(RCurl)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)


## load the data
df <- readRDS("ovA_NCSS_peds.rds")

#### Extract names of environmental covariates ####
cov <- colnames(df)
cov <- cov[19:length(cov)]


# This code can be repeated for all soil properties - 2 soil properties provided here
####Clay####
formulaStringClay = as.formula(paste('clay ~', paste(cov, collapse="+")))
df.clay <- df[,all.vars(formulaStringClay)]
#### Add location identifier and build matrix for cv_numeric ####
df.clay$LOC_ID <- df$LOC_ID
sel.na.clay <- complete.cases(df.clay)
summary(sel.na.clay)
df.clay = df.clay[sel.na.clay,]
#### Run cross validation ####
clay.prop <- cv_numeric(formulaStringClay, rmatrix=df.clay, nfold=10, idcol="LOC_ID")
saveRDS(clay.prop, file = "clay.prop2.Rds")
str(clay.prop)
#Read saved cross validation results
#clay.prop <- readRDS("clay.prop2.Rds")


####SOC####
formulaStringSOC= as.formula(paste('soc ~ ', paste(cov, collapse="+")))
df.soc = df[,all.vars(formulaStringSOC)]
#### Add location identifier and build matrix for cv_numeric ####
df.soc$LOC_ID <- df$LOC_ID
df.soc$LATWGS84 <- df$LATWGS84
df.soc$LONWGS84 <- df$LONWGS84
sel.na.soc <- complete.cases(df.soc)
summary(sel.na.soc)
df.soc = df.soc[sel.na.soc,]

soc.prop <- cv_numeric(formulaStringSOC, rmatrix=df.soc, nfold=  10, idcol="LOC_ID", Log = TRUE)
saveRDS(soc.prop, file = "soc.prop3.Rds")
str(soc.prop)

#### Load Models and plot results ####
ph.prop <- readRDS("ph.prop2.Rds")
n_tot.prop <- readRDS("n_tot2.prop.Rds")
bd.prop <- readRDS("bd.prop2.Rds")
soc.prop <- readRDS("soc.prop3.Rds")
sand.prop <- readRDS("sand.prop1.Rds")
clay.prop <- readRDS("clay.prop1.Rds")

#### Plot CV results:####

plt_clay <- xyplot(clay.prop[[1]]$Predicted~clay.prop[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("#404788FF", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (machine learning)")


## Hexbin plot - Clay ##

d.meas.clay <- min(clay.prop[[1]]$Observed, na.rm=TRUE)
pred.clay <- clay.prop[[1]]$Predicted+ifelse(d.meas.clay==0, 1, d.meas.clay)
meas.clay <- clay.prop[[1]]$Observed+ifelse(d.meas.clay==0, 1, d.meas.clay)
lim.clay <- range(clay.prop[[1]]$Observed, na.rm=TRUE)

## Hexbin plot - SOC ##

d.meas.soc <- min(soc.prop[[1]]$Observed, na.rm=TRUE)
pred.soc <- soc.prop[[1]]$Predicted+ifelse(d.meas.soc==0, 1, d.meas.soc)
meas.soc <- soc.prop[[1]]$Observed+ifelse(d.meas.soc==0, 1, d.meas.soc)
lim.soc <- range(soc.prop[[1]]$Observed, na.rm=TRUE)
lim.soc[1] <- ifelse(d.meas.soc==0,lim.soc[1]<-1, lim.soc[1])


# Set color scheme
viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF")

plt.clay <- hexbinplot(pred.clay~meas.clay,colramp=colorRampPalette(rev(viri)), main="Percent Clay", xlab="Measured", ylab="Predicted", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, asp=1, xbins=25, density=40, xlim=lim.clay, ylim=lim.clay, panel=pfun)

plt.soc <- hexbinplot(pred.soc~meas.soc,colramp=colorRampPalette(rev(viri)), main="log(Percent SOC)", xlab="Measured", ylab="Predicted", type="g", scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), lwd=1, lcex=8, inner=.2, cex.labels=.8, asp=1, xbins=25, density=40, xlim=lim.soc, ylim=lim.soc, panel=pfun)

# log plot for n and soc.

#### ggplots with log scale for colors, along with perfect fit and best fit (least square method) ####

my_breaks <- c(2,10,50,200,1000,6000)

# Clay
df_clay <- clay.prop[[1]]
clay.coef <- coef(lm(df_clay$Observed ~ df_clay$Predicted))

plt.clay.2 <-   ggplot(clay.prop[[1]], aes(Observed, Predicted))  + stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1)   + xlim(lim.clay) + ylim(lim.clay)   + theme(axis.text=element_text(size=6), legend.text=element_text(size=6),axis.title=element_text(size=6), plot.title=element_text(size=6), legend.title=element_text(size=6))  + xlab("Measured") + ylab("Predicted")  + scale_fill_gradientn(name = "log(Count)", trans = "log", breaks = my_breaks, labels = my_breaks, colours = rev(plasma(256)))+ geom_abline(intercept = clay.coef[1], slope = clay.coef[2], linetype = 2) + ggtitle("Clay % (CV"~R^{2}~"=0.72)")

# SOC
df_soc <- soc.prop[[1]]
df_soc$log_obs <- log10(df_soc$Observed)
df_soc$log_pred <- log10(df_soc$Predicted)
df_soc <- subset(df_soc,df_soc$log_obs < 3) 
df_soc <- subset(df_soc,df_soc$log_obs > -3)
soc.coef <- coef(lm(df_soc$log_obs ~ df_soc$log_pred))


plt.soc.2 <-   ggplot(df_soc, aes(log_obs, log_pred)) + stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1)  + xlim(-3,2) + ylim(-3,2) + theme(axis.text=element_text(size=6), legend.text=element_text(size=6),axis.title=element_text(size=6), plot.title=element_text(size=6)) + xlab("Measured") + ylab("Predicted") + scale_fill_gradientn(name = "log(Count)", trans = "log", breaks = my_breaks, labels = my_breaks, colours = rev(plasma(256))) + ggtitle("log(SOC) (CV"~R^{2}~"=0.63)")+ geom_abline(intercept = soc.coef[1], slope = soc.coef[2], linetype = 2)

# Create shared legend for multiple plots
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

#plot panel plots of correlation graphs
corplts<- grid_arrange_shared_legend(plt.clay.2, plt.sand.2,plt.ph.2, plt.bd.2, plt.soc.2,plt.n_tot.2, ncol = 3, nrow = 2)
ggsave('correlation5.png', plot = corplts, device = "png", dpi = 600, limitsize = TRUE)








