## Cross-validation of soil properties adapted from https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/cv/cv_functions.R
## Tom.hengl@gmail.com and Amanda Ramcharan <a.m.ramcharan@gmail.com>

list.of.packages <- c("nnet", "plyr", "ROCR", "randomForest", "plyr", "parallel", "psych", "mda", "h2o", "dismo", "grDevices", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


## --------------------------------------------------------------
## Build cross validation function for soil properties:
## --------------------------------------------------------------

## predict soil properties in parallel:
predict_parallelP <- function(j, sel, varn, formulaString, rmatrix, idcol, method, cpus, Nsub=1e4, remove_duplicates=FALSE){
  s.train <- rmatrix[!sel==j,]
  if(remove_duplicates==TRUE){
    ## TH: optional - check how does model performs without the knowledge of the 3D dimension
    sel.dup = !duplicated(s.train[,idcol])
    s.train <- s.train[sel.dup,]
  }
  s.test <- rmatrix[sel==j,]
  n.l <- dim(s.test)[1]
  if(missing(Nsub)){ Nsub = length(all.vars(formulaString))*50 }
  if(Nsub>nrow(s.train)){ Nsub = nrow(s.train) }
  if(method=="h2o"){
    ## select only complete point pairs
    train.hex <- as.h2o(s.train[complete.cases(s.train[,all.vars(formulaString)]),all.vars(formulaString)], destination_frame = "train.hex")
    gm1 <- h2o.randomForest(y=1, x=2:length(all.vars(formulaString)), training_frame=train.hex) 
    gm2 <- h2o.deeplearning(y=1, x=2:length(all.vars(formulaString)), training_frame=train.hex)
    test.hex <- as.h2o(s.test[,all.vars(formulaString)], destination_frame = "test.hex")
    v1 <- as.data.frame(h2o.predict(gm1, test.hex, na.action=na.pass))$predict
    gm1.w = gm1@model$training_metrics@metrics$r2
    v2 <- as.data.frame(h2o.predict(gm2, test.hex, na.action=na.pass))$predict
    gm2.w = gm2@model$training_metrics@metrics$r2
    ## mean prediction based on accuracy:
    pred <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
    gc()
    h2o.removeAll()
  }
  if(method=="caret"){
    test = s.test[,all.vars(formulaString)]
    ## tuning parameters:
    cl <- makeCluster(cpus)
    registerDoParallel(cl)
    ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
    gb.tuneGrid <- expand.grid(eta = c(0.3,0.4), nrounds = c(50,100), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
    rf.tuneGrid <- expand.grid(mtry = seq(4,length(all.vars(formulaString))/3,by=2))
    ## fine-tune RF parameters:
    t.mrfX <- caret::train(formulaString, data=s.train[sample.int(nrow(s.train), Nsub),], method="rf", trControl=ctrl, tuneGrid=rf.tuneGrid)
    gm1 <- ranger(formulaString, data=s.train, write.forest=TRUE, mtry=t.mrfX$bestTune$mtry)
    gm1.w = 1/gm1$prediction.error
    gm2 <- caret::train(formulaString, data=s.train, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid)
    gm2.w = 1/(min(gm2$results$RMSE, na.rm=TRUE)^2)
    v1 <- predict(gm1, test, na.action=na.pass)$predictions
    v2 <- predict(gm2, test, na.action=na.pass)
    pred <- rowSums(cbind(v1*gm1.w, v2*gm2.w))/(gm1.w+gm2.w)
  }
  if(method=="ranger"){
    gm <- ranger(formulaString, data=s.train, write.forest=TRUE, num.trees=85)
    pred <- predict(gm, s.test, na.action = na.pass)$predictions 
  }
  obs.pred <- as.data.frame(list(s.test[,varn], pred))
  names(obs.pred) = c("Observed", "Predicted")
  obs.pred[,idcol] <- s.test[,idcol]
  obs.pred$fold = j
  return(obs.pred)
}

cv_numeric <- function(formulaString, rmatrix, nfold, idcol, cpus, method="ranger", Log=FALSE, LLO=TRUE){     
  varn = all.vars(formulaString)[1]
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting method ", method," ...", sep=""))
  if(nfold > nrow(rmatrix)){ 
    stop("'nfold' argument must not exceed total number of points") 
  }
  if(sum(duplicated(rmatrix[,idcol]))>0.5*nrow(rmatrix)){
    if(LLO==TRUE){
      ## TH: Leave whole locations out
      ul <- unique(rmatrix[,idcol])
      sel.ul <- dismo::kfold(ul, k=nfold)
      sel <- lapply(1:nfold, function(o){ data.frame(row.names=which(rmatrix[,idcol] %in% ul[sel.ul==o]), x=rep(o, length(which(rmatrix[,idcol] %in% ul[sel.ul==o])))) })
      sel <- do.call(rbind, sel)
      sel <- sel[order(as.numeric(row.names(sel))),]
      message(paste0("Subsetting observations by unique location"))
    } else {
      sel <- dismo::kfold(rmatrix, k=nfold, by=rmatrix[,idcol])
      message(paste0("Subsetting observations by '", idcol, "'"))
    }
  } else {
    sel <- dismo::kfold(rmatrix, k=nfold)
    message(paste0("Simple subsetting of observations using kfolds"))
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
    snowfall::sfInit(parallel=TRUE, cpus=ifelse(nfold>cpus, cpus, nfold))
    snowfall::sfExport("predict_parallelP","idcol","formulaString","rmatrix","sel","varn","method")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="ranger", character.only=TRUE)
    out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, method=method)})
    snowfall::sfStop()
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE) 
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  #R.squared = 1-sum((out$Observed - out$Predicted)^2, na.rm=TRUE)/(var(out$Observed, na.rm=TRUE)*sum(!is.na(out$Observed)))
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
  closeAllConnections()
}

## correlation plot:
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.abline(0,1,lty=1,lw=2,col="black")
  panel.abline(0+RMSE,1,lty=2,lw=2,col="black")
  panel.abline(0-RMSE,1,lty=2,lw=2,col="black")
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








