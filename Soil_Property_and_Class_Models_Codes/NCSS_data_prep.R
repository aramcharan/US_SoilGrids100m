#### Prepare soil properties for soilgrids : soc,bd, sand, silt, clay, pH_h2o, n_tot ####


rm(list = ls())

library(plyr)
#setwd("C:/Users/amr418/Google Drive/Dissertation Research/Final DATA")


#### Extract csv files from Microsoft Access ####
Peds_MA <- read.csv("Q_SOC_BD_samplelist12.20.2016.csv")
Salinity <- read.csv("Salt12.20.2016.csv")
ph <- read.csv("pH_and_Carbonates.csv")
Peds <- join(Peds_MA,Salinity, type = "left", match="first")
Peds <- join(Peds,ph, by = "labsampnum", type = "left", match="first")

####Format and clean dataframe####

Peds <-unique(Peds, incomparables = FALSE)
names(Peds) <- tolower(names(Peds))

#Edit variable names 
names(Peds)[names(Peds)== "latitude_decimal_degrees"]<- "lat"
names(Peds)[names(Peds)== "longitude_decimal_degrees"]<- "long"
names(Peds)[names(Peds)== "layer_sequence"]<- "hzn_seq"

names(Peds)[names(Peds)== "carb"]<- "carbonate"
names(Peds)[names(Peds)== "hor_thick"]<- "hzn_thick"
names(Peds)[names(Peds)== "model_desg"]<- "hzn_desg"
names(Peds)[names(Peds)== "tex_psda"]<- "text_class"
names(Peds)[names(Peds)== "vfg2"]<- "rock_percent"
names(Peds)[names(Peds)== "avgofsand_tot_psa"]<- "sand"
names(Peds)[names(Peds)== "avgofsilt_tot_psa"]<- "silt"
names(Peds)[names(Peds)== "avgofclay_tot_psa"]<- "clay"

#### Correct missing data and zeros in data ####

Peds$cn_ratio[Peds$cn_ratio > 600 | Peds$cn_ratio == 0] <- NA    # C/N ratio of spruce sawdust Brady & Weil p. 329 
Peds$n_tot[Peds$n_tot == 0] <- NA           # %N of spruce sawdust Brady & Weil p.329

Peds$bd[Peds$bd == 0] <- NA                 # Indicate empty data cells
Peds$bd[Peds$bd > 2.65] <- NA               # BD higher that quartz mineral is incorrect Brady & Weil p.108

Peds$silt[Peds$silt < 0] <- NA
Peds$clay[Peds$clay < 0] <- NA
Peds$sand[Peds$sand < 0] <- NA

Peds$clay[Peds$clay == 0 & Peds$sand == 0] <- NA
Peds$sand <- ifelse(is.na(Peds$clay), NA, Peds$sand)
Peds$silt <- ifelse(is.na(Peds$clay), NA, Peds$silt)

Peds$rock_percent[Peds$rock_percent >= 100] <- NA           # Replace instances where rock volume content is more than 100%
Peds$soc[Peds$soc < 0] <- NA

#### Remove data that has no lat long info ####
Peds <- subset(Peds, !is.na(Peds$lat))


#### Save file ####
saveRDS(Peds, file = "NCSS_Peds12202016.Rds")


