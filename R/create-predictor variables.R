# This script downloads environmental variables of interest and links it with
# complexes to create the predictor variables for the analysis

#devtools::install_github("sean-rohan-NOAA/akgfmaps", build_vignettes = TRUE)
#devtools::install_github("afsc-gap-products/coldpool")

library(data.table)
library(stringr)
library(coldpool)
library(rerddap)
library(tidyverse)
library(sf)
library(readxl)
library(here)
library(dataRetrieval)
library(stars)
library(akgfmaps)
library(ggplot2)

foYears<-1987:2022
foMonths<-7:11
input.path<-file.path(here(), "data-raw")
afpath1<-file.path(here(),"data","Adult female 95% UD.rds")
coords<-readRDS(afpath1) |>
  sf::st_transform(crs=4326)|>
  st_bbox()

#-----------------------------------------------------------------------------------
# 1. Read in environmental data of interest
#-----------------------------------------------------------------------------------

# Temperature data from survey --------------------------------------------
bottomTempS<-coldpool::ebs_bottom_temperature |>
  terra::rast()

names(bottomTempS)<-parse_number(names(bottomTempS))

surfaceTempS<-ebs_surface_temperature |>
  terra::rast()

names(surfaceTempS)<-parse_number(names(surfaceTempS))

# Temperature data from remote sensing --------------------------------------------

oisst<-read_ncdf("http://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/sst.mon.mean.nc",var="sst",
                 ncsub = cbind(start = c(610, 570, 72), 
                               count = c(200, 40,421))) 

dates<-oisst |>
  st_get_dimension_values("time") |>
  as.data.frame()|>
  data.table::setnames("time")|>
  filter(month(time)>=7 & month(time)<=11)

datesIndex<-which(st_get_dimension_values(oisst, "time")  %in% dates$time)

oisst2<-oisst |>
  slice(time,datesIndex)|>
  as("SpatRaster")

# Fish data from survey (from https://www.fisheries.noaa.gov/foss)  -----
fish<-data.table::fread(file.path(input.path,"EBS Trawl Survey Data.csv")) |>
  dplyr::select(Year, `Common Name`,`Longitude Dd`,`Latitude Dd`,`Weight Kg`,Count,`Area Swept Ha`,`Cpue Kgkm2`,`Cpue Nokm2`) |>
  filter(`Common Name`=="Walleye Pollock" | `Common Name`=="Pacific Herring" | `Common Name` == "Atka Mackerel" | `Common Name` ==
           "Sablefish")

#-----------------------------------------------------------------------------------
# 2. Link local variables with each complex
#-----------------------------------------------------------------------------------

af95<-readRDS(afpath1) |>
  sf::st_transform(crs=4326) 

af950360<-af95 |>
  st_shift_longitude()

# Clip to just the shelf
bsSurvey <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = st_crs(af95))$survey.area

af95Shelf<-st_intersection(af95, bsSurvey) |>
  st_shift_longitude()

# Remotely sensed data ----------------------------------------------------

oisstExtract<-list()
oisstExtractShelf<-list()

oisstMonth<-str_extract_all(names(oisst2), "\\d+", simplify=T)[,2]
oisstYear<-str_extract_all(names(oisst2), "\\d+", simplify=T)[,1]

for (i in 1:length(names(oisst2))){
  
  oisstExtract[[i]]<- exactextractr::exact_extract(oisst2[[i]], af950360,"mean") |>
    data.frame(Complex=af950360$id, var="sstRe", Month=oisstMonth[i], Year=oisstYear[i]) %>%
    data.table::setnames(1,"value")
  
  oisstExtractShelf[[i]]<- exactextractr::exact_extract(oisst2[[i]], af95Shelf,"mean") |>
    data.frame(Complex=af95Shelf$id, var="sstReShelf", Month=oisstMonth[i], Year=oisstYear[i]) %>%
    data.table::setnames(1,"value")
}

sstRe<-data.table::rbindlist(oisstExtract)
sstReShelf<-data.table::rbindlist(oisstExtractShelf)
remoteSST<-rbind(sstRe, sstReShelf) 

# Survey data - temperature -----------------------------------------------

sstExtract<-list()
btExtract<-list()
surfaceTempS<-project(surfaceTempS,af95)
bottomTempS<-project(bottomTempS, af95)

for (i in 1:length(names(bottomTempS))){

  sstExtract[[i]]<- exactextractr::exact_extract(surfaceTempS[[i]], af95,"mean") %>%
    data.frame(Complex=af95$id, var="sstS", Year=names(surfaceTempS)[i]) %>%
    data.table::setnames(1,"value")
  
  btExtract[[i]]<- exactextractr::exact_extract(bottomTempS[[i]], af95,"mean") %>%
    data.frame(Complex=af95$id, var="bottS", Year=names(bottomTempS)[i]) %>%
    data.table::setnames(1,"value")
}

sstS<-data.table::rbindlist(sstExtract) |>
  filter(Year %in% foYears)

bottS<-data.table::rbindlist(btExtract)|>
  filter(Year %in% foYears)

# Survey data - fish ------------------------------------------------------

# Convert to sf 
fishsp<-st_as_sf(fish, coords = c("Longitude Dd", "Latitude Dd"), crs = 4326) 

# Join with complex UD
fishsp2 <- st_join(fishsp, af95, join = st_within)

# Create summaries
fishSum <-fishsp2 |>
  st_drop_geometry() |>
  group_by(id, Year, `Common Name`) |>
  dplyr::summarise(MeanCPUEkg=mean(`Cpue Kgkm2`, na.rm=T), MeanCPUENo=mean(`Cpue Nokm2`, na.rm=T),
                   SumCPUEkg1=(sum(`Weight Kg`)/sum(`Area Swept Ha`))*100,SumCPUENo1=(sum(`Count`)/sum(`Area Swept Ha`))*100,
                   SumCPUEkg2=sum(`Cpue Kgkm2`, na.rm=T), SumCPUENo2=sum(`Cpue Nokm2`, na.rm=T)) |>
  filter(Year %in% foYears)|>
  pivot_longer(-c(id, Year, `Common Name`),values_to="value", names_to="variable") |>
  unite("var",`Common Name`:variable) |>
  data.table::setnames('id',"Complex")

# Put all values that are static for the whole season
surveyVars<-plyr::rbind.fill(fishSum, bottS, sstS)|>
  filter(!is.na(Complex))

#-----------------------------------------------------------------------------------
# 4. Save output
#-----------------------------------------------------------------------------------

fwrite(surveyVars,file.path(here(),"data", "Independent predictor variables - Survey.csv"))
fwrite(remoteSST, file.path(here(),"data", "Independent predictor variables - Remote SST.csv"))
