
# Creates pup metric summaries by complex and year for the diet analysis

library(here)
library(tidyverse)
library(readxl)
library(corrr)

# Set paths and read in data ----------------------------------------------

input.path<-here("data-raw")
output.path<-here("data")

pupWeight<-data.table::fread(file.path(input.path,"pupWTLiz_2023.csv"))
pupMort<-data.table::fread(file.path(input.path,"liz_mort_2023.csv"))

complexes<-readxl::read_xlsx(file.path(input.path,"Rookery complexes.xlsx")) |>
  group_by(rkname)|>
  dplyr::summarise(Complex=Complex[1])

complexes2<-readxl::read_xlsx(file.path(input.path,"Rookery complexes.xlsx")) 

# Pup weight --------------------------------------------------------------

pupWeight2<-pupWeight |>
  mutate(rookery=str_trim(rookery))|>
  left_join(complexes, by=c("rookery"="rkname"))|>
  mutate(Date=as.POSIXct(paste("2021",month,day,sep="-")),Complex=factor(Complex),
         Sex=factor(sex,labels=c("U","F","F","M","U")), length=ifelse(plen<0, NA, plen),
         lengthF=ifelse(plen<0 | Sex=="M",NA, plen),lengthM=ifelse(plen<0 | Sex=="F",NA, plen))|>
  filter((Sex=="M" | Sex=="F") & month==8 & !is.na(weight) & weight>2) |>
  group_by(Complex,Sex)|>
  nest()|>
  mutate(lm=map(data, ~lm(weight ~ day, data=.x) |> broom::tidy())) |>
  unnest(lm) |>
  unnest(data) |>
  filter(term=="day") |>
  mutate(weightAdj=weight+((25-day)*estimate))|>
  dplyr::select(Complex,rookery,Sex,weight,weightAdj,day,year)|>
  dplyr::rename("Year"="year")|>
  filter(day>=19 & day<=31)

pupWeightSum<-pupWeight2 |>
  group_by(Year,Complex,Sex)|>
  dplyr::summarise(mWeight=mean(weightAdj),mWeightUnAdj=mean(weight),medWeight=median(weightAdj),qWeight25=quantile(weightAdj, 0.25),qWeight75=quantile(weightAdj,0.75),sdWeight=sd(weightAdj),sdWeightUnAdj=sd(weight), mDay=mean(day), nPup=length(!is.na(weightAdj)))
 
pupWeightSumRook<-pupWeight2 |>
  group_by(Year,Complex,rookery,Sex)|>
  dplyr::summarise(mWeight=mean(weightAdj),sdWeight=sd(weightAdj), nPup=length(weightAdj))

# Save output
saveRDS(pupWeightSum,file.path(output.path, "Pup weights by complex.rds"))
saveRDS(pupWeightSumRook, file.path(output.path, "Pup weights by rookery.rds"))

# Pup mortality -----------------------------------------------------------

pupMort2<-pupMort|>
  filter(!is.na(deadPups))|>
  dplyr::left_join(complexes,by=c("rcod"="rkname"))|>
  dplyr::rename("Year"="year")|>
  dplyr::mutate(MMort=deadPups/pupsBorn)
  
pupMortSum<-pupMort2|>
  group_by(Year,Complex)|>
  dplyr::summarise(MMort=mean(deadPups/pupsBorn,na.rm=T))|>
  ungroup()|>
  filter(!is.na(Complex) & !is.nan(MMort))

# Save output
saveRDS(pupMortSum, file.path(output.path, "Pup mortality by complex.rds"))
saveRDS(pupMort2, file.path(output.path, "Pup mortality by rookery.rds"))
