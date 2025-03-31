
library(data.table)
library(tidyverse)

output.path<-file.path(here::here(), "output")
input.path<-file.path(here::here(), "data")

# Read in response and predictor variables ---------------------------------

# Diet data
diet<-fread(file.path(input.path, "Diet response variables for analysis_Focal prey species.csv"))|>
  mutate(FO=freq/N)

# Pup metrics
pupWeight<-fread(file.path(input.path,"pupWTLiz_2023.csv"))
pupMort<-fread(file.path(input.path,"liz_mort_2023.csv"))
complexes<-readxl::read_xlsx(file.path(input.path,"Rookery complexes.xlsx")) |>
  group_by(rkname)|>
  dplyr::summarise(Complex=Complex[1])
complexes2<-readxl::read_xlsx(file.path(input.path,"Rookery complexes.xlsx")) 


# Create linked data for each complex -------------------------------------

dietVars<- diet |>
  mutate(Complex=as.factor(Complex), KLPreyGroup=as.factor(KLPreyGroup), Island=as.factor(Island)) |>
  group_by(Complex, KLPreyGroup)|>
  dplyr::mutate(mFO=mean(FO), anamFO=FO-mFO, sdFO=sd(FO))|>
  ungroup()|>
  tidytable::group_split(Complex, .named=T)|>
  map(droplevels)

dietVarsUn<-rbindlist(dietVars)


# Reorganize and filter - pup weight ----------------------------------------------

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
  dplyr::select(Complex,Sex,weight,weightAdj,day,year)|>
  dplyr::rename("Year"="year")|>
  filter(day>=19 & day<=31)
  
  pupWeightSum<-pupWeight2 |>
  group_by(Year,Complex,Sex)|>
  dplyr::summarise(mWeight=mean(weightAdj),mWeightUnAdj=mean(weight),medWeight=median(weightAdj),qWeight25=quantile(weightAdj, 0.25),qWeight75=quantile(weightAdj,0.75),sdWeight=sd(weightAdj),sdWeightUnAdj=sd(weight), mDay=mean(day), nPup=length(!is.na(weightAdj)))|>
  ungroup() 
  
# Reorganize and filter - pup mortality ------------------------------------

pupMortSum<-pupMort|>
  filter(!is.na(deadPups))|>
  dplyr::left_join(complexes,by=c("rcod"="rkname"))|>
  dplyr::rename("Year"="year")|>
  group_by(Year,Complex)|>
  dplyr::summarise(MMort=mean(deadPups/pupsBorn,na.rm=T))|>
  ungroup()|>
  filter(!is.na(Complex) & !is.nan(MMort))

# Join with diet data -----------------------------------------------------

# Pup weight calculated as a mean - has multiple joins because of pup sex
dietWeight<-dietVarsUn |>
  left_join(pupWeightSum) |>
  dplyr::select(Complex, Year, Month, FO,anamFO,sdFO,KLPreyGroup, Sex, mWeight, medWeight,qWeight25,qWeight75,mWeightUnAdj,sdWeight, sdWeightUnAdj,sdFO, nPup, mDay) |>
  filter(!is.na(mWeight)) 

dietWeightW<-dietWeight |>
  pivot_wider(names_from=KLPreyGroup, values_from=FO,id_cols = c(Complex,Year,Sex,mWeight,medWeight,mWeightUnAdj,qWeight25,qWeight75,nPup, mDay,sdWeight, sdWeightUnAdj))|>
  mutate(Complex=as.factor(Complex))

# Prey groups of Pollock and not pollock - sum
dietWeightW2<-dietWeight |>
  mutate(KLPreyGroup2=ifelse(KLPreyGroup=="Pollock", "Pollock", "Other"))|>
  group_by(Complex, Year, Sex,KLPreyGroup2,mWeight,anamWeight, sdWeight)|>
  summarise(sumFO=sum(FO))|>
  pivot_wider(names_from=KLPreyGroup2, values_from=sumFO,id_cols = c(Complex,Year,Sex,mWeight,anamWeight,sdWeight))|>
  mutate(Complex=as.factor(Complex), Exclude=case_when(Year==1994 & (Complex=="SGNorth" | Complex=="SGSouth")~"Yes",.default="No"))

#Mortality
dietMort<-dietVarsUn |>
  left_join(pupMortSum)|>
  dplyr::select(Complex, Year,FO,KLPreyGroup,MMort,sdFO) |>
  filter(!is.na(MMort))

dietMortW<-dietMort |>
  pivot_wider(names_from=KLPreyGroup, values_from=FO,  id_cols = c(Complex,Year,MMort)) |>
  mutate(Complex=as.factor(Complex), YearComplex=interaction(Year,Complex))

dietMortW2<-dietMort |>
  mutate(KLPreyGroup2=ifelse(KLPreyGroup=="Pollock", "Pollock", "Other"))|>
  group_by(Complex, Year,MMort, KLPreyGroup2)|>
  summarise(sumFO=sum(FO))|>
  pivot_wider(names_from=KLPreyGroup2, values_from=sumFO,id_cols = c(Complex,Year,MMort))|>
  mutate(Complex=as.factor(Complex))


# Analysis ----------------------------------------------------------------

library(mgcv)
library(gratia)
library(DHARMa)

# Pup weights - mean --------------------------------------------------------
fitWeight1<-gam(mWeight~Sex+s(Complex, bs="re")+s(Year)+s(Other)+s(Pollock), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")
fitWeight2<-gam(mWeight~Sex+s(Complex, bs="re")+s(Pollock) + s(Year), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")
fitWeight3<-gam(mWeight~Sex+s(Complex, bs="re")+s(Year) + s(Other), data=subset(dietWeightW2, Exclude=="No"),select=T, method="REML")
fitWeight4<-gam(mWeight~Sex+s(Complex, bs="re")+s(Pollock) + s(Other), data=subset(dietWeightW2, Exclude=="No"),select=T, method="REML")

gam.check(fitWeight1)
appraise(fitWeight1)

concurvity(fitWeight1)
concurvity(fitWeight2)
concurvity(fitWeight3)
concurvity(fitWeight4)

# To see how smooths change
draw(compare_smooths(fitWeight1,fitWeight2,fitWeight3,fitWeight4))

# With year
pupWeightResid1<-simulateResiduals(fitWeight1, plot=F)
plot(pupWeightResid1)
plotResiduals(pupWeightResid1, form =dietWeightW2$Year[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW2$Pollock[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW2$Other[dietWeightW2$Exclude=="No"])

# Without year
pupWeightResid2<-simulateResiduals(fitWeight4, plot=F)
plotResiduals(pupWeightResid2, form =dietWeightW2$Year[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW2$Pollock[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW2$Other[dietWeightW2$Exclude=="No"])

# Pup weights sd --------------------------------------------------------

fitWeightsd1<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Pollock)+s(Year)+s(Other), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")
fitWeightsd2<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Pollock)+s(Other), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")
fitWeightsd3<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Year)+s(Other), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")
fitWeightsd4<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Year)+s(Pollock), data=subset(dietWeightW2, Exclude=="No"), select=T, method="REML")

gam.check(fitWeightsd1)
gam.check(fitWeightsd2)
gam.check(fitWeightsd3)

concurvity(fitWeightsd1)
concurvity(fitWeightsd2)
concurvity(fitWeightsd3)

appraise(fitWeightsd1)
appraise(fitWeightsd2)
appraise(fitWeightsd3)

draw(compare_smooths(fitWeightsd1,fitWeightsd2,fitWeightsd3,fitWeightsd4))

pupWeightSDResid1<-simulateResiduals(fitWeightsd1, plot=F)
pupWeightSDResid2<-simulateResiduals(fitWeightsd2, plot=F)

plot(pupWeightSDResid1)
plotResiduals(pupWeightSDResid1, form =dietWeightW2$Year[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW2$Other[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW2$Pollock[dietWeightW2$Exclude=="No"])

plot(pupWeightSDResid2)
plotResiduals(pupWeightSDResid2, form =dietWeightW2$Year[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW2$Other[dietWeightW2$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW2$Pollock[dietWeightW2$Exclude=="No"])


# Pup Mortality --------------------------------------------------------
fitMort1<-gam(MMort~s(Complex, bs="re")+s(Pollock)+s(Other)+s(Year), data=dietMortW2,select=T, method="REML")
fitMort2<-gam(MMort~s(Complex, bs="re")+s(Pollock)+s(Year), data=dietMortW2, select=T, method="REML")
fitMort3<-gam(MMort~s(Complex, bs="re")+s(Pollock)+s(Other), data=dietMortW2, select=T, method="REML")
fitMort4<-gam(MMort~s(Complex, bs="re")+s(Other)+s(Year), data=dietMortW2, select=T, method="REML")

concurvity(fitMort1)
concurvity(fitMort2)
concurvity(fitMort3)
concurvity(fitMort4)

draw(compare_smooths(fitMort1,fitMort2,fitMort3,fitMort4))

# With year
pupMortResid1<-simulateResiduals(fitMort1, plot=F)
plot(pupMortResid1)
plotResiduals(pupMortResid1, form=dietMortW2$Complex)
plotResiduals(pupMortResid1, form=dietMortW2$Year)
plotResiduals(pupMortResid1, form=dietMortW2$Pollock)

# Without year
pupMortResid2<-simulateResiduals(fitMort3, plot=F)
plot(pupMortResid2)
plotResiduals(pupMortResid2, form=dietMortW2$Complex)
plotResiduals(pupMortResid2, form=dietMortW2$Year)
plotResiduals(pupMortResid2, form=dietMortW2$Pollock)
plotResiduals(pupMortResid2, form=dietMortW2$Year)

# Save final model outputs ------------------------------------------------

saveRDS(fitMort1, file.path(output.path,"Pup mortality gam output.rds"))
saveRDS(fitWeightsd1, file.path(output.path,"Pup weight sd gam output.rds"))
saveRDS(fitWeight1, file.path(output.path,"Pup weight gam output.rds"))

# Save final data frames
saveRDS(dietMortW, file.path(output.path,"Pup mortality data.rds"))
saveRDS(dietWeightW, file.path(output.path,"Pup weight data.rds"))

