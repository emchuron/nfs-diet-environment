# Runs gams to look at the influence of diet on August pup weights and early pup mortality

library(data.table)
library(tidyverse)
library(here)

output.path<-here("output")
input.path<-here("data")

# Read in response and predictor variables ---------------------------------

# Diet data
diet<-fread(file.path(input.path, "Diet response variables for analysis_Focal prey species.csv"))|>
  mutate(FO=freq/N)

# Pup metrics
pupMort<-readRDS(file.path(input.path, "Pup mortality by rookery.rds"))
pupWeight<-readRDS(file.path(input.path, "Pup weights by rookery.rds"))

# Join with diet data -----------------------------------------------------

dietVars<- diet |>
  mutate(Complex=as.factor(Complex), KLPreyGroup=as.factor(KLPreyGroup), Island=as.factor(Island)) |>
  group_by(Complex, KLPreyGroup)|>
  dplyr::mutate(mFO=mean(FO), anamFO=FO-mFO, sdFO=sd(FO))|>
  ungroup()

# Pup weight - has multiple joins because of pup sex
dietWeight<-dietVars |>
  left_join(pupWeight) |>
  dplyr::select(Complex, rookery,Year, Month, FO,anamFO,sdFO,KLPreyGroup, Sex, mWeight,sdWeight,nPup) |>
  filter(!is.na(mWeight)) 

# Prey groups of Pollock and not pollock - Sum
dietWeightW<-dietWeight |>
  mutate(KLPreyGroup2=ifelse(KLPreyGroup=="Pollock", "Pollock", "Other"))|>
  group_by(Complex, rookery,Year, Sex,KLPreyGroup2,mWeight,sdWeight,nPup)|>
  summarise(sumFO=sum(FO))|>
  pivot_wider(names_from=KLPreyGroup2, values_from=sumFO,id_cols = c(Complex,rookery,Year,Sex,mWeight,sdWeight,nPup))|>
  mutate(Complex=as.factor(Complex), rookery=factor(rookery), Exclude=case_when(Year==1994 & (Complex=="SGNorth" | Complex=="SGSouth")~"Yes",.default="No"))

#Mortality
dietMort<-dietVars |>
  left_join(pupMort)|>
  mutate(Mort=deadPups/pupsBorn)|>
  dplyr::select(Complex, rcod,Year,FO,KLPreyGroup,Mort,sdFO) |>
  filter(!is.na(Mort))

dietMortW<-dietMort |>
  mutate(KLPreyGroup2=ifelse(KLPreyGroup=="Pollock", "Pollock", "Other"))|>
  group_by(Complex,rcod, Year,Mort, KLPreyGroup2)|>
  summarise(sumFO=sum(FO))|>
  pivot_wider(names_from=KLPreyGroup2, values_from=sumFO,id_cols = c(Complex,rcod,Year,Mort))|>
  mutate(Complex=as.factor(Complex), rookery=factor(rcod))

# Summary ----------------------------------------------------------------
dietWeightSum<-dietWeightW|>
  filter(Exclude=="No")|>
  group_by(Year,Complex,Sex)|>
  summarise(TPup=sum(nPup))

dietMortW|>
  ungroup()|>
  summarise(MMort=mean(Mort))

summary(dietMortW)

# Analysis ----------------------------------------------------------------

library(mgcv)
library(gratia)
library(DHARMa)

# Notes
# The random effect of rookery was included in R1 instead of pre-analysis
# averaging. Now there are two random effects (Complex and rookery). I could 
# have included a third random effect s(Complex, rookery, bs="re"), which 
# estimates the deviations from the rookery and Complex averages but this 
# didn't seem necessary and smooths were nearly identical between the two.
# This addition reflects the patterns you see in the plots, where within a complex
# rookeries are not that variable in weight (and there is more consistent sampling)
# of the same rookeries, so there was really no effect of adding it here. It
# had a slight change to the mortality conclusions, as a plot of the data show
# a lot more variability within a complex, but the general pattern of
# the results was still the same. 

# Pup weights - mean --------------------------------------------------------
# Note - fitting the final model to individual pup sexes produces generally the same patterns
fitWeight1<-gam(mWeight~Sex+s(Year, by=Sex)+s(Complex, bs="re")
                + s(rookery,bs="re")+s(Other, by=Sex)+s(Pollock, by=Sex), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")
fitWeight1a<-gam(mWeight~Sex+s(Complex, bs="re")+s(rookery, bs="re")+s(Year,by=Sex)+s(Other,by=Sex)+s(Pollock,by=Sex), data=subset(dietWeightW), select=T, method="REML") # all the data - no exclusions
fitWeight2<-gam(mWeight~Sex+s(Complex, bs="re")+s(Pollock, by=Sex)+s(rookery, bs="re") + s(Year, by=Sex), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")
fitWeight3<-gam(mWeight~Sex+s(Complex, bs="re")+s(Year, by=Sex) +s(rookery, bs="re")+ s(Other, by=Sex), data=subset(dietWeightW, Exclude=="No"),select=T, method="REML")
fitWeight4<-gam(mWeight~Sex+s(Complex, bs="re")+s(Pollock, by=Sex) + s(Other, by=Sex)+s(rookery, bs="re"), data=subset(dietWeightW, Exclude=="No"),select=T, method="REML")

gam.check(fitWeight1)
gam.check(fitWeight2)
gam.check(fitWeight3)
gam.check(fitWeight4)

appraise(fitWeight1)
appraise(fitWeight2)
appraise(fitWeight3)
appraise(fitWeight4)

concurvity(fitWeight1)
concurvity(fitWeight2)
concurvity(fitWeight3)
concurvity(fitWeight4)

# To see how smooths change
draw(compare_smooths(fitWeight1,fitWeight2,fitWeight3,fitWeight4))

# With year
pupWeightResid1<-simulateResiduals(fitWeight1, plot=F)
plot(pupWeightResid1)
plotResiduals(pupWeightResid1, form =dietWeightW$Year[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW$Pollock[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW$Other[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW$Complex[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid1, form =dietWeightW$rookery[dietWeightW$Exclude=="No"])

# Without year
pupWeightResid2<-simulateResiduals(fitWeight4, plot=F)
plotResiduals(pupWeightResid2, form =dietWeightW$Year[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW$Pollock[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW$Other[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW$Complex[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightResid2, form =dietWeightW$rookery[dietWeightW$Exclude=="No"])

# Pup weights sd --------------------------------------------------------

fitWeightsd1<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Pollock, by=Sex)+s(Year, by=Sex)+s(Other, by=Sex)+s(rookery, bs="re"), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")
fitWeightsd1a<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Pollock, by=Sex)+s(Year, by=Sex)+s(Other, by=Sex)+s(rookery, bs="re"), data=subset(dietWeightW), select=T, method="REML")
fitWeightsd2<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Pollock, by=Sex)+s(Other, by=Sex)+s(rookery, bs="re"), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")
fitWeightsd3<-gam(sdWeight~Sex+s(Complex, bs="re")+s(Year, by=Sex)+s(Other, by=Sex)+s(rookery, bs="re"), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")
fitWeightsd4<-gam(sdWeight~Sex+s(Complex, bs="re")+s(rookery, bs="re")+s(Year, by=Sex)+s(Pollock, by=Sex), data=subset(dietWeightW, Exclude=="No"), select=T, method="REML")

gam.check(fitWeightsd1)
gam.check(fitWeightsd2)
gam.check(fitWeightsd3)
gam.check(fitWeightsd4)

concurvity(fitWeightsd1)
concurvity(fitWeightsd2)
concurvity(fitWeightsd3)
concurvity(fitWeightsd4)

appraise(fitWeightsd1)
appraise(fitWeightsd2)
appraise(fitWeightsd3)
appraise(fitWeightsd4)

draw(compare_smooths(fitWeightsd1,fitWeightsd2,fitWeightsd3,fitWeightsd4))

pupWeightSDResid1<-simulateResiduals(fitWeightsd1, plot=F)
pupWeightSDResid2<-simulateResiduals(fitWeightsd2, plot=F)

plot(pupWeightSDResid1)
plotResiduals(pupWeightSDResid1, form =dietWeightW$Year[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW$Other[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW$Pollock[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW$Complex[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid1, form =dietWeightW$rookery[dietWeightW$Exclude=="No"])

plot(pupWeightSDResid2)
plotResiduals(pupWeightSDResid2, form =dietWeightW$Year[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW$Other[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW$Pollock[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW$Complex[dietWeightW$Exclude=="No"])
plotResiduals(pupWeightSDResid2, form =dietWeightW$rookery[dietWeightW$Exclude=="No"])

# Pup Mortality --------------------------------------------------------
fitMort1<-gam(Mort~s(Complex, bs="re")+s(Pollock)+s(rookery, bs="re")+s(Other)+s(Year), data=dietMortW,select=T, method="REML")
fitMort2<-gam(Mort~s(Complex, bs="re")+s(Pollock)+s(Year)+s(rookery, bs="re"), data=dietMortW, select=T, method="REML")
fitMort3<-gam(Mort~s(Complex, bs="re")+s(Pollock)+s(Other)+s(rookery, bs="re"), data=dietMortW, select=T, method="REML")
fitMort4<-gam(Mort~s(Complex, bs="re")+s(Other)+s(Year)+s(rookery, bs="re"), data=dietMortW, select=T, method="REML")

concurvity(fitMort1)
concurvity(fitMort2)
concurvity(fitMort3)
concurvity(fitMort4)

draw(compare_smooths(fitMort1,fitMort2,fitMort3,fitMort4))

# With year
pupMortResid1<-simulateResiduals(fitMort1, plot=F)
plot(pupMortResid1)
plotResiduals(pupMortResid1, form=dietMortW$Complex)
plotResiduals(pupMortResid1, form=dietMortW$Year)
plotResiduals(pupMortResid1, form=dietMortW$Pollock)
plotResiduals(pupMortResid1, form=dietMortW$rookery)

# Without year
pupMortResid2<-simulateResiduals(fitMort3, plot=F)
plot(pupMortResid2)
plotResiduals(pupMortResid2, form=dietMortW$Complex)
plotResiduals(pupMortResid2, form=dietMortW$Year)
plotResiduals(pupMortResid2, form=dietMortW$Pollock)
plotResiduals(pupMortResid2, form=dietMortW$rookery)

# Save final model outputs ------------------------------------------------

saveRDS(fitMort1, file.path(output.path,"Pup mortality gam output.rds"))
saveRDS(fitWeightsd1, file.path(output.path,"Pup weight sd gam output.rds"))
saveRDS(fitWeight1, file.path(output.path,"Pup weight gam output.rds"))

# Save final data frames
saveRDS(dietMortW, file.path(output.path,"Pup mortality data with diet.rds"))
saveRDS(dietWeightW, file.path(output.path,"Pup weight data with diet.rds"))

