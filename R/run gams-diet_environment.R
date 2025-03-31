

library(mgcv)
library(gratia)
library(data.table)
library(tidyverse)

output.path<-file.path(here::here(), "output")
input.path<-file.path(here::here(), "data")

# Read in response and predictor variables --------------------------------

# Diet data
diet<-fread(file.path(input.path, "Diet response variables for analysis_Focal prey species.csv"))|>
  mutate(FO=freq/N)

# Predictor variables
surveyVars<-fread(file.path(input.path, "Independent predictor variables - Survey.csv"))
remoteSST<-fread(file.path(input.path, "Independent predictor variables - Remote SST.csv"))

# Rearrange and combine
surveyVarsW<-surveyVars |>
  pivot_wider(values_from=value,names_from=var) %>%
  janitor::clean_names()|>
  setnames(c("complex","year"), c("Complex","Year")) 

remoteSSTW<-remoteSST |>
  pivot_wider(values_from=value,names_from=var)

allVars<-list(remoteSSTW, surveyVarsW)|>
 reduce(left_join)

# Create linked datasets for each complex ---------------------------------

dietVars<- dplyr::left_join(diet,allVars) |>
  mutate(Complex=as.factor(Complex), KLPreyGroup=as.factor(KLPreyGroup), Island=as.factor(Island)) |>
  tidytable::group_split(Complex, .named=T)|>
  map(droplevels)

dietVarsUn<-rbindlist(dietVars)

dietSum<-dietVarsUn |>
  group_by(Complex, Year) |>
  dplyr::summarise(n=unique(N))

# Correlations among predictor variables ----------------------------------

# At the complex level
corEast<-cor(dietVars$SPEast[dietVars$SPEast$KLPreyGroup=="Pollock",c(9:36)], method="pearson")
corEBay<-cor(dietVars$SPEnglishBay[dietVars$SPEnglishBay$KLPreyGroup=="Pollock",c(9:36)], method="pearson")
corReef<-cor(dietVars$SPReefPoint[dietVars$SPReefPoint$KLPreyGroup=="Pollock",c(9:36)], method="pearson")
corNorth<-cor(dietVars$SGNorth[dietVars$SGNorth$KLPreyGroup=="Pollock",c(9:36)], method="pearson")
corSouth<-cor(dietVars$SGSouth[dietVars$SGSouth$KLPreyGroup=="Pollock",c(9:36)], method="pearson")

# At the broader level
corAll<-cor(dietVarsUn[dietVarsUn$KLPreyGroup=="Pollock",c(9:36)], method="pearson")


# Fit the analysis --------------------------------------------------------
library(DHARMa)

fitSpecies<-dietVarsUn |>
  mutate(Year2=as.factor(Year))|>
  group_by(KLPreyGroup)|>
  mutate(.row = row_number())|>
  nest() |>
  mutate(fitNullNoYear=purrr::map(data, ~  gam(FO ~ s(Complex, bs="re"),
                                               data = .x, binomial("logit"), method="REML", select=T, weights=N)),
    fitNull= purrr::map(data, ~  gam(FO ~ s(Complex, bs="re")+s(Year),
                                          data = .x, binomial("logit"), method="REML", select=T, weights=N)),
    fitG = purrr::map(data, ~  gam(FO ~ s(Complex, bs="re")+s(bott_s)+s(sstReShelf)+s(walleye_pollock_mean_cpue_no) +s(Year), 
                                   data = .x, binomial("logit"), method="REML", select=T, weights=N)),
        fitGNoYear = purrr::map(data, ~  gam(FO ~ s(Complex, bs="re")+s(bott_s)+s(sstReShelf)+s(walleye_pollock_mean_cpue_no),
                                    data = .x, binomial("logit"), method="REML", select=T, weights=N)))

  
# Look at assumptions -----------------------------------------------------

fitAssump<-fitSpecies |>
  mutate(partCorr=map(fitG, ~pacf(resid(.x))),
         autoCorr=map(fitG, ~acf(resid(.x))),
         simResid=map(fitG, simulateResiduals)) |>
  dplyr::select(-c(fitNull))

# Residuals plots from dHARMA 
plot(fitAssump$simResid[[1]])
testDispersion(fitAssump$simResid[[1]])
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$bott_s)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$sstReShelf)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$Complex)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$Year)

plot(fitAssump$simResid[[2]])
testDispersion(fitAssump$simResid[[2]])
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$bott_s)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$sstReShelf)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$Complex)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$Year)

plot(fitAssump$simResid[[3]])
testOverdispersion(fitAssump$simResid[[3]])
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$bott_s)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$sstReShelf)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$Complex)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$Year)

plot(fitAssump$simResid[[4]])
testDispersion(fitAssump$simResid[[4]])
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$bott_s)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$sstReShelf)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$Complex)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$Year)

plot(fitAssump$simResid[[5]])
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$bott_s)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$sstReShelf)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$Complex)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$Year)

plot(fitAssump$simResid[[6]])
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$bott_s)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$sstReShelf)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$Complex)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$Year)

plot(fitAssump$simResid[[7]])
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$bott_s)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$sstReShelf)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$Complex)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$Year)

plot(fitAssump$simResid[[8]])
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$bott_s)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$sstReShelf)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$Complex)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$Year)

# Autocorrelation plots
fitAssump$autoCorr
fitAssump$partCorr

# Look at assumptions - No Year-----------------------------------------------------

fitAssump<-fitSpecies |>
  mutate(partCorr=map(fitGNoYear, ~pacf(resid(.x))),
         autoCorr=map(fitGNoYear, ~acf(resid(.x))),
         simResid=map(fitGNoYear, simulateResiduals)) |>
  dplyr::select(-c(fitNull))

# Residuals plots from dHARMA 
plot(fitAssump$simResid[[1]])
testOverdispersion(fitAssump$simResid[[1]])
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$bott_s)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$sstReShelf)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$Complex)
plotResiduals(fitAssump$simResid[[1]], form =fitSpecies$data[[1]]$Year)

plot(fitAssump$simResid[[2]])
testOverdispersion(fitAssump$simResid[[2]])
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$bott_s)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$sstReShelf)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$Complex)
plotResiduals(fitAssump$simResid[[2]], form =fitSpecies$data[[2]]$Year)

plot(fitAssump$simResid[[3]])
testOverdispersion(fitAssump$simResid[[3]])
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$bott_s)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$sstReShelf)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$Complex)
plotResiduals(fitAssump$simResid[[3]], form =fitSpecies$data[[3]]$Year)

plot(fitAssump$simResid[[4]])
testOverdispersion(fitAssump$simResid[[4]])
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$bott_s)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$sstReShelf)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$Complex)
plotResiduals(fitAssump$simResid[[4]], form =fitSpecies$data[[4]]$Year)

plot(fitAssump$simResid[[5]])
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$bott_s)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$sstReShelf)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$Complex)
plotResiduals(fitAssump$simResid[[5]], form =fitSpecies$data[[5]]$Year)

plot(fitAssump$simResid[[6]])
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$bott_s)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$sstReShelf)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$Complex)
plotResiduals(fitAssump$simResid[[6]], form =fitSpecies$data[[6]]$Year)

plot(fitAssump$simResid[[7]])
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$bott_s)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$sstReShelf)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$Complex)
plotResiduals(fitAssump$simResid[[7]], form =fitSpecies$data[[7]]$Year)

plot(fitAssump$simResid[[8]])
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$bott_s)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$sstReShelf)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$walleye_pollock_mean_cpue_no)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$Complex)
plotResiduals(fitAssump$simResid[[8]], form =fitSpecies$data[[8]]$Year)

# Comparing year vs no year smooths ---------------------------------------

draw(compare_smooths(fitSpecies$fitG[[1]], fitSpecies$fitGNoYear[[1]]))  # No real diff
draw(compare_smooths(fitSpecies$fitG[[2]], fitSpecies$fitGNoYear[[2]])) # No big differences
draw(compare_smooths(fitSpecies$fitG[[3]], fitSpecies$fitGNoYear[[3]])) # No big differences
draw(compare_smooths(fitSpecies$fitG[[4]], fitSpecies$fitGNoYear[[4]])) # No big differences
draw(compare_smooths(fitSpecies$fitG[[5]], fitSpecies$fitGNoYear[[5]])) # No real differences
draw(compare_smooths(fitSpecies$fitG[[6]], fitSpecies$fitGNoYear[[6]])) # No real differences
draw(compare_smooths(fitSpecies$fitG[[7]], fitSpecies$fitGNoYear[[7]])) # A few differences
draw(compare_smooths(fitSpecies$fitG[[8]], fitSpecies$fitGNoYear[[8]])) # No difference

# Autocorrelation plots
fitAssump$autoCorr
fitAssump$partCorr


# Alternate model set up --------------------------------------------------

#https://stats.stackexchange.com/questions/487135/how-to-use-a-gam-to-predict-the-probability-in-binomial-data-as-a-function-of-pr

# This gives the same results as above. Just doing it to verify that how I set up the analysis doesn't affect the results
fitSpeciesRaw<-dietVarsUn |>
  group_by(KLPreyGroup)|>
  mutate(.row = row_number())|>
  nest() |>
  mutate(fitG = purrr::map(data, ~  gam(cbind(freq,N-freq) ~ s(Complex, bs="re")+s(bott_s)+s(sstReShelf)+s(walleye_pollock_mean_cpue_no) +s(Year), 
                                        data = .x, family="binomial", method="REML", select=T)))

draw(compare_smooths(fitSpecies$fitG[[1]], fitSpeciesRaw$fitG[[1]])) 


# Save output -------------------------------------------------------------

saveRDS(fitSpecies, file.path(output.path,"Diet-environmental gam output by species.rds"))

