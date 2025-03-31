# Results and figures for the analysis between pup metrics and diet

library(tidyverse)
library(patchwork)
library(gratia)
library(gam.hp)

output.path<-file.path(here::here(), "output")
figure.path<-file.path(here::here(), "figures")

# Load in data outputs -----------------------------------------------------

# gam outputs
fitmort<-readRDS(file.path(output.path,"Pup mortality gam output.rds"))
fitweightsd<-readRDS(file.path(output.path,"Pup weight sd gam output.rds"))
fitweight<-readRDS(file.path(output.path,"Pup weight gam output.rds"))

# data
mort<-readRDS(file.path(output.path,"Pup mortality data.rds")) |>
  mutate(Other=GBBM+GMGM+Herring+HexSable+Salmon+Sandlance+Smoothtongue) 

weight<-readRDS(file.path(output.path,"Pup weight data.rds")) |>
  mutate(Exclude=case_when(Year==1994 & (Complex=="SGNorth" | Complex=="SGSouth")~"Yes",.default="No"),
         Other=GBBM+GMGM+Herring+HexSable+Salmon+Sandlance+Smoothtongue) |> 
  filter(Exclude=="No")

# Contribution of variables -----------------------------------------------
gam.hp(fitmort)
gam.hp(fitweightsd)
gam.hp(fitweight)

# Summaries ---------------------------------------------------------------

weightSum<-weight |>
  group_by(Complex, Sex) |>
  summarise(weight=mean(mWeight), sdweight=sd(mWeight),minweight=min(mWeight),maxweight=max(mWeight),
            totN=sum(nPup), nYear=length(unique(Year)))

mortSum<-mort |>
  group_by(Complex)|>
  summarise(Mort=mean(MMort), sdMort=sd(MMort), minMort=min(MMort), maxMort=max(MMort),
            n=length(MMort))
# Plots -------------------------------------------------------------------

mortSmooths<-smooth_estimates(fitmort) |> add_confint() 
weightSmooths<-smooth_estimates(fitweight) |> add_confint()
weightsdSmooths<-smooth_estimates(fitweightsd) |> add_confint()
mort<-add_partial_residuals(mort, fitmort)
weight2<-add_partial_residuals(weight, fitweight)

GMortOther<-ggplot(subset(mortSmooths, !is.na(Other)), aes(x=Other, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2)+
  geom_rug(data=mort,sides="b", inherit.aes=F, aes(x=Other,y=0),linewidth=0.25)+
  #geom_text(data=mort, aes(label=Year, color=Complex,y=`s(Other)`))+
  ggthemes::theme_few()+
  ylab("Partial effect - mortality")+
  xlab("Other prey FO")

GMortYear<-ggplot(subset(mortSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2)+
  geom_rug(data=mort,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  ggthemes::theme_few()+
  ylab("Partial effect - mortality")+
  xlab("Year")+
  theme(plot.title=element_text(size=12, hjust=0.5))+
  scale_x_continuous(expand=c(0.01,0))

GWeightsOther<-ggplot(subset(weightSmooths,!is.na(Other)), aes(x=Other, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2)+
  geom_rug(data=weight2,sides="b", inherit.aes=F, aes(x=Other,y=0),linewidth=0.25)+
  #geom_text(data=weight2, aes(label=Year, color=Complex,y=`s(Other)`))+
  ggthemes::theme_few()+
  ylab("Partial effect - mass")+
  xlab("Other prey FO")+
  scale_x_continuous(expand=c(0.01,0))

GWeightsYear<-ggplot(subset(weightSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2,)+
  geom_rug(data=weight,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  ggthemes::theme_few()+
  ylab("Partial effect - mass")+
  xlab("Year")+
  scale_x_continuous(expand=c(0.01,0))+
  theme(plot.title=element_text(size=12, hjust=0.5))

GWeightsSDYear<-ggplot(subset(weightsdSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2)+
  geom_rug(data=weight,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  ggthemes::theme_few()+
  ylab("Partial effect - mass sd")+
  xlab("Year")+
  scale_x_continuous(expand=c(0.01,0))+
  theme(plot.title=element_text(size=12, hjust=0.5))

GAllYear<-(GMortYear/GWeightsYear/GWeightsSDYear)+plot_layout(axis_title="collect")+
  plot_annotation(title="a")& 
  theme(plot.title = element_text(size = 12,face="bold"))

GAllDiet<-GMortOther/GWeightsOther/plot_spacer()+
  plot_layout(axis_title="collect")+ plot_annotation(title="b")& 
  theme(plot.title = element_text(size = 12,face="bold"))

GAll<-wrap_elements(GAllYear)|wrap_elements(GAllDiet)

#ggsave(file.path(figure.path, "Figure_5.tiff"), GAll,width=7, height=8, compression="lzw")
ggsave(file.path(figure.path, "Figure_5.pdf"), GAll,width=7, height=8,device=cairo_pdf,family="Arial",dpi=600)

