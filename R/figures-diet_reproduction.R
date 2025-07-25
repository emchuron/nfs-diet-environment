# Results and figures for the analysis between pup metrics and diet

library(tidyverse)
library(here)
library(patchwork)
library(gratia)
library(gam.hp)

output.path<-here("output")
figure.path<-here("figures")

# Load in data outputs -----------------------------------------------------

# gam outputs
fitmort<-readRDS(file.path(output.path,"Pup mortality gam output.rds"))
fitweightsd<-readRDS(file.path(output.path,"Pup weight sd gam output.rds"))
fitweight<-readRDS(file.path(output.path,"Pup weight gam output.rds"))

# data
mort<-readRDS(file.path(output.path,"Pup mortality data with diet.rds")) 

weight<-readRDS(file.path(output.path,"Pup weight data with diet.rds")) |>
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
  summarise(MMort=mean(Mort), sdMort=sd(Mort), minMort=min(Mort), maxMort=max(Mort),
            n=length(unique(Year)))
# Plots -------------------------------------------------------------------

mortSmooths<-smooth_estimates(fitmort) |> add_confint() 
weightSmooths<-smooth_estimates(fitweight) |> add_confint()
weightsdSmooths<-smooth_estimates(fitweightsd) |> add_confint()
mort<-add_partial_residuals(mort, fitmort)
weight2<-add_partial_residuals(weight, fitweight)

GMortYear<-ggplot(subset(mortSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci), alpha=0.2)+
  geom_rug(data=mort,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  ggthemes::theme_few()+
  ylab("Partial effect - mortality")+
  xlab("Year")+
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text=element_text(size=10))+
  scale_x_continuous(expand=c(0.01,0))

GWeightsOther<-ggplot(subset(weightSmooths,!is.na(Other)), aes(x=Other, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line(aes(color=Sex))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci,fill=Sex), alpha=0.2)+
  geom_rug(data=weight2,sides="b", inherit.aes=F, aes(x=Other,y=0),linewidth=0.25)+
  scale_color_manual(values=c("red","blue"), labels=c("Female","Male"))+
  scale_fill_manual(values=c("red","blue"),labels=c("Female","Male"))+
  #geom_text(data=weight2, aes(label=Year, color=Complex,y=`s(Other)`))+
  ggthemes::theme_few()+
  ylab("Partial effect - mass")+
  xlab("Other prey FO")+
  scale_x_continuous(expand=c(0.01,0))+
  theme(legend.position="none")

GWeightsYear<-ggplot(subset(weightSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line(aes(color=Sex))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=Sex), alpha=0.2,)+
  geom_rug(data=weight,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  scale_color_manual(values=c("red","blue"), labels=c("Female","Male"))+
  scale_fill_manual(values=c("red","blue"),labels=c("Female","Male"))+
  ggthemes::theme_few()+
  ylab("Partial effect - mass")+
  xlab("Year")+
  scale_x_continuous(expand=c(0.01,0))+
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text=element_text(size=10),
        legend.text=element_text(size=10), legend.title=element_blank(),legend.position="inside",legend.position.inside=c(0.2,0.2))

GWeightsSDYear<-ggplot(subset(weightsdSmooths,!is.na(Year)), aes(x=Year, y=.estimate))+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  geom_line(aes(color=Sex))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=Sex), alpha=0.2)+
  geom_rug(data=weight,sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  scale_color_manual(values=c("red","blue"), labels=c("Female","Male"))+
  scale_fill_manual(values=c("red","blue"),labels=c("Female","Male"))+
  ggthemes::theme_few()+
  ylab("Partial effect - mass sd")+
  xlab("Year")+
  scale_x_continuous(expand=c(0.01,0))+
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text=element_text(size=10),
        legend.position="none")

GWeight<-GWeightsYear+GWeightsOther+plot_layout(axes="collect")

GAll<-GWeight/(GWeightsSDYear+GMortYear)+plot_annotation(tag_levels=list(c("a","","b","c")))& 
  theme(plot.tag = element_text(size = 10,face="bold"))

ggsave(file.path(figure.path, "Figure_5.pdf"), GAll,width=7, height=6,device=cairo_pdf,family="Arial",dpi=600)

