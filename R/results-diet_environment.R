
library(data.table)
library(tidyverse)
library(gratia)
library(patchwork)
library(gam.hp)

# Paths
fig.path<-file.path(here::here(), "figures")
input.path<-file.path(here::here(), "data")
output.path<-file.path(here::here(), "output")

# Functions
source(file.path(here::here(), 'R',"gam_hp_weight.R"))

dev_exp<-dev_expl <- function(model){
  if(is.null(model$null.deviance) | is.null(model$deviance)){
    stop("model$null.deviance or model$deviance is NULL.")
  }
  return((model$null.deviance - model$deviance)/model$null.deviance * 100)
}

# GAM outputs
fitSpecies<-readRDS(file.path(output.path,"Diet-environmental gam output by species.rds"))

# Diet data
diet<-fread(file.path(input.path, "Diet response variables for analysis_Focal prey species.csv"))|>
  mutate(FO=freq/N)

# Model summaries ---------------------------------------------------------
linkTransform<-inv_link(fitSpecies$fitG[[1]])

fitSum<-fitSpecies |>
  mutate(tidied=map(fitG, broom::tidy),
         deviance=map(fitG, dev_exp)) |>
  unnest(c(tidied,deviance)) |>
  dplyr::select(-c(fitNullNoYear,fitNull,fitG))|>
  mutate(p.value=ifelse(p.value>0.05, NA, p.value)) |>
  filter(!is.na(p.value)) |>
  arrange(p.value)  |>
  filter(p.value<0.04 & !(KLPreyGroup=="Sandlance" & term=="s(bott_s)")& !(KLPreyGroup=="GBBM" & term=="s(sstReShelf)"))

table(fitSum$term)

fitSum |>
  group_by(KLPreyGroup) |>
  summarise(dev=deviance[1])

# Relative contribution of variables to deviance explained
fitSumIndDev<-fitSpecies |>
  mutate(deviance=map(fitG, ~gam_hp_weight(.x) |> pluck(2)),
         devianceCommon=map(fitG, ~gam_hp_weight(.x, commonality=T) |> pluck(2)))|>
  dplyr::select(KLPreyGroup, data, deviance, devianceCommon)

# Diet summaries ---------------------------------------------------------

dietSum<-diet |>
  group_by(Complex, KLPreyGroup)|>
  dplyr::summarise(mFO=round(mean(FO), digits=3),sdFO=sd(FO), minFO=round(min(FO),digits=3), maxFO=round(max(FO), digits=3),
                   n=mean(N), sdn=sd(N),nSum=sum(N), minN=min(N), maxN=max(N),
                   nYear=length(unique(Year))) |>
  ungroup()


# Environmental effects with partial residuals ----------------------------

complexPal<-c("#00346E", "#007CBF", "#06ABDF", "#EDD03E", "#F5A200", "#6D8600", "#424D0C") # from(palettes_d$feathers$bee_eater)

fitSpecies<-fitSpecies |>
  mutate(pred=map(fitG,~smooth_estimates(.x)),
         pred2=map(fitG,~smooth_estimates(.x) |> pivot_longer(bott_s:walleye_pollock_mean_cpue_no) |> filter(!is.na(value))),
         conf=map(pred, ~add_confint(.x) |> pivot_longer(bott_s:walleye_pollock_mean_cpue_no) |> filter(!is.na(value))),
         resid=map2(data,fitG, ~add_partial_residuals(.x,.y)),
         resid2=map2(resid,pred, ~dplyr::select(.x,Year,Complex,Island,which(colnames(.x) %in% colnames(.y)), "s(walleye_pollock_mean_cpue_no)","s(bott_s)","s(sstReShelf)","s(Year)")),
         resid3=map(resid2, ~pivot_longer(.x,cols=sstReShelf:bott_s, names_to="name",values_to = "data_value")),
         resid4=map(resid3, ~ pivot_longer(.x,cols=c("s(sstReShelf)","s(walleye_pollock_mean_cpue_no)","s(bott_s)","s(Year)"), names_to="smooth", values_to="smooth_value")|>
                      mutate(smooth=substring(smooth,3, nchar(smooth)-1))|>
                      filter(name==smooth)))

partialSpecies<- fitSpecies |> 
  dplyr::select(KLPreyGroup,conf)|>
  unnest(conf)

partialResidSpecies<- fitSpecies |> 
  dplyr::select(KLPreyGroup,resid4)|>
  unnest(resid4) |>
  mutate(smooth=paste0("s(",smooth, ")"))

# Bottom temperature
GBottPR<-ggplot(subset(partialSpecies, name=="bott_s" & paste0(KLPreyGroup,.smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=value, y=.estimate))+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, x=value),alpha = 0.2, fill="black")+
  geom_point(data=subset(partialResidSpecies, name=="bott_s" & paste0(KLPreyGroup,smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=data_value, y=smooth_value, color=Complex,fill=Complex),pch=21)+
  scale_color_manual(values=colorspace::darken(complexPal,0.2), name="",labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  scale_fill_manual(values=complexPal, name="", labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  facet_wrap(~KLPreyGroup, nrow=1, scales="free_y")+
  ggthemes::theme_few()+
  labs(y = "Partial effect", x=bquote("Mean bottom temperature ("*degree*"C)"))+
  theme(legend.position="bottom", legend.text=element_text(size=12), strip.text=element_text(size=12),
        axis.text=element_text(size=12))+
  scale_x_continuous(expand=c(0.02,0))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# Pollock
GPollPR<-ggplot(subset(partialSpecies, name=="walleye_pollock_mean_cpue_no" & paste0(KLPreyGroup,.smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=value, y=.estimate))+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, x=value),alpha = 0.2, fill="black")+
  geom_point(data=subset(partialResidSpecies, name=="walleye_pollock_mean_cpue_no" & paste0(KLPreyGroup,smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=data_value, y=smooth_value, color=Complex,fill=Complex), pch=21)+
  scale_color_manual(values=colorspace::darken(complexPal,0.2), name="",labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  scale_fill_manual(values=complexPal, name="", labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  facet_wrap(~KLPreyGroup, nrow=1, scales="free_y")+
  ggthemes::theme_few()+
  labs(y = "Partial effect", x=bquote("Mean pollock abundance (no"~km^-2*")"))+
  theme(legend.position="bottom", legend.text=element_text(size=12), strip.text=element_text(size=12),
        axis.text=element_text(size=12))+
  scale_x_continuous(expand=c(0.02,0))+
  guides(color = guide_legend(override.aes = list(size = 5)))

# SST
GSurfPR<-ggplot(subset(partialSpecies, name=="sstReShelf" & paste0(KLPreyGroup,.smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=value, y=.estimate))+
  geom_line()+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, x=value),alpha = 0.2, fill="black")+
  geom_point(data=subset(partialResidSpecies, name=="sstReShelf" & paste0(KLPreyGroup,smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=data_value, y=smooth_value, color=Complex,fill=Complex),pch=21)+
  #geom_text(data=subset(partialResidSpecies, name=="sstReShelf" & paste0(KLPreyGroup,smooth) %in% paste0(fitSum$KLPreyGroup, fitSum$term)), aes(x=data_value, y=smooth_value, color=Complex,fill=Complex,label=Year))+
  scale_color_manual(values=colorspace::darken(complexPal,0.2), name="",labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  scale_fill_manual(values=complexPal, name="", labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South"))+
  facet_wrap(~KLPreyGroup, nrow=1, scales="free_y")+
  ggthemes::theme_few()+
  theme(legend.position="bottom", legend.text=element_text(size=12), strip.text=element_text(size=12),
        axis.text=element_text(size=12))+
  scale_x_continuous(expand=c(0.02,0))+
  labs(y = "Partial effect", x=bquote("Mean surface temperature ("*degree*"C)"))+
  guides(color = guide_legend(override.aes = list(size = 5)))
  
layout <- "
AAAAA
BBBB#
CCCC#
"
GAllR<-GBottPR/GSurfPR/GPollPR +
  plot_layout(design = layout) &
  plot_annotation(tag_levels="a") & 
  theme(plot.tag = element_text(size = 12,face="bold"))


