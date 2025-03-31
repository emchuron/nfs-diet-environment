# Create the figures associated with the diet - environmental analysis

library(marginaleffects)
library(paletteer)
library(here)
library(tidyverse)
library(ggbreak)
library(ggokabeito)

fig.path<-file.path(here::here(), "figures")
input.path<-file.path(here::here(), "data")
output.path<-file.path(here::here(), "output")

dev_exp<-dev_expl <- function(model){
  if(is.null(model$null.deviance) | is.null(model$deviance)){
    stop("model$null.deviance or model$deviance is NULL.")
  }
  return((model$null.deviance - model$deviance)/model$null.deviance * 100)
}

# Read in environmental data ------------------------------------------------------------

# Predictor variables
surveyVars<-fread(file.path(input.path, "Independent predictor variables - Survey.csv"))
remoteSST<-fread(file.path(input.path, "Independent predictor variables - Remote SST.csv"))

surveyVarsW<-surveyVars |>
  pivot_wider(values_from=value,names_from=var) %>%
  janitor::clean_names()|>
  setnames(c("complex","year"), c("Complex","Year")) 

remoteSSTW<-remoteSST |>
  pivot_wider(values_from=value,names_from=var)

predVars<-list(filter(remoteSSTW, Month==8), surveyVarsW)|>
  reduce(left_join) |>
  dplyr::select(Year, Complex, sstReShelf,bott_s,walleye_pollock_mean_cpue_no) |>
  mutate(walleye_pollock_mean_cpue_no=walleye_pollock_mean_cpue_no/1000)|>
  pivot_longer(-c(Complex,Year)) |>
  mutate(Prediction=ifelse(Year>2012, "Yes", "No"),
         Complex=factor(Complex, levels=c("SPEast","SPEnglishBay","SPReefPoint","SGNorth","SGSouth"),
                        labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South")))

# Diet model output ------------------------------------------------------------

fitSpecies<-readRDS(file.path(output.path,"Diet-environmental gam output by species.rds"))
linkTransform<-inv_link(fitSpecies$fitG[[1]])

fitSpecies<-fitSpecies|>
  mutate(sumPoll=purrr::map(fitG,~plot_predictions(.x, condition=c("walleye_pollock_mean_cpue_no"), draw=FALSE,type="link",transform=linkTransform)),
         sumBott=purrr::map(fitG,~plot_predictions(.x, condition=c("bott_s"), draw=FALSE,type="link",transform=linkTransform)),
         sumSurf=purrr::map(fitG,~plot_predictions(.x, condition=c("sstReShelf"), draw=FALSE,type="link",transform=linkTransform)))

fo<-fitSpecies |>
  dplyr::select(KLPreyGroup, data)|>
  unnest(data) 

fitSum<-fitSpecies |>
  mutate(tidied=purrr::map(fitG, broom::tidy),
         deviance=purrr::map(fitG, dev_exp),
         devianceNull=purrr::map(fitNull, dev_exp),
         devianceNullNo=purrr::map(fitNullNoYear, dev_exp)) |>
  unnest(c(tidied,deviance, devianceNull, devianceNullNo)) |>
  dplyr::select(-c(fitNullNoYear,fitNull,fitG))|>
  mutate(p.value=ifelse(p.value>0.05, NA, p.value)) |>
  filter(!is.na(p.value)) |>
  arrange(p.value)  |>
  filter(p.value<0.04 & !(KLPreyGroup=="Sandlance" & term=="s(bott_s)")& !(KLPreyGroup=="GBBM" & term=="s(sstReShelf)"))


# General plot stuff ------------------------------------------------------

complexPal<-c("#00346E", "#007CBF", "#06ABDF", "#EDD03E", "#F5A200", "#6D8600", "#424D0C")
speciesCol<-grafify:::graf_palettes$fishy |>as.character()

# Figure 1 - Foraging ranges ---------------------------------------------------------

af95<-file.path(input.path, "Adult female 95% UD.rds")|>
  readRDS()|>
  st_transform(" +proj=longlat +datum=WGS84 +no_defs ") |>
  st_shift_longitude() |>
  mutate(Complex=factor(id, levels=c("SPEast","SPEnglishBay","SPReefPoint","SGNorth","SGSouth"),
                        labels=c("SP - East","SP - English Bay", "SP - Reef Point", "SG - North","SG - South")))

# Alaska map
alaska <-  rnaturalearth::ne_states(country="United States of America", returnclass="sf") |>
  filter(name=="Alaska")|>
  st_transform(" +proj=longlat +datum=WGS84 +no_defs ")|>
  st_shift_longitude()

# Bathymetry
bathy <- akgfmaps::get_base_layers(select.region = "bs.all", set.crs = st_crs(af95))$bathymetry|>
  st_shift_longitude()

pribs<-data.frame(lat=c(57.1208,56.6016), lon=c(189.717,190.4519),
                  island=c("St. Paul","St. George"))

GPribs<-ggplot(NULL)+
  ggthemes::theme_few()+
  geom_sf(data=bathy, color="black", linewidth=0.2)+
  geom_sf(data=af95, alpha=0.5,aes(fill=Complex, color=Complex))+
  geom_sf(data=alaska, fill="gray60",color="black")+
  coord_sf(ylim = c(52,61.5), xlim = c(180, 200))+
  ylab(NULL)+
  xlab(NULL)+
  ggthemes::theme_few()+
  theme(axis.text=element_text(size=12), legend.position="bottom")+
  scale_fill_manual(values=complexPal,name="")+
  scale_color_manual(values=complexPal,name="")+
  annotate("text",label=c("Oceanic", "Outer", "Middle", "Inner"), x=c(181.5,184,187,192), y=c(58,60.5,61,61))+
  guides(fill = guide_legend(override.aes = list(alpha=0.8)))+
  ggnewscale::new_scale_fill()+
  #geom_sf(data=pribs, aes(fill=island),pch=21)
  ggstar::geom_star(data=pribs,aes(fill=island, x=lon, y=lat), size=3) +
  scale_fill_manual(values=c("yellow","blue"))+
  guides(fill = "none")

#ggsave(file.path(fig.path, "Figure 1 - Complex map.tiff"), GPribs, width=6.5, height=5, dpi=300, compression="lzw")
ggsave(file.path(fig.path, "Figure_1.pdf"), GPribs, width=6.5, height=5, dpi=600,device=cairo_pdf,family="Arial")


# Figure 2 - Predictor variables -------------------------------------------
modelData<-fitSpecies$data[[1]] |>
  dplyr::select(Year, Complex, sstReShelf, bott_s,walleye_pollock_mean_cpue_no) |>
  mutate(walleye_pollock_mean_cpue_no=walleye_pollock_mean_cpue_no/1000)|>
  pivot_longer(-c(Year,Complex))|>
  mutate(Complex=factor(Complex, levels=c("SPEast","SPEnglishBay","SPReefPoint","SGNorth","SGSouth"),
                        labels=c("SP - East", "SP - English Bay", "SP - Reef Point","SG - North", "SG - South")))


GEnviro<-ggplot(data=predVars, aes(x=Year, y=value))+
  geom_point(aes(color=Complex), alpha=0.2)+
  geom_line(aes(color=Complex), alpha=0.4)+
  geom_point(data=modelData, aes(color=Complex), size=2)+
  facet_wrap(~name, scales="free", labeller=as_labeller(c("bott_s"="Bottom temperature","sstReShelf"="Surface temperature","walleye_pollock_mean_cpue_no"="Pollock abundance index")))+
  scale_color_manual(values=complexPal, name="")+
  scale_x_continuous(breaks=seq(1988,2022, by=4))+
  ggthemes::theme_few()+
  geom_vline(xintercept=2012.5, lty=2)+
  theme(legend.position="bottom", axis.text.x=element_text(angle=30, hjust=1,size=12),
        axis.text.y=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=12),strip.text=element_text(size=12))+
  xlab(NULL)+
  ylab("Predictor value")

#ggsave(file.path(fig.path, "Figure_2.tiff"),GEnviro, width=12, height=5, compression="lzw", dpi=300)
ggsave(file.path(fig.path, "Figure_2.pdf"),GEnviro, width=12, height=5, device=cairo_pdf,family="Arial", dpi=600)

# Figure 3 - Year effects ------------------------------------------------------------
yearSmooths<-compare_smooths(fitSpecies$fitG[[1]], fitSpecies$fitG[[2]],fitSpecies$fitG[[3]],fitSpecies$fitG[[4]],fitSpecies$fitG[[5]],fitSpecies$fitG[[6]],fitSpecies$fitG[[7]], select="s(Year)") |>
  mutate(.model=factor(`.model`, labels=c("Gb/Bm", "Gm/Gm","Pacific herring", "Hexagrammid/sablefish","Walleye pollock","Salmon", "Sand lance")))

GYear<-draw(yearSmooths)+
  ggthemes::theme_few()+
  scale_x_continuous(expand=c(0.01,0), breaks=seq(1988, 2012, by=4))+
  geom_rug(data=fitSpecies$data[[1]],sides="b", inherit.aes=F, aes(x=Year,y=0),linewidth=0.25)+
  ggtitle(NULL)+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=12))+
  ylab("Partial effect")+
  scale_color_manual(name="", values=speciesCol[1:7])+
  scale_fill_manual(name="",values=speciesCol[1:7])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~.model, scales="free", nrow=2)+
  theme(legend.position="none")+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)+
  xlab(NULL)

#ggsave(file.path(fig.path, "Figure_3.tiff"),GYear, width=10, height=5.5, compression="lzw", dpi=300)
ggsave(file.path(fig.path, "Figure_3.pdf"),GYear, width=10, height=5.5, device=cairo_pdf,family="Arial", dpi=600)

# On the response scale
yearSmooths2<-fitSpecies |>
  filter(KLPreyGroup %in% fitSum$KLPreyGroup[fitSum$term=="s(Year)"]) |>
  mutate(estimates=map(fitG,~smooth_estimates(.x) |> add_confint()),
         constant=map(fitG, ~coef(.x)[1])) |>
  unnest(c(estimates,constant)) |>
  dplyr::select(-c(fitNullNoYear, fitNull, fitG, sumPoll, sumSurf, sumBott)) |>
  mutate(across(tidyselect::all_of(c(".estimate", ".lower_ci", ".upper_ci")),
                .fns = \(x) linkTransform(x + constant))) |>
  filter(.smooth=="s(Year)") 

GYear2<-ggplot(data=yearSmooths2, aes(y=.estimate, x=Year))+
  geom_line(aes(color=KLPreyGroup))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=KLPreyGroup), alpha=0.2)+
  ggthemes::theme_few()+
  scale_x_continuous(expand=c(0.01,0), breaks=seq(1988, 2012, by=4))+
  geom_rug(data=subset(fo, KLPreyGroup %in% yearSmooths2$KLPreyGroup),sides="b", inherit.aes=F, aes(x=Year,y=NULL),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none",axis.text.x=element_text(angle=30, hjust=1, size=12))+
  labs(y = "Partial effect", x=NULL)+
  scale_color_manual(name="", values=speciesCol[1:7])+
  scale_fill_manual(name="",values=speciesCol[1:7])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~KLPreyGroup, scales="free", nrow=2)

ggsave(file.path(fig.path, "Fig. 3_Year effects_response scale.tiff"),GYear2, width=10, height=5.5, compression="lzw", dpi=300)

# Figure 4c - Pollock effects ------------------------------------------------------------

pollockSmooths<-compare_smooths(fitSpecies$fitG[[3]],fitSpecies$fitG[[4]],fitSpecies$fitG[[5]],fitSpecies$fitG[[6]], select="s(walleye_pollock_mean_cpue_no)")|>
  mutate(.model=factor(`.model`, labels=c("Pacific herring", "Hexagrammid/sablefish","Walleye pollock","Salmon")))

GPollP<-draw(pollockSmooths)+
  ggthemes::theme_few()+
  geom_rug(data=fo,sides="b", inherit.aes=F, aes(x=walleye_pollock_mean_cpue_no,y=0),linewidth=0.25)+
  ggtitle(NULL)+
  theme(axis.text.x=element_text(hjust=1))+
  scale_x_continuous(breaks=seq(0,50000, by=10000), labels=c(0,10,20,30,40,50))+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean pollock abundance (1000s"~km^-2*")"))+
  scale_color_manual(values=speciesCol[c(3:6)])+
  scale_fill_manual(values=speciesCol[c(3:6)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~.model, scales="free", nrow=1)+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)

#ggsave(file.path(fig.path, "Fig. 4_Pollock effects.tiff"),GPollP, width=10, height=3.5, compression="lzw", dpi=300)

# On the response scale
pollockSmooths2<-fitSpecies |>
  filter(KLPreyGroup %in% fitSum$KLPreyGroup[fitSum$term=="s(walleye_pollock_mean_cpue_no)"]) |>
  mutate(estimates=map(fitG,~smooth_estimates(.x) |> add_confint()),
         constant=map(fitG, ~coef(.x)[1])) |>
  unnest(c(estimates,constant)) |>
  dplyr::select(-c(fitNullNoYear, fitNull, fitG, sumPoll, sumSurf, sumBott)) |>
  mutate(across(tidyselect::all_of(c(".estimate", ".lower_ci", ".upper_ci")),
                .fns = \(x) linkTransform(x + constant))) |>
  filter(.smooth=="s(walleye_pollock_mean_cpue_no)")

GPollP2<-ggplot(data=pollockSmooths2, aes(y=.estimate, x=walleye_pollock_mean_cpue_no))+
  geom_line(aes(color=KLPreyGroup))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=KLPreyGroup), alpha=0.2)+
  ggthemes::theme_few()+
  scale_x_continuous(breaks=seq(0,50000, by=10000), labels=c(0,10,20,30,40,50))+
  geom_rug(data=subset(fo, KLPreyGroup %in% pollockSmooths2$KLPreyGroup),sides="b", inherit.aes=F, aes(x=walleye_pollock_mean_cpue_no),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean pollock abundance (1000s"~km^-2*")"))+
  scale_color_manual(values=speciesCol[c(3:6)])+
  scale_fill_manual(values=speciesCol[c(3:6)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~KLPreyGroup, scales="free", nrow=2)

#ggsave(file.path(fig.path, "Fig. 4_Pollock effects_response.tiff"),GPollP2, width=8, height=5.5, compression="lzw", dpi=300)

# Figure 4a - Bottom Temperature effects ------------------------------------------------------------
bottomSmooths<-compare_smooths(fitSpecies$fitG[[2]],fitSpecies$fitG[[3]], fitSpecies$fitG[[4]],fitSpecies$fitG[[5]],fitSpecies$fitG[[6]],select="s(bott_s)") |>
  mutate(.model=factor(`.model`, labels=c("Gm/Gm","Pacific herring", "Hexagrammid/sablefish","Walleye pollock","Salmon")))


GBottP<-draw(bottomSmooths)+
  ggthemes::theme_few()+
  geom_rug(data=fo,sides="b", inherit.aes=F, aes(x=bott_s,y=0),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean bottom temperature ("*degree*"C)"))+
  scale_color_manual(values=speciesCol[c(2:6)])+
  scale_fill_manual(values=speciesCol[c(2:6)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~.model, scales="free", nrow=1)+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)

#ggsave(file.path(fig.path, "Fig. 5_Bottom Temp effects.tiff"),GBottP, width=12, height=3.5, compression="lzw", dpi=300)

# On the response scale
bottSmooths2<-fitSpecies |>
  filter(KLPreyGroup %in% fitSum$KLPreyGroup[fitSum$term=="s(bott_s)"]) |>
  mutate(estimates=map(fitG,~smooth_estimates(.x) |> add_confint()),
         constant=map(fitG, ~coef(.x)[1])) |>
  unnest(c(estimates,constant)) |>
  dplyr::select(-c(fitNullNoYear, fitNull, fitG, sumPoll, sumSurf, sumBott)) |>
  mutate(across(tidyselect::all_of(c(".estimate", ".lower_ci", ".upper_ci")),
                .fns = \(x) linkTransform(x + constant))) |>
  filter(.smooth=="s(bott_s)") 

GBottP2<-ggplot(data=bottSmooths2, aes(y=.estimate, x=bott_s))+
  geom_line(aes(color=KLPreyGroup))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=KLPreyGroup), alpha=0.2)+
  ggthemes::theme_few()+
  geom_rug(data=subset(fo, KLPreyGroup %in% bottSmooths2$KLPreyGroup),sides="b", inherit.aes=F, aes(x=bott_s),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean bottom temperature ("*degree*"C)"))+
  scale_color_manual(values=speciesCol[c(3:6)])+
  scale_fill_manual(values=speciesCol[c(3:6)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~KLPreyGroup, scales="free", nrow=2,
             labeller=as_labeller(c("GBBM"="Gb/Bm","Herring"="Herring","HexSable"=
                                      "Hexagrammid/sablefish","Pollock"="Walleye pollock","Salmon"= "Salmon")))


#ggsave(file.path(fig.path, "Fig. 5_Bottom Temp effects_response scale.tiff"),GBottP2, width=6, height=6, compression="lzw", dpi=300)

# Figure 4b - Surface Temperature effects ------------------------------------------------------------

surfSmooths<-compare_smooths(fitSpecies$fitG[[2]],fitSpecies$fitG[[3]],fitSpecies$fitG[[5]],fitSpecies$fitG[[8]],select="s(sstReShelf)") |>
  mutate(.model=factor(`.model`, labels=c("Gm/Gm","Pacific herring","Walleye pollock","N.smoothtongue")))

GSurfP<-draw(surfSmooths)+
  ggthemes::theme_few()+
  #scale_x_continuous(expand=c(0.01,0), breaks=seq(1988, 2012, by=2))+
  geom_rug(data=fo,sides="b", inherit.aes=F, aes(x=sstReShelf,y=0),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean surface temperature ("*degree*"C)"))+
  scale_color_manual(values=speciesCol[c(2:3,5,8)])+
  scale_fill_manual(values=speciesCol[c(2:3,5,8)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~.model, scales="free", nrow=1)+
  geom_hline(yintercept=0, lty=2, color="gray40", linewidth=0.5)

#ggsave(file.path(fig.path, "Fig. 6_Surface Temp effects.tiff"),GSurfP, width=10, height=3.5, compression="lzw", dpi=300)

# On the response scale
surfSmooths2<-fitSpecies |>
  filter(KLPreyGroup %in% fitSum$KLPreyGroup[fitSum$term=="s(sstReShelf)"]) |>
  mutate(estimates=map(fitG,~smooth_estimates(.x) |> add_confint()),
         constant=map(fitG, ~coef(.x)[1])) |>
  unnest(c(estimates,constant)) |>
  dplyr::select(-c(fitNullNoYear, fitNull, fitG, sumPoll, sumSurf, sumBott)) |>
  mutate(across(tidyselect::all_of(c(".estimate", ".lower_ci", ".upper_ci")),
                .fns = \(x) linkTransform(x + constant))) |>
  filter(.smooth=="s(sstReShelf)") 

GSurfP2<-ggplot(data=surfSmooths2, aes(y=.estimate, x=sstReShelf))+
  geom_line(aes(color=KLPreyGroup))+
  geom_ribbon(aes(ymin=.lower_ci, ymax=.upper_ci, fill=KLPreyGroup), alpha=0.2)+
  ggthemes::theme_few()+
  geom_rug(data=subset(fo, KLPreyGroup %in% surfSmooths2$KLPreyGroup),sides="b", inherit.aes=F, aes(x=sstReShelf),linewidth=0.25)+
  ggtitle(NULL)+
  theme(legend.position="none")+
  labs(y = "Partial effect", x=bquote("Mean surface temperature ("*degree*"C)"))+
  scale_color_manual(values=speciesCol[c(3:6)])+
  scale_fill_manual(values=speciesCol[c(3:6)])+
  guides(fill=guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~KLPreyGroup, scales="free", nrow=2)
  

#ggsave(file.path(fig.path, "Fig. 6_Surface Temp effects.tiff"),GSurfP, width=8, height=6, compression="lzw", dpi=300)


# Figure 4 all - All three plots together ------------------------------------------------

layout <- "
AAAAA
BBBB#
CCCC#
"
GAll<-GBottP/GSurfP/GPollP +
     plot_layout(design = layout) &
    plot_annotation(tag_levels="a") & 
  theme(plot.tag = element_text(size = 12,face="bold"))

#ggsave(file.path(fig.path,"Figure_4.tiff"), GAll, width=12, height=10, dpi=300, compression="lzw")
ggsave(file.path(fig.path,"Figure_4.pdf"), GAll, width=12, height=10,device=cairo_pdf,family="Arial",dpi=600)

# Complex effects (not a figure in manuscript) ------------------------------------------------------------
complexSmooths<-fitSpecies |>
  filter(KLPreyGroup %in% fitSum$KLPreyGroup[fitSum$term=="s(Complex)"]) |>
  mutate(estimates=map(fitG,~smooth_estimates(.x, select="s(Complex)") |> add_confint()),
         constant=map(fitG, ~coef(.x)[1])) |>
  unnest(c(estimates,constant)) |>
  dplyr::select(-c(fitNullNoYear, fitNull, fitG, sumPoll, sumSurf, sumBott)) |>
  mutate(across(tidyselect::all_of(c(".estimate", ".lower_ci", ".upper_ci")),
                .fns = \(x) linkTransform(x + constant)))

# On the response scale
GComplex<-ggplot(data=complexSmooths, aes(y=.estimate, x=Complex))+
  geom_point(aes(color=Complex))+
  geom_pointrange(aes(ymin=.lower_ci, ymax=.upper_ci, color=Complex))+
  ggthemes::theme_few()+
  ggtitle(NULL)+
  theme(legend.position="bottom",axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  labs(y = "Partial effect", x=NULL)+
  scale_color_manual(values=complexPal, name="", labels=c("SP - East","SP - English Bay","SP - Reef Point",
                                                          "SG - North", "SG - South"))+
  facet_grid(~KLPreyGroup,  labeller=as_labeller(c("GBBM"="Gb/Bm","GMGM" ="Gm/Gm","HexSable"=
                                                     "Hexagrammid/sablefish","Pollock"="Walleye pollock","Salmon"= "Salmon","Sandlance" ="Sand lance","Smoothtongue" = "N. smoothtongue")))

ggsave(file.path(fig.path, "Fig. X_Complex effects_response scale.tiff"),GComplex, width=12, height=4, compression="lzw", dpi=300)


