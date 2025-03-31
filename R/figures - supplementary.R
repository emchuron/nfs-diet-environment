# Creates supplementary figures

library(sf)
library(readxl)
library(here)
library(patchwork)


complexPal<-c("#00346E", "#007CBF", "#06ABDF", "#EDD03E", "#F5A200", "#6D8600", "#424D0C") # from(palettes_d$feathers$bee_eater)
complex<-read_xlsx(file.path(here(),"data","Rookery complexes.xlsx")) 
input.path<-file.path(here::here(), "data")


# Complex shapefiles ------------------------------------------------------

spComplex<-st_read(file.path(input.path,"Rookery complexes","SNP_NFSRooks.shp"))|>
  mutate(Complex=c("SPReefPoint","SPReefPoint","SPEast","SPEast","SPEnglishBay",rep("SPEast",times=4), "SPReefPoint",
                   "SPEnglishBay","SPEast",rep("SPEnglishBay", times=2)))

sgComplex<-st_read(file.path(input.path,"Rookery complexes","SNG_NFSRooks.shp"))|>
  left_join(complex)

sg<-st_read(file.path(input.path,"Rookery complexes","SNG.shp"))
sp<-st_read(file.path(input.path,"Rookery complexes","SNP.shp"))

GStGeorge<-ggplot(NULL)+
  geom_sf(data=sgComplex,aes(color=Complex), linewidth=3)+
  geom_sf(data=sg)+
  ggthemes::theme_map()+
  scale_color_manual(values=complexPal[4:5], name="", labels=c("North","South"))+
  theme(legend.position=c(0.8,0.05),plot.background=element_rect("white"),
        legend.text=element_text(size=12))+
  ggspatial::annotation_scale(pad_x=unit(0.15,"cm"), pad_y=unit(0.15,"cm"))+
  ggtitle("St. George Island")

GStPaul<-ggplot(NULL)+
  geom_sf(data=spComplex,aes(color=Complex), linewidth=3)+
  geom_sf(data=sp)+
  ggthemes::theme_map()+
  scale_color_manual(values=complexPal[1:3],name="", labels=c("East","English Bay","Reef Point"))+
  theme(legend.position=c(0.75,0.05), plot.background=element_rect("white"),legend.text=element_text(size=12))+
  ggspatial::annotation_scale(pad_x=unit(0.15,"cm"), pad_y=unit(0.15,"cm"))+
  ggtitle("St. Paul Island")

GRookeries<-GStGeorge/GStPaul
#ggsave(file.path(here(),"figures","Figure_S1.tiff"), GRookeries,width=6, height=7.5, dpi=300,compression="lzw")
ggsave(file.path(here(),"figures","Figure_A1.pdf"), GRookeries,width=6, height=7.5, device=cairo_pdf,family="Arial",dpi=600)
