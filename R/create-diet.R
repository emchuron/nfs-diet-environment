# This script creates the diet file that is used in the GAMs

library(here)
library(tidyverse)
library(readxl)
library(patchwork)

output.path<-file.path(here(), "data")
input.path<-file.path(here(), "data-raw")
fig.path<-file.path(here::here(), "figures")

source(file.path(here(),"R","diet_fun.R"))

# Read in diet data -------------------------------------------------------

samplecut<-40
taxon<-readxl::read_xlsx(file.path(input.path, "Taxon_Name.xlsx"))
Complex<-readxl::read_xlsx(file.path(input.path, "Rookery complexes.xlsx"))

diet<-readxl::read_xlsx(file.path(input.path,"Species_Composition.xlsx")) |>
  dplyr::left_join(dplyr::select(taxon, -Taxon5)) |>
  dplyr::left_join(dplyr::select(Complex, -rkname))|>
  filter(Island!="Bogoslof" & TaxonCount!=Taxon5Count) |>
  mutate(CollectionDate=as.POSIXct(paste(CollectionYear, CollectionMonth,CollectionDay, sep="-")))|>
  filter(Haulout==0 & TaxonName!=Taxon5)

# Summary statistics ------------------------------------------------------

# Sample sizes by year, island, complex, and collection month
nSum<-diet |>
  group_by(CollectionYear, Island,Complex, CollectionMonth)|>
  dplyr::summarise(nSample=length(unique(SampleID)), nLabel=length(unique(SampleLabel)))

# Sample sizes by year, island, complex, collection month, and sample type
nSum2<-diet |>
  group_by(CollectionYear, Island,Complex, CollectionMonth, SampleType)|>
  dplyr::summarise(nSample=length(unique(SampleID)), nLabel=length(unique(SampleLabel)))

# Proportion of sample types
propType<-diet |>
  group_by(CollectionYear, Island,Complex, CollectionMonth, SampleType)|>
  dplyr::summarise(n=length(unique(SampleLabel))) |>
  group_by(CollectionYear,Island,Complex, CollectionMonth) |>
  dplyr::mutate(N=sum(n), prop=n/N)|>
  ungroup()|>
  filter(N>=samplecut)
  
# Distribution of samples among rookeries
nSumRook<-diet |> 
  group_by(CollectionYear, Island,Complex, Rookery,CollectionMonth)|>
  summarise(n=length(unique(SampleLabel))) |>
  ungroup()|>
  group_by(CollectionYear, Island,Complex, CollectionMonth) |>
  mutate(N=sum(n))|>
  ungroup()

# Calculate FO and identify key prey groups -------------------------------

fo<-freqOcc2(x=diet, groupid="KLPreyGroup", sampleid="SampleLabel", spaceid=c("Complex","Island"),
              timeid=c("CollectionYear","CollectionMonth"), otherid=NULL) |>
  data.table::setnames(c("CollectionYear","CollectionMonth"),c("Year","Month"))

# Summary of FO by prey group - to assess important prey groups
foSum<-fo |>
  group_by(KLPreyGroup, Complex) |>
  filter(N>=samplecut & Month==8)|>
  dplyr::summarise(MeanFO=mean(freq/N), MinFO=min(freq/N),MaxFO=max(freq/N),quantileFO50=quantile(freq/N,0.50),quantileFO75=quantile(freq/N,0.75)) |>
  ungroup() |>
  filter(MeanFO>=0.05) |>
  mutate(Group=interaction(KLPreyGroup, Complex, sep="_"))|>
  arrange(Complex)|>
  filter(KLPreyGroup!="Gadus" & KLPreyGroup!="Gonatus" & KLPreyGroup!="Squid" & KLPreyGroup!="Fish")

# Subset to August, key prey groups, and years with sufficient sample sizes
foSub<-fo |>
 filter(N>=samplecut & Month==8 &KLPreyGroup %in% foSum$KLPreyGroup)

# Raw version of the data (on a by sample basis instead of summarized by year)
foRaw<-diet |>
  group_by(CollectionMonth, CollectionYear,Complex)|>
  mutate(N=length(unique(SampleLabel))) |>
   ungroup()|>
  filter(N>=samplecut & CollectionMonth==8 &KLPreyGroup %in% foSum$KLPreyGroup)|>
  mutate(n=1)|>
  tidyr::complete(SampleLabel,KLPreyGroup, fill=list(n=0))|>
  group_by(SampleLabel)|>
  fill(SampleID,SampleType,Rookery,Haulout,CollectionMonth,CollectionDate,CollectionYear,CollectionDay,Complex,Exclude,Island,N, .direction = 'updown') |>
  ungroup() |>
  data.table::setnames(c("CollectionYear","CollectionMonth"),c("Year","Month"))
  
# Effects of sample size on FO --------------------------------------------

# Look at how sample size affects FO estimates based on repeated resampling. This is how sampleCut size was determined

sampleFO<-foSize(diet, samples=100) |>
  bind_rows()

sampleFOSum<-sampleFO |>
  group_by(KLPreyGroup, Complex, CollectionYear,N)|>
  dplyr::summarise(mFO=mean(FO), sdFO=sd(FO))|>
  ungroup()

saveRDS(sampleFO, file=file.path(output.path, "Resampled data for sample size determination.rds"))

# Save output -------------------------------------------------------------

data.table::fwrite(fo, file.path(output.path, "Diet response variables for analysis_All prey species.csv"))
data.table::fwrite(foSub, file.path(output.path, "Diet response variables for analysis_Focal prey species.csv"))
data.table::fwrite(foRaw, file.path(output.path, "Diet response variables for analysis_Focal prey species_raw.csv"))

# Plots -------------------------------------------------------------------

complexPal<-c("#00346E", "#007CBF", "#06ABDF", "#EDD03E", "#F5A200", "#6D8600", "#424D0C") # from(palettes_d$feathers$bee_eater)


# Sample size effects ----------------------------------------------------------------

GSample1<-ggplot(subset(sampleFO, KLPreyGroup%in% foSub2$KLPreyGroup  & CollectionYear==2002 & KLPreyGroup=="Pollock"), aes(x=N, y=FO))+
  geom_boxplot(outliers=F, aes(group=N, fill=Complex, color=Complex))+
  facet_wrap(~Complex, labeller=labeller(Complex=c("SGNorth"="SG - North", "SGSouth"="SG - South",
                                                   "SPEast"="SP - East", "SPEnglishBay"="SP - English Bay",
                                                   "SPReefPoint"="SP - Reef Point")))+
  geom_vline(xintercept=c(40), lty=2, color="red")+
  scale_y_continuous(lim=c(0,1), label=scales::percent_format())+
  ggthemes::theme_few()+
  geom_point(data=subset(fo, CollectionYear==2002 & KLPreyGroup=="Pollock"), color="red", aes(y=freq/N)) +
  ylab("Frequency of occurrence")+
  xlab("Number of diet samples")+
  scale_fill_manual(values=complexPal[c(4:5,1:3)], name="", labels=c("SG - North", "SG - South", "SP - East","SP - English Bay","SP - Reef Point"))+
  scale_color_manual(values=complexPal[c(4:5,1:3)], name="", labels=c("SG - North", "SG - South", "SP - East","SP - English Bay","SP - Reef Point"))+
  theme(legend.position="none")

GSample2<-ggplot(subset(sampleFO, CollectionYear==2002 & KLPreyGroup=="Salmon"), aes(x=N, y=FO))+
  geom_boxplot(outliers=F, aes(group=N, fill=Complex, color=Complex))+
  facet_wrap(~Complex, labeller=labeller(Complex=c("SGNorth"="SG - North", "SGSouth"="SG - South",
                                                   "SPEast"="SP - East", "SPEnglishBay"="SP - English Bay",
                                                   "SPReefPoint"="SP - Reef Point")))+
  geom_vline(xintercept=c(40), lty=2, color="red")+
  scale_y_continuous(lim=c(0,1),label=scales::percent_format())+
  ggthemes::theme_few()+
  geom_point(data=subset(fo, CollectionYear==2002 & KLPreyGroup=="Salmon"), color="red", aes(y=freq/N))+
  ylab("Frequency of occurrence")+
  xlab("Number of diet samples")+
  scale_fill_manual(values=complexPal[c(4:5,1:3)], name="", labels=c("SG - North", "SG - South", "SP - East","SP - English Bay","SP - Reef Point"))+
  scale_color_manual(values=complexPal[c(4:5,1:3)], name="", labels=c("SG - North", "SG - South", "SP - East","SP - English Bay","SP - Reef Point"))+
  theme(legend.position="none")

GSample<-GSample1/GSample2+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(face="bold", size=12))

ggsave(GSample, file=file.path(fig.path, "Sample size effects on FO.png"), width=11, height=9)

# Summary plots -----------------------------------------------------------

GRook1<-ggplot(subset(nSumRook, CollectionMonth==8), aes(x=CollectionYear, y=n))+
  geom_bar(aes(fill=Rookery), stat="identity")+
  facet_wrap(~Complex, nrow=1)+
  ggthemes::theme_few()+
  theme(axis.text.x=element_text(angle=30, hjust=1))+
  geom_hline(yintercept=40, lty=2)+
  scale_y_continuous(expand=c(0,0))+
  ylab("Sample size")+
  xlab(NULL)

GRook2<-ggplot(subset(nSumRook, CollectionMonth==8), aes(x=CollectionYear, y=n/N))+
  geom_bar(aes(fill=Rookery), stat="identity")+
  facet_wrap(~Complex, nrow=1)+
  ggthemes::theme_few()+
  theme(axis.text.x=element_text(angle=30, hjust=1))+
  scale_y_continuous(labels=scales::percent_format(), lim=c(0,1.02), expand=c(0,0))+
  ylab("Relative contribution")+
  xlab(NULL)

GRook<-GRook1/GRook2+plot_layout(guides = "collect")
#ggsave(file.path(fig.path, "Rookery contributions.png"), GRook, width=12, height=10)

