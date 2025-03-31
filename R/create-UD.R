# This script creates utlization distributions for each complex based on satellite telemetry data

library(adehabitatHR)
library(sp)
library(here)

input.path<-file.path(here(), "data-raw")
output.path<-file.path(here(), "data")


# Load in files -----------------------------------------------------------

tracks<-data.table::fread(file.path(input.path, "nfs.csv")) |>
  filter(sex=="female" & ageClass=="Adult" & (island=="St. George" | island=="St. Paul"))|>
  mutate(trip=paste(dbid, tripno, sep="_"))

trackbad<-readxl::read_xlsx(file.path(input.path, "Trips to exclude.xlsx"))
Complex<-readxl::read_xlsx(file.path(input.path, "Rookery complexes.xlsx"))

# Create UD by complex ----------------------------------------------------

tracksSub<-tracks |>
  filter(!(tracks$trip %in% trackbad$trip) & tripno!=999)|>
  dplyr::left_join(Complex,by=c("rookery"="Rookery")) |>
  filter(Complex!="Bogoslof" & month=="August")|>
  mutate(Complex=factor(Complex)) |>
  dplyr::select(mu.x, mu.y, Complex)

afProj<-st_as_sf(tracksSub, coords=c("mu.x", "mu.y"), crs=4326) |>
  st_transform(crs=3338) |>
  as("Spatial")

afUD<-adehabitatHR::kernelUD(afProj[,"Complex"], grid=1000, same4all = T)

afUD95<-adehabitatHR::getverticeshr(afUD, percent=95) #Extract the 95% UD contours

afUD95.sf<-st_as_sf(afUD95)

saveRDS(afUD95.sf, file=file.path(output.path, "Adult female 95% UD.rds"))
#saveRDS(afUD90.sf, file=file.path(output.path, "Adult female 90% UD.rds"))
#saveRDS(afUD50.sf, file=file.path(output.path, "Adult female 50% UD.rds"))
