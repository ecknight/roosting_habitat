# ---
# title: CONI roost habitat selection - home range & diffuse stopover area estimation
# author: Elly Knight
# created: November 9, 2021
# updated: June 9, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(suncalc) # time since sunrise
library(sf) # working with shps
library(sp) # working with shps and adehabitat
library(ctmm) # spatial autocorrelation
library(adehabitatHR) # home range estimation
library(data.table) # list handling
library(MuMIn) # AIC wrappers
library(lme4) # mixed effects models

#1. Read in data----
ids.exclude <- read.csv("Data/ExclusionYearForBirdsWith2BreedingYears.csv") 
load("Interim/CONIMCP_Habitat.Rdata")

#2. Remove 1 year of breeding for birds with 2 years of breeding data----
dat.hab <- dat.country |> 
  anti_join(ids.exclude |> 
              rename(Year = PtYear)) |> 
  dplyr::mutate(doy = yday(DateTime),
                Season = case_when(PinpointID==2217 & doy %in% c(157:160) ~ "SpringMig",
                                   !is.na(Season) ~ Season))

#3. Identify individuals that don't have enough points for KDE ----
#remove females for breeding season (nesting locations)
dat.use <- dat.hab   |> 
  dplyr::filter(Season %in% c("Breed", "Winter")) |> 
  dplyr::filter(Type != "Band") |> 
  dplyr::filter(!(Sex=="F" & Season=="Breed"))

ids <- data.frame(table(dat.use$PinpointID, dat.use$Season, dat.use$Winter)) |> 
  dplyr::filter(Freq > 5) |> 
  rename(PinpointID=Var1, Season=Var2, Winter=Var3) |> 
  mutate(PinpointID = as.integer(as.character(PinpointID)),
         Winter = as.integer(as.character(Winter)))

#4. Separate data into KDE and not----
#remove banding points for these
dat.kde <- dat.use |> 
  inner_join(ids) |> 
  mutate(ID=paste0(PinpointID, "-", Season, "-", Winter)) |> 
  dplyr::select(ID, PinpointID, Season, Population, Sex, Winter, Long, Lat, DateTime)

#4. Add stopover data----
dat.mig <- dat.hab |> 
  dplyr::filter(Season %in% c("FallMig", "SpringMig")) |> 
  group_by(PinpointID, Season, cluster) |> 
  mutate(count = n()) |> 
  ungroup() |> 
  mutate(kde = ifelse(count >= 5 & !is.na(cluster), 1, 0)) |> 
  mutate(ID = paste0(PinpointID,"-",Season,"-", cluster))

#save this for an appendix figure
write.csv(dat.mig, "Interim/MigrationData.csv", row.names = FALSE)

#5. Put together----
dat.kde.m <- dat.kde |> 
  dplyr::select(ID, Lat, Long, DateTime, Season) |> 
  rbind(dat.mig |> 
          dplyr::filter(kde==1) |> 
          dplyr::select(ID, Lat, Long, DateTime, Season)) |> 
  rename(tag.local.identifier=ID,
         location.lat=Lat,
         location.long=Long,
         timestamp=DateTime) |> 
  unique() |> 
  arrange(tag.local.identifier, timestamp)

write.csv(dat.kde.m, "Interim/KDEData.csv", row.names = FALSE)

#6. Format for ctmm----
datt <- as.telemetry(dat.kde.m, timezone="UTC")

#7. Inspect variograms----
ids <- unique(dat.kde.m$tag.local.identifier)
for(i in 1:length(ids)){
  datt.i <- datt[[ids[i]]]
  SVF <- variogram(datt.i)
  
  jpeg(paste0("Figures/Variograms/Variogram_", ids[i], ".jpeg"))
  plot(SVF)
  dev.off()
  
}

#8. Model & calculate area----
m.sum <- data.frame()
m.area <- data.frame()
for(i in 1:length(ids)){
  m.guess <- ctmm.guess(datt[[ids[i]]], interactive=FALSE)
  m.fits <- ctmm.select(datt[[ids[i]]], m.guess, verbose=TRUE, cores=2)
  m.sum <- data.frame(dAIC=summary(m.fits)[,1],
                      ) |> 
    mutate(mod = row.names(data.frame(summary(m.fits))),
           ID=ids[i])  |> 
    rbind(m.sum)
  m.use <- ctmm.fit(datt[[ids[i]]], m.fits[[1]])
  m.kde <- akde(datt[[ids[i]]], m.use, weights=FALSE)
  m.area <- summary(m.kde)$CI |> 
    data.frame() |> 
    mutate(ID=ids[i]) |> 
    rbind(m.area)
  
  cat(i, "   ")
}

#Save for appendix
write.csv(m.sum, "Results/KDEModelSelection.csv", row.names = FALSE)

#9. Wrangle results----
m.ns <- dat.kde.m |> 
  group_by(tag.local.identifier) |> 
  summarize(n = n()) |> 
  ungroup() |> 
  rename(ID=tag.local.identifier)

m.area.clean <- m.area |> 
  mutate(units = str_sub(row.names(m.area), 7, gregexpr(row.names(m.area), pattern=")"))) |> 
  mutate(est.km = case_when(units=="hectares)" ~ est*.01,
                            units=="square kilometers)" ~ est,
                            units=="square meters)" ~ est*0.000001)) |> 
  separate(ID, into=c("PinpointID", "Season", "cluster"), remove=FALSE) |> 
  mutate(PinpointID = as.numeric(PinpointID)) |> 
  left_join(m.ns) |> 
  left_join(dat.hab |> 
              dplyr::select(PinpointID, Sex) |> 
              unique())

#Save for figures
write.csv(m.area.clean, "Results/KDEArea.csv", row.names = FALSE)

#10. Test for effects of things----
#10a. Stationary seasons----
m.area.stat <- m.area.clean |> 
  dplyr::filter(Season %in% c("Breed", "Winter")) |> 
  group_by(PinpointID, Season) |> 
  sample_n(1) |> 
  ungroup()

#dredge global model
lm.stat <- lm(est.km ~ Season*n + Sex*Season, data=m.area.stat, na.action="na.fail")
dredge.stat <- data.frame(dredge(lm.stat, rank="AICc"))
dredge.stat

#save for appendix
write.csv(dredge.stat, "Results/KDETest_Stationary.csv", row.names = FALSE)

#refit most parsimonious model
lm.stat.use <- lm(est.km ~ Season + (1|PinpointID), data=m.area.stat, na.action="na.fail")
summary(lm.stat.use)

#10b. Migration----
m.area.mig <- m.area.clean |> 
  dplyr::filter(!Season %in% c("Breed", "Winter"))

#dredge global model
lm.mig <- lm(est.km ~ Season*n, data=m.area.mig, na.action="na.fail")
dredge.mig <- data.frame(dredge(lm.mig, rank="AICc"))
dredge.mig #nada

#save for appendix
write.csv(dredge.mig, "Results/KDETest_Migration.csv", row.names = FALSE)

#11. Final summary----
m.area.sum <- m.area.clean |> 
  mutate(Season2 = ifelse(Season %in% c("FallMig", "SpringMig"), "Migration", Season)) |> 
  group_by(Season2) |> 
  summarize(area.mean = mean(est.km),
            area.sd = sd(est.km),
            n = n()) |> 
  mutate(radius.mean = sqrt(area.mean/3.1416),
         radius.sd = sqrt(area.sd/3.1416))
m.area.sum

#Save for covariate extraction
write.csv(m.area.sum, "Interim/KDEAreaMean.csv", row.names = FALSE)

#12. Save out examples for Figure 2----
dat.kde.i <- dat.kde.m |> 
  dplyr::filter(tag.local.identifier %in% c("826-Winter-1", "442-Breed-0", "2217-FallMig-7"))

datt.i <- as.telemetry(dat.kde.i, timezone="UTC")

ids.i <- unique(dat.kde.i$tag.local.identifier)
m.sum <- data.frame()
for(i in 1:length(ids.i)){
  m.guess <- ctmm.guess(datt.i[[ids.i[i]]], interactive=FALSE)
  m.fits <- ctmm.select(datt.i[[ids.i[i]]], m.guess, verbose=TRUE, cores=2)
  m.sum <- data.frame(dAIC=summary(m.fits)[,1]) |> 
    mutate(mod = row.names(data.frame(summary(m.fits))),
           ID=ids.i[i])  |> 
    rbind(m.sum)
  m.use <- ctmm.fit(datt.i[[ids.i[i]]], m.fits[[1]])
  m.akde <- akde(datt.i[[ids.i[i]]],m.use)
  if(i==1){
    shps.i <- ctmm::as.sf(m.akde, level.UD=c(0.05, 0.25, 0.5, 0.75, 0.96))
  }
  else{
    shps.i <- ctmm::as.sf(m.akde, level.UD=c(0.05, 0.25, 0.5, 0.75, 0.96)) |> 
      rbind(shps.i)
  }

}

shps <- shps.i |> 
  separate(name, into=c("id", "iso", "ci"), sep=" ") |> 
  separate(id, into=c("PinpointID", "Season", "id"), remove=FALSE) 

write_sf(shps, "Results/Shapefiles/ExampleKDE.shp")
write.csv(dat.kde.i, "Results/Shapefiles/ExampleKDEData.csv", row.names = FALSE)
