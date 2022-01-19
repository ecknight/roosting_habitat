library(tidyverse)
library(sf)
library(jagsUI)
library(coda)
library(data.table)

options(scipen = 9999)

#1. Read in tracking data & put together----
raw <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  dplyr::select(ptID, country, region, GCD, Sex, Mass, Population, Lat, Long, Season)

#3. Calculate distance to gulf----
gom <- read_sf("Shapefiles/GulfOfMexico.shp") %>% 
  st_make_valid()

dat.sf <- raw %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326)

dat.dist <- data.frame(distance = as.numeric(st_distance(dat.sf, gom))/1000) %>% 
  cbind(raw) %>% 
  mutate(position = if_else(Lat > 25, "north", "south")) %>% 
  dplyr::select(ptID, country, region, GCD, Sex, Mass, Population, Season, distance, position, Lat, Long) %>% 
  unique() %>% 
  mutate(distance.rd = round(distance, -3),
         distance.rd = ifelse(distance.rd > 4000, 4000, distance.rd))

#Visualize
ggplot() +
  geom_sf(data=gom) +
  geom_point(data=dat.dist, aes(x=Long, y=Lat, colour=distance)) +
  facet_grid(Season ~ position) +
  scale_colour_viridis_c()

#4. Put together with covariate data----
pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID, 
         evi = EVI)

hr <- read.csv("Data/Covariates_hr.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(read.csv("Data/Covariates_hr.csv") %>% 
          dplyr::select(-X, -Y))

dat <- rbind(pt, hr) %>% 
  mutate(scale = Radius,
         scale = factor(scale, levels=c("pt", "hr"))) %>% 
  mutate(water2 = water^2,
         water2.s = water.s^2) %>% 
  left_join(dat.dist) %>% 
  mutate(distance.rd.dir = case_when(Season=="FallMig" & position=="north" ~ -distance.rd,
                                     Season=="SpringMig" & position=="south" ~ -distance.rd,
                                     !is.na(distance.rd) ~ distance.rd),
         distance.dir = case_when(Season=="FallMig" & position=="north" ~ -distance,
                                  Season=="SpringMig" & position=="south" ~ -distance,
                                  !is.na(distance) ~ distance),
         distance.s = scale(distance),
         distance.dir.s = scale(distance.dir))

#5. Visualize some things----
ggplot(dat %>% dplyr::filter(Season %in% c("FallMig", "SpringMig"))) +
  geom_smooth(aes(x=water, y=used, colour=Season)) +
  facet_grid(scale ~ distance.rd.dir, scales="free")

#6.Set up loop----
loop <- dat %>% 
  dplyr::filter(Season %in% c("FallMig", "SpringMig")) %>% 
  dplyr::select(scale, Season) %>% 
  unique()

#objects to save results
outM <- list()
model.list <- list()
summary.list <- list()
betas.list <- list()
fit.list <- list()

for(i in 1:nrow(loop)){
  
  #3. Subset data----
  scale.i <- loop$scale[i]
  season.i <- loop$Season[i]
  
  dat.i <- dat %>% 
    dplyr::filter(Season==season.i,
                  scale==scale.i)
  
  #4. Define model variables----
  
  #number of choice sets
  nsets <- length(unique(dat.i$ptID))
  
  #a vector identifying how many alternatives (including the chosen one) are available for each choice set.
  nchoices <- rep(26, nsets)
  
  #a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
  y <- matrix(dat.i$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  #a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
  pt.breed.bird <- dat.i %>% 
    dplyr::select(BirdID, ptID) %>% 
    unique()
  bird <- pt.breed.bird$BirdID
  
  #a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
  X1 <- matrix(dat.i$evi.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X2 <- matrix(dat.i$distance.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X3 <- matrix(dat.i$evi.s*dat.i$distance.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  #max value in the bird vector above (i.e., an integer representing how many unique birds there are)
  nbirds <- max(bird)
  
  #5. Model specification----
  
  sink("Mixed model.txt")
  cat("model{    
#Priors
mu.beta1 ~ dnorm(0, 0.01)
tau.beta1 ~ dgamma(10, 0.1)
sigma.beta1 <- 1/sqrt(tau.beta1)

mu.beta2 ~ dnorm(0, 0.01)
tau.beta2 ~ dgamma(10, 0.1)
sigma.beta2 <- 1/sqrt(tau.beta2)

mu.beta3 ~ dnorm(0, 0.01)
tau.beta3 ~ dgamma(10, 0.1)
sigma.beta3 <- 1/sqrt(tau.beta3)

for(b in 1:nbirds){    
beta1[b] ~ dnorm(mu.beta1, tau.beta1)    
beta2[b] ~ dnorm(mu.beta2, tau.beta2)    
beta3[b] ~ dnorm(mu.beta3, tau.beta3)    
}    

#Likelihood   
    for(i in 1:nsets){    
    y[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    ysim[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    
    for(j in 1:nchoices[i]){    
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j] + beta2[bird[i]]*X2[i,j]+ beta3[bird[i]]*X3[i,j]
    
    p[i,j] <- phi[i,j]/(sum(phi[i,1:nchoices[i]]))    
    
    D1[i,j] <- pow(y[i,j] - p[i,j], 2)    
    D1sim[i,j] <- pow(ysim[i,j] - p[i,j], 2)    
    }    
    
    D2[i] <- sum(D1[i,1:nchoices[i]])    
    D2sim[i] <- sum(D1sim[i,1:nchoices[i]])
    }    
    
    fit.data <- sum(D2[1:nsets])    
    fit.sim <- sum(D2sim[1:nsets])    
    bpv <- fit.data - fit.sim    
    }",fill = TRUE)
  
  sink()
  
  #6. JAGS parameters----
  #Specify data for JAGS
  win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, nbirds=nbirds)

  #Specify initial values
  inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))
  
  #Specify parameters to track
  params = c("mu.beta1", "sigma.beta1", "beta1",  
             "mu.beta2", "sigma.beta2", "beta2",           
             "mu.beta3", "sigma.beta3", "beta3",           
             "p", "fit.data", "fit.sim", "bpv")
  
  #Number of chains, iterations, burnin, and thinning 
  nc=3; ni=200000; nb=100000; nt=50; na=1000 
  
  #7. Run JAGS----
  outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T, n.adapt = 1000)

  #9. Save outputs----
  model.list[[i]] <- data.frame(DIC = outM$DIC,
                                time = outM$mcmc.info$elapsed.mins,
                                scale = scale.i,
                                season = season.i)
  
  summary.list[[i]] <- outM$summary %>% 
    data.frame() %>% 
    mutate(scale=scale.i,
           season=season.i)
  summary.list[[i]]$val <- row.names(outM$summary)
  
  betas.list[[i]] <- data.frame(beta1=outM$sims.list$mu.beta1,
                                beta2=outM$sims.list$mu.beta2,
                                beta3=outM$sims.list$mu.beta3,
                                scale=scale.i,
                                season=season.i)
  
  fit.list[[i]] <- data.frame(fit.data=outM$sims.list$fit.data,
                              fit.sim=outM$sims.list$fit.sim,
                              bpv=outM$sims.list$bpv,
                              scale=scale.i,
                              season=season.i)
  
  saveRDS(object=outM, file=paste0("Models/Migration/DiscreteChoice_", scale.i,"_", season.i,".RDS"))
  
  print(paste0("Finished model ", i, " of ", nrow(loop), " models in ", outM$mcmc.info$elapsed.mins, " minutes"))
  
}

#10. Collapse outputs----
model <- rbindlist(model.list)
summary <- rbindlist(summary.list)
betas <- rbindlist(betas.list) %>% 
  pivot_longer(beta1:beta3, names_to="beta", values_to="value")
fit <- rbindlist(fit.list)

#12. Beta overlap----
overlap <- summary %>% 
  dplyr::filter(val %in% c("mu.beta1", "mu.beta2", "mu.beta3")) %>% 
  mutate(cov = case_when(val=="mu.beta1" ~ "evi",
                         val=="mu.beta2" ~ "distance",
                         val=="mu.beta3" ~ "evi*distance"))

overlap$scale <- factor(overlap$scale, levels=c("pt", "hr"))

overlap.0 <- overlap %>% 
  dplyr::filter(overlap0==0) %>% 
  arrange(season, scale, cov) %>% 
  dplyr::select(season, scale, cov, mean ,sd, 'X97.5.', 'X2.5.') %>% 
  rename(upper = 'X97.5.',
         lower = 'X2.5.')

table(overlap.0$season, overlap.0$scale)

betas.0 <- betas %>% 
  right_join(overlap.0)

ggplot(betas.0) +
  geom_density(aes(x=value, colour=season)) +
  geom_vline(aes(xintercept=0)) +
  facet_grid(scale ~ cov, scales="free")

save.image("CONIRoosting_WorkSpace_Migration.Rdata")
#load("CONIRoosting_WorkSpace.Rdata")