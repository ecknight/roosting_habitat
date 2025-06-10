# ---
# title: CONI roost habitat selection - discrete choice models for habitat selection
# author: Elly Knight
# created: November 9, 2021
# updated: June 9, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(jagsUI) # jags
library(coda) # model evaluation
library(data.table) # list handling
library(sf) # working with shps
library(usdm) # vif
library(corrplot) # correlation plots
library(ggridges) # visualization

#1. Read in data & put together----
load("Data/CONI_CleanData.Rdata")
pop <- dat |> 
  dplyr::select(PinpointID, Population) |> 
  unique() |> 
  rbind(data.frame(PinpointID = 2217, Population = 8))

pt <- read.csv("Interim/Covariates_local.csv") |> 
  rename(scale = Radius) |> 
  dplyr::select(PinpointID, ptID, scale, used, Season, evi, alan, hmi, crop, water, tree) |> 
  left_join(pop)

hr <- read.csv("Interim/Covariates_landscape.csv") |> 
  st_as_sf(coords=c("X", "Y"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_coordinates() |> 
  cbind(read.csv("Interim/Covariates_landscape.csv") |> 
          dplyr::select(-X, -Y)) |> 
  rename(scale = Radius) |> 
  dplyr::select(PinpointID, ptID, scale, used, Season, evi, alan, hmi, crop, water, tree) |> 
  left_join(pop)

#2. Put together and scale-----
dat <- rbind(pt, hr) |> 
  mutate(alan = (alan - min(alan))/(max(alan) - min(alan)),
         evi = (evi - min(evi))/(max(evi) - min(evi)))

#3. Visualize covs----

#3a. EVI----
ggplot(dat, aes(x=evi, fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=evi, y=used, colour=scale)) + 
  geom_jitter(size=0.5) +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#Use - consider polynomial

#3b. Probability treed----
ggplot(dat, aes(x=tree, fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=tree, y=used, colour=scale)) + 
  geom_jitter(size=0.5) +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#use - no polynomial

#3c. ALAN -----
ggplot(dat, aes(x=log(alan), fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=log(alan), y=used, colour=scale)) + 
  geom_jitter(size=0.5) +
  geom_smooth(method="lm") +
  facet_wrap(~Season, scales="free")
#use - consider polynomial

#3d. Human modification ----
ggplot(dat, aes(x=hmi, fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=hmi, y=used, colour=scale)) + 
#  geom_jitter(size=0.5) +
  geom_smooth(method="lm") +
  facet_wrap(~Season, scales="free")
#use - no polynomial

#3e. Crop probability ----
ggplot(dat, aes(x=crop, fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=crop, y=used, colour=scale)) + 
  geom_jitter(size=0.5) +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#use - no polynomial

#3f. Water probability -----
ggplot(dat, aes(x=log(water), fill=scale)) +
  geom_histogram() + 
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=log(water), y=used, colour=scale)) + 
  geom_jitter(size=0.5) +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

#4. Test for covariance----
loop <- dat |> 
  dplyr::select(scale, Season) |> 
  unique()

vif.out <- data.frame()
corr.out <- data.frame()
for(i in 1:nrow(loop)){
  
  covs <- dat |> 
    dplyr::filter(scale==loop$scale[i],
                  Season==loop$Season[i]) |> 
    dplyr::select(evi, alan, hmi, tree, crop, water)
  
  vif.out <- data.frame(vif(covs)) |> 
                    mutate(scale = loop$scale[i],
                          Season = loop$Season[i]) |> 
    rbind(vif.out)
  
  corr.out <- data.frame(cor(covs)) |> 
    mutate(scale = loop$scale[i],
           Season = loop$Season[i]) |> 
    rbind(corr.out)
  
}

#5. Set up loop----

#objects to save results
outM <- list()
model.list <- list()
summary.list <- list()
betas.list <- list()
betas.ind.list <- list()
fit.list <- list()

for(i in 2:nrow(loop)){
  
  #6. Subset data----
  scale.i <- loop$scale[i]
  season.i <- loop$Season[i]
  
  dat.i <- dat |> 
    dplyr::filter(Season==season.i,
                  scale==scale.i) |> 
    arrange(ptID)
  
  #7. Define model variables----
  
  #number of choice sets
  nsets <- length(unique(dat.i$ptID))
  
  #a vector identifying how many alternatives (including the chosen one) are available for each choice set.
  nchoices <- rep(21, nsets)
  
  #a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
  y <- matrix(dat.i$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  #a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
  dat.bird <- dat.i |> 
    dplyr::select(PinpointID) |> 
    unique() |> 
    arrange(PinpointID) |> 
    mutate(BirdID = row_number()) |> 
    left_join(dat.i |> 
                dplyr::select(PinpointID, ptID) |> 
                unique()) |> 
    arrange(ptID)
  bird <- dat.bird$BirdID
  
  #a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
  X1 <- matrix(dat.i$evi, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X2 <- matrix(dat.i$tree, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X3 <- matrix(dat.i$alan, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X4 <- matrix(dat.i$hmi, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X5 <- matrix(dat.i$water, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X6 <- matrix(dat.i$crop, nrow=nsets, ncol=nchoices, byrow=TRUE)

  #max value in the bird vector above (i.e., an integer representing how many unique birds there are)
  nbirds <- max(bird)
#  npops <- max(pop)
  
  #8. Model specification----
  sink("Mixed model.txt")
  cat("model{    
  
#Priors
mu.beta1 ~ dnorm(0, 0.01)
tau.beta1 ~ dgamma(0.1, 0.1)
sigma.beta1 <- 1/sqrt(tau.beta1)

mu.beta2 ~ dnorm(0, 0.01)
tau.beta2 ~ dgamma(0.1, 0.1)
sigma.beta2 <- 1/sqrt(tau.beta2)

mu.beta3 ~ dnorm(0, 0.01)
tau.beta3 ~ dgamma(0.1, 0.1)
sigma.beta3 <- 1/sqrt(tau.beta3)

mu.beta4 ~ dnorm(0, 0.01)
tau.beta4 ~ dgamma(0.1, 0.1)
sigma.beta4 <- 1/sqrt(tau.beta4)

mu.beta5 ~ dnorm(0, 0.01)
tau.beta5 ~ dgamma(0.1, 0.1)
sigma.beta5 <- 1/sqrt(tau.beta5)

mu.beta6 ~ dnorm(0, 0.01)
tau.beta6 ~ dgamma(0.1, 0.1)
sigma.beta6 <- 1/sqrt(tau.beta6)

for(b in 1:nbirds){
beta1[b] ~ dnorm(mu.beta1, tau.beta1)
beta2[b] ~ dnorm(mu.beta2, tau.beta2)    
beta3[b] ~ dnorm(mu.beta3, tau.beta3)    
beta4[b] ~ dnorm(mu.beta4, tau.beta4)   
beta5[b] ~ dnorm(mu.beta5, tau.beta5)
beta6[b] ~ dnorm(mu.beta6, tau.beta6) 
}    

#Likelihood   
    for(i in 1:nsets){    
    y[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    ysim[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    
    for(j in 1:nchoices[i]){    
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j] + beta2[bird[i]]*X2[i,j]+ beta3[bird[i]]*X3[i,j] + beta4[bird[i]]*X4[i,j] + beta5[bird[i]]*X5[i,j]
    
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

  #9. JAGS parameters----
  #Specify data for JAGS
  win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, nbirds=nbirds)
  
  #Specify initial values
  inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))
  
  #Specify parameters to track
  params = c("mu.beta1", "sigma.beta1",  "beta1",
             "mu.beta2", "sigma.beta2", "beta2",           
             "mu.beta3", "sigma.beta3", "beta3",           
             "mu.beta4", "sigma.beta4", "beta4", 
             "mu.beta5", "sigma.beta5", "beta5", 
             "mu.beta6", "sigma.beta6", "beta6", 
             "p", "fit.data", "fit.sim", "bpv")

  #Number of chains, iterations, burnin, and thinning 
  nc=3; ni=250000; nb=50000; nt=50; na=1000 
  
  #10. Run JAGS----
  outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T, n.adapt = 1000)
  
  #11. Save outputs----
  model.list[[i]] <- data.frame(DIC = outM$DIC,
                           time = outM$mcmc.info$elapsed.mins,
                           scale = scale.i,
                           season = season.i)
  
  summary.list[[i]] <- outM$summary |> 
    data.frame() |> 
    mutate(scale=scale.i,
           season=season.i)
  summary.list[[i]]$val <- row.names(outM$summary)
  
  betas.list[[i]] <- data.frame(beta1=outM$sims.list$mu.beta1,
                                beta2=outM$sims.list$mu.beta2,
                                beta3=outM$sims.list$mu.beta3,
                                beta4=outM$sims.list$mu.beta4,
                                beta5=outM$sims.list$mu.beta5,
                                beta6=outM$sims.list$mu.beta6,
                                scale=scale.i,
                                season=season.i)
  
  betas.ind.list[[i]] <- rbind(data.frame(outM$sims.list$beta1) |> 
                                 mutate(beta="beta1"),
                               data.frame(outM$sims.list$beta2) |> 
                                 mutate(beta="beta2"),
                               data.frame(outM$sims.list$beta3) |> 
                                 mutate(beta="beta3"),
                               data.frame(outM$sims.list$beta4) |> 
                                 mutate(beta="beta4"),
                               data.frame(outM$sims.list$beta5) |> 
                                 mutate(beta="beta5"),
                               data.frame(outM$sims.list$beta6) |> 
                                 mutate(beta="beta6"))
  
  fit.list[[i]] <- data.frame(fit.data=outM$sims.list$fit.data,
                         fit.sim=outM$sims.list$fit.sim,
                         bpv=outM$sims.list$bpv,
                         scale=scale.i,
                         season=season.i)
  
  save(outM, file=paste0("Results/Models/DiscreteChoice_", scale.i,"_", season.i,".Rdata"))
  
  print(paste0("Finished model ", i, " of ", nrow(loop), " models in ", outM$mcmc.info$elapsed.mins, " minutes"))
  
}

#12. Collapse outputs----
model <- rbindlist(model.list)
summary <- rbindlist(summary.list)
betas <- rbindlist(betas.list, fill=TRUE) |> 
  pivot_longer(beta1:beta6, names_to="beta", values_to="value") 
fit <- rbindlist(fit.list) |> 
   mutate(p = ifelse(bpv>=0, 1, 0))

#13. Save workspace----
save.image("CONIRoosting_WorkSpace.Rdata")
load("CONIRoosting_WorkSpace.Rdata")

#14. Save out traceplots----
for(i in 1:nrow(loop)){
  
  load(paste0("Results/Models/DiscreteChoice_", loop$scale[i], "_", loop$Season[i], ".Rdata"))
   
  jpeg(paste0("Figures/Traceplots/Traceplot_beta_",loop$scale[i], "_", loop$Season[i], ".jpeg"))
  jagsUI::traceplot(outM, parameters = c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5"))
  dev.off()
  
  jpeg(paste0("Figures/Traceplots/Densityplot_beta_",loop$scale[i], "_", loop$Season[i], ".jpeg"))
  jagsUI::densityplot(outM, parameters = c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5"))
  dev.off()
  
  jpeg(paste0("Figures/Traceplots/Traceplot_sigma_",loop$scale[i], "_", loop$Season[i], ".jpeg"))
  jagsUI::traceplot(outM, parameters = c("sigma.beta1", "sigma.beta2", "sigma.beta3", "sigma.beta4", "mu.beta5"))
  dev.off()
  
  jpeg(paste0("Figures/Traceplots/Densityplot_sigma_",loop$scale[i], "_", loop$Season[i], ".jpeg"))
  jagsUI::densityplot(outM, parameters = c("sigma.beta1", "sigma.beta2", "sigma.beta3", "sigma.beta4", "mu.beta5"))
  dev.off()
  
  print(paste0("Finished traceplots ", i, " of ", nrow(loop)))
  
}

#15. Beta overlap----
overlap <- summary |> 
  dplyr::filter(val %in%  c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5", "mu.beta6")) |> 
  mutate(cov = case_when(val=="mu.beta1" ~ "evi",
                         val=="mu.beta2" ~ "tree",
                         val=="mu.beta3" ~ "alan",
                         val=="mu.beta4" ~ "hmi",
                         val=="mu.beta5" ~ "water",
                         val=="mu.beta6" ~ "crop"))

overlap$scale <- factor(overlap$scale, levels=c("pt", "hr"))

overlap.0 <- overlap |> 
  dplyr::filter(overlap0==0,
                cov!="water2") |> 
  arrange(season, scale, cov) |> 
  dplyr::select(season, scale, cov, mean ,sd, 'X97.5.', 'X2.5.') |> 
  rename(upper = 'X97.5.',
         lower = 'X2.5.')

table(overlap.0$season, overlap.0$scale)

#16. Density plots----
betas$cov <-  case_when(betas$beta=="beta1" ~ "evi",
                        betas$beta=="beta2" ~ "tree",
                        betas$beta=="beta3" ~ "alan",
                        betas$beta=="beta4" ~ "hmi",
                        betas$beta=="beta5" ~ "water",
                        betas$beta=="beta6" ~ "crop")

betas$scale <- factor(betas$scale, levels=c("pt", "hr"))

betas.ci <- betas |> 
  full_join(overlap) |> 
  rename(lower = 'X2.5.', upper = 'X97.5.') |> 
  dplyr::filter(value > lower, value < upper) |> 
  mutate(cov = as.factor(cov))

write.csv(betas.ci, "Results/betas.csv", row.names = FALSE)
betas.ci <- read.csv("Results/betas.csv")

ggplot(betas.ci) +
  geom_density_ridges(aes(x=value, y=season, fill=season, alpha=factor(overlap0)), show.legend = FALSE) +
  geom_vline(aes(xintercept=0)) +
  facet_grid(scale ~ cov, scales="free") +
  scale_fill_viridis_d() +
  scale_alpha_manual(values=c(1, 0.3))

#17. Looking for nonconvergence----
summary.rhat <- summary |> 
  dplyr::filter(Rhat > 1.1)

#Check for betas that didn't converge
summary.rhat.beta <- summary.rhat |> 
  dplyr::filter(str_sub(val, 1, 4)=="beta")
nrow(summary.rhat.beta)

table(summary.rhat$season, summary.rhat$scale)  

#18. Model fit----
fit.sum <- fit |> 
  group_by(scale, season) |> 
  summarize(p = mean(p),
            bpv = mean(bpv)) |> 
  ungroup()
View(fit.sum)

#19. Prior predictive check----

sink("Mixed model prior check.txt")
cat("model{    
  
#Priors
mu.beta ~ dnorm(0, 0.01)
sigma.beta ~ dunif(0, 5)
tau.beta <- pow(sigma.beta, -2)

for(b in 1:nbirds){    
beta[b] ~ dnorm(mu.beta, tau.beta)    
}    

#Likelihood   
    for(i in 1:nsets){    
    y[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    
    for(j in 1:nchoices[i]){    
    log(phi[i,j]) <- beta[bird[i]]*X1[i,j]
    
    p[i,j] <- phi[i,j]/(sum(phi[i,1:nchoices[i]]))    
    
    D1[i,j] <- pow(y[i,j] - p[i,j], 2)    
    }    
    
    D2[i] <- sum(D1[i,1:nchoices[i]])    
    }    
    
    fit.data <- sum(D2[1:nsets])    
    }",fill = TRUE)

sink()

inits = function()list(mu.beta=rnorm(1), sigma.beta=runif(1))

win.data.check = list(nsets=nsets, nchoices=nchoices, bird=bird, X1=X1, nbirds=nbirds)

params.check = c("mu.beta", "sigma.beta")

outMcheck = jags(win.data.check, inits, params.check, "Mixed model prior check.txt", n.chains=nc, n.thin=nt, n.iter=200000, n.burnin=nb, parallel=T, n.adapt = 1000)

jagsUI::traceplot(outMcheck)
jagsUI::densityplot(outMcheck)

#20. Cov overlap----
overlap.cov <- overlap |> 
  rename(lower = 'X2.5.', upper = 'X97.5.') |> 
  dplyr::select(season, scale, cov, mean, lower, upper) |> 
  arrange(scale, cov, season)

#21. Individual-level betas----
beta1 <- data.frame(outM$sims.list$beta1) |> 
  mutate(cov="evi") |> 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta2 <- data.frame(outM$sims.list$beta2) |> 
  mutate(cov="alan") |> 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta3 <- data.frame(outM$sims.list$beta3) |> 
  mutate(cov="tree") |> 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta4 <- data.frame(outM$sims.list$beta4) |> 
  mutate(cov="water") |> 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta5 <- data.frame(outM$sims.list$beta5) |> 
  mutate(cov="crop") |> 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta.id <- rbind(beta1, beta2, beta3, beta4, beta5) |> 
  mutate(BirdID = as.integer(str_sub(bird, 2, 3))) |> 
  left_join(dat.bird) |> 
  left_join(pop)

ggplot(beta.id) +
  geom_density_ridges(aes(x=beta, y=bird, fill=factor(Population))) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
# xlim(c(-10, 10)) +
  facet_wrap(~cov, scales="free_x")
