library(tidyverse)
library(jagsUI)
library(coda)
library(data.table)
library(sf)
library(ggridges)

options(scipen = 9999)

#1. Read in data & put together----
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  dplyr::select(PinpointID, Population) %>% 
  unique() %>% 
  rbind(data.frame(PinpointID = 2217, Population = 8))

pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID, 
         evi = EVI,
         alan = ALAN,
         hmi = HMI,
         scale = Radius) %>% 
  dplyr::select(PinpointID, ptID, scale, used, Season, evi, grass, shrub, patch, alan, hmi, cropdw, waterdw, treedw) %>% 
  left_join(pop)

hr <- read.csv("Data/Covariates_hr.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(read.csv("Data/Covariates_hr.csv") %>% 
          dplyr::select(-X, -Y)) %>% 
  rename(scale = Radius) %>% 
  dplyr::select(PinpointID, ptID, scale, used, Season, evi, grass, shrub, patch, alan, hmi, cropdw, waterdw, treedw) %>% 
  left_join(pop)

#2. Put together and scale-----
dat <- rbind(pt, hr) %>% 
  mutate(grass = grass/100,
         shrub = shrub/100,
         alan = (alan - min(alan))/(max(alan) - min(alan)),
         evi = (evi - min(evi))/(max(evi) - min(evi)))

#2. Visualize covs for polynomials----
ggplot(dat, aes(x=evi, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#Use - no polynomial

ggplot(dat, aes(x=treedw, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#use - consider polynomial

ggplot(dat, aes(x=patch, y=used, colour=scale)) + 
#  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#gah, wish this was running....

ggplot(dat, aes(x=alan, y=used, colour=scale)) + 
  geom_jitter() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#use - consider polynomial

ggplot(dat, aes(x=hmi, y=used, colour=scale)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~Season, scales="free")
#use - no polynomial

ggplot(dat, aes(x=cropdw, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")
#use - no polynomial

ggplot(dat, aes(x=waterdw, y=used, colour=scale)) + 
    geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

#2. Set up loop----
loop <- dat %>% 
  dplyr::select(scale, Season) %>% 
  unique()

#objects to save results
outM <- list()
model.list <- list()
summary.list <- list()
betas.list <- list()
betas.ind.list <- list()
fit.list <- list()

#for(i in 1:nrow(loop)){
for(i in c(1:1)){  
  
  #3. Subset data----
  scale.i <- loop$scale[i]
  season.i <- loop$Season[i]
  
  dat.i <- dat %>% 
    dplyr::filter(Season==season.i,
                  scale==scale.i) %>% 
    arrange(ptID)
  
  #4. Define model variables----
  
  #number of choice sets
  nsets <- length(unique(dat.i$ptID))
  
  #a vector identifying how many alternatives (including the chosen one) are available for each choice set.
  nchoices <- rep(26, nsets)
  
  #a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
  y <- matrix(dat.i$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  #a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
  dat.bird <- dat.i %>% 
    dplyr::select(PinpointID) %>% 
    unique() %>% 
    arrange(PinpointID) %>% 
    mutate(BirdID = row_number()) %>% 
    left_join(dat.i %>% 
                dplyr::select(PinpointID, ptID) %>% 
                unique()) %>% 
    arrange(ptID)
  bird <- dat.bird$BirdID
  
  #a vector that is the same length as the number of choice sets, specifying a numeric population id
  # dat.pop <- dat.i %>% 
  #   dplyr::select(Population) %>% 
  #   unique() %>% 
  #   arrange(Population) %>% 
  #   mutate(PopID = row_number()) %>% 
  #   left_join(dat.i %>% 
  #               dplyr::select(Population, BirdID) %>% 
  #               unique()) %>% 
  #   arrange(BirdID)
  # pop <- dat.pop$PopID
  
  #a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
  X1 <- matrix(dat.i$evi, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X2 <- matrix(dat.i$treedw, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X3 <- matrix(dat.i$alan, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X4 <- matrix(dat.i$waterdw, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X5 <- matrix(dat.i$cropdw, nrow=nsets, ncol=nchoices, byrow=TRUE)

  #max value in the bird vector above (i.e., an integer representing how many unique birds there are)
  nbirds <- max(bird)
#  npops <- max(pop)
  
  #5. Model specification----
  sink("Mixed model.txt")
  cat("model{    
  
#Priors
mu.beta1 ~ dnorm(0, 0.01)
# tau.beta1 ~ dexp(10)
# sigma.beta1 <- 1/sqrt(tau.beta1)
# sigma.beta1 ~ dunif(0, 5)
# tau.beta1 <- pow(sigma.beta1, -2)
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

for(b in 1:nbirds){
beta1[b] ~ dnorm(mu.beta1, tau.beta1)
beta2[b] ~ dnorm(mu.beta2, tau.beta2)    
beta3[b] ~ dnorm(mu.beta3, tau.beta3)    
beta4[b] ~ dnorm(mu.beta4, tau.beta4)   
beta5[b] ~ dnorm(mu.beta5, tau.beta5) 
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

  #6. JAGS parameters----
  #Specify data for JAGS
  win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, nbirds=nbirds)
  
  #Specify initial values
  inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))
  
  #Specify parameters to track
  params = c("mu.beta1", "sigma.beta1",  "beta1",
             "mu.beta2", "sigma.beta2", "beta2",           
             "mu.beta3", "sigma.beta3", "beta3",           
             "mu.beta4", "sigma.beta4", "beta4", 
             "mu.beta5", "sigma.beta5", "beta5", 
             "p", "fit.data", "fit.sim", "bpv")

  #Number of chains, iterations, burnin, and thinning 
  nc=3; ni=250000; nb=50000; nt=50; na=1000 
  
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
                                beta4=outM$sims.list$mu.beta4,
                                beta5=outM$sims.list$mu.beta5,
                                scale=scale.i,
                                season=season.i)
  
  betas.ind.list[[i]] <- rbind(data.frame(outM$sims.list$beta1) %>% 
                                 mutate(beta="beta1"),
                               data.frame(outM$sims.list$beta2) %>% 
                                 mutate(beta="beta2"),
                               data.frame(outM$sims.list$beta3) %>% 
                                 mutate(beta="beta3"),
                               data.frame(outM$sims.list$beta4) %>% 
                                 mutate(beta="beta4"),
                               data.frame(outM$sims.list$beta5) %>% 
                                 mutate(beta="beta5"))
  
  fit.list[[i]] <- data.frame(fit.data=outM$sims.list$fit.data,
                         fit.sim=outM$sims.list$fit.sim,
                         bpv=outM$sims.list$bpv,
                         scale=scale.i,
                         season=season.i)
  
  saveRDS(object=outM, file=paste0("Models/DiscreteChoice_", scale.i,"_", season.i,".RDS"))
  
  print(paste0("Finished model ", i, " of ", nrow(loop), " models in ", outM$mcmc.info$elapsed.mins, " minutes"))
  
}

#10. Collapse outputs----
model <- rbindlist(model.list)
summary <- rbindlist(summary.list)
betas <- rbindlist(betas.list, fill=TRUE) %>% 
  pivot_longer(beta1:beta5, names_to="beta", values_to="value") 
fit <- rbindlist(fit.list) %>% 
   mutate(p = ifelse(bpv>=0, 1, 0))

#11. Save workspace----
#save.image("CONIRoosting_WorkSpace_V2.Rdata")
#load("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/roosting_habitat/CONIRoosting_WorkSpace_V2.Rdata")

#12. Save out traceplots----
for(i in 1:nrow(loop)){
  
  outM <- readRDS(paste0("Models/DiscreteChoice_", loop$scale[i], "_", loop$Season[i], ".RDS"))
   
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

#13. Beta overlap----
overlap <- summary %>% 
  dplyr::filter(val %in%  c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5")) %>% 
  mutate(cov = case_when(val=="mu.beta1" ~ "evi",
                         val=="mu.beta2" ~ "alan",
                         val=="mu.beta3" ~ "tree",
                         val=="mu.beta4" ~ "water",
                         val=="mu.beta5" ~ "crop"))

overlap$scale <- factor(overlap$scale, levels=c("pt", "hr"))

overlap.0 <- overlap %>% 
  dplyr::filter(overlap0==0,
                cov!="water2") %>% 
  arrange(season, scale, cov) %>% 
  dplyr::select(season, scale, cov, mean ,sd, 'X97.5.', 'X2.5.') %>% 
  rename(upper = 'X97.5.',
         lower = 'X2.5.')

table(overlap.0$season, overlap.0$scale)

#14. Density plots----
betas$cov <-  case_when(betas$beta=="beta1" ~ "evi",
                        betas$beta=="beta2" ~ "alan",
                        betas$beta=="beta3" ~ "tree",
                        betas$beta=="beta4" ~ "water",
                        betas$beta=="beta5" ~ "crop")

betas$scale <- factor(betas$scale, levels=c("pt", "hr"))

betas.ci <- betas %>% 
  full_join(overlap) %>% 
  rename(lower = 'X2.5.', upper = 'X97.5.') %>% 
  dplyr::filter(value > lower, value < upper) %>% 
  mutate(cov = as.factor(cov))

write.csv(betas.ci, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/roosting_habitat/betas.csv", row.names = FALSE)

betas.ci <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/roosting_habitat/betas.csv")

ggplot(betas.ci) +
  geom_density_ridges(aes(x=value, y=season, fill=season, alpha=factor(overlap0)), show.legend = FALSE) +
  geom_vline(aes(xintercept=0)) +
  facet_grid(scale ~ cov, scales="free") +
  scale_fill_viridis_d() +
  scale_alpha_manual(values=c(1, 0.3))

#15. Looking for nonconvergence----
summary.rhat <- summary %>% 
  dplyr::filter(Rhat > 1.1)

#Check for betas that didn't converge
summary.rhat.beta <- summary.rhat %>% 
  dplyr::filter(str_sub(val, 1, 4)=="beta")
nrow(summary.rhat.beta)

table(summary.rhat$season, summary.rhat$scale)  

#16. Model fit----
fit.sum <- fit %>% 
  group_by(scale, season) %>% 
  summarize(p = mean(p),
            bpv = mean(bpv)) %>% 
  ungroup()
View(fit.sum)

#17. Prior predictive check----

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

#18. Cov overlap----
overlap.cov <- overlap %>% 
  rename(lower = 'X2.5.', upper = 'X97.5.') %>% 
  dplyr::select(season, scale, cov, mean, lower, upper) %>% 
  arrange(scale, cov, season)

#19. Individual-level betas----
beta1 <- data.frame(outM$sims.list$beta1) %>% 
  mutate(cov="evi") %>% 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta2 <- data.frame(outM$sims.list$beta2) %>% 
  mutate(cov="alan") %>% 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta3 <- data.frame(outM$sims.list$beta3) %>% 
  mutate(cov="tree") %>% 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta4 <- data.frame(outM$sims.list$beta4) %>% 
  mutate(cov="water") %>% 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta5 <- data.frame(outM$sims.list$beta5) %>% 
  mutate(cov="crop") %>% 
  pivot_longer(cols=X1:X23, names_to="bird", values_to = "beta")
beta.id <- rbind(beta1, beta2, beta3, beta4, beta5) %>% 
  mutate(BirdID = as.integer(str_sub(bird, 2, 3))) %>% 
  left_join(dat.bird) %>% 
  left_join(pop)

ggplot(beta.id) +
  geom_density_ridges(aes(x=beta, y=bird, fill=factor(Population))) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
# xlim(c(-10, 10)) +
  facet_wrap(~cov, scales="free_x")
