library(tidyverse)
library(jagsUI)
library(coda)
library(data.table)

options(scipen = 9999)

#1. Read in data----
dat <- read.csv("CONIMCP_CleanDataAll_Habitat_Roosting_Covs.csv")

#2. Set up loop----
loop <- dat %>% 
  dplyr::select(scale, Season) %>% 
  unique()

#objects to save results
outM <- list()
model.list <- list()
summary.list <- list()
betas.list <- list()
fit.list <- list()

#for(i in 1:nrow(loop)){
for(i in 5:12){  
  
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
  nchoices <- rep(21, nsets)
  
  #a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
  y <- matrix(dat.i$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  #a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
  pt.breed.bird <- dat.i %>% 
    dplyr::select(BirdID, ptID) %>% 
    unique()
  bird <- pt.breed.bird$BirdID
  
  #a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
  X1 <- matrix(dat.i$tree.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X2 <- matrix(dat.i$crops.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X3 <- matrix(dat.i$bare.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  X4 <- matrix(dat.i$water.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
  if(scale.i=="land"){
    X5 <- matrix(dat.i$water2.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
  }

  
  #max value in the bird vector above (i.e., an integer representing how many unique birds there are)
  nbirds <- max(bird)
  
  #5. Model specification----
  
  if(scale.i=="land"){
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

mu.beta4 ~ dnorm(0, 0.01)
tau.beta4 ~ dgamma(10, 0.1)
sigma.beta4 <- 1/sqrt(tau.beta4)  

mu.beta5 ~ dnorm(0, 0.01)
tau.beta5 ~ dgamma(10, 0.1)
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
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j] + beta2[bird[i]]*X2[i,j]+ beta3[bird[i]]*X3[i,j]+ beta4[bird[i]]*X4[i,j] + beta5[bird[i]]*X5[i,j]
    
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
  }
  else {
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

mu.beta4 ~ dnorm(0, 0.01)
tau.beta4 ~ dgamma(10, 0.1)
sigma.beta4 <- 1/sqrt(tau.beta4)  

for(b in 1:nbirds){    
beta1[b] ~ dnorm(mu.beta1, tau.beta1)    
beta2[b] ~ dnorm(mu.beta2, tau.beta2)    
beta3[b] ~ dnorm(mu.beta3, tau.beta3)    
beta4[b] ~ dnorm(mu.beta4, tau.beta4)    
}    

#Likelihood   
    for(i in 1:nsets){    
    y[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    ysim[i,1:nchoices[i]] ~ dmulti(p[i,1:nchoices[i]],1)    
    
    for(j in 1:nchoices[i]){    
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j] + beta2[bird[i]]*X2[i,j]+ beta3[bird[i]]*X3[i,j]+ beta4[bird[i]]*X4[i,j]
    
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
  }
  
  #6. JAGS parameters----
  #Specify data for JAGS
  if(scale.i=="land"){
    win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, nbirds=nbirds)
  }
  else{
    win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, nbirds=nbirds)
  }
  
  #Specify initial values
  inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))
  
  #Specify parameters to track
  if(scale.i=="land"){
    params = c("mu.beta1", "sigma.beta1", "beta1",  
               "mu.beta2", "sigma.beta2", "beta2",           
               "mu.beta3", "sigma.beta3", "beta3",           
               "mu.beta4", "sigma.beta4", "beta4", 
               "mu.beta5", "sigma.beta5", "beta5", 
               "p", "fit.data", "fit.sim", "bpv")
  }
  else{
    params = c("mu.beta1", "sigma.beta1", "beta1",  
               "mu.beta2", "sigma.beta2", "beta2",           
               "mu.beta3", "sigma.beta3", "beta3",           
               "mu.beta4", "sigma.beta4", "beta4", 
               "p", "fit.data", "fit.sim", "bpv")
  }

  
  #Number of chains, iterations, burnin, and thinning 
  nc=3; ni=300000; nb=100000; nt=50; na=1000 
  
  #7. Run JAGS----
  outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T, n.adapt = 1000)
  
  #8. Predictions----
#  outmat <- as.matrix(outM)
  
#  newdat.tree <- expand.grid(X1=seq(min(dat.i$tree.s), max(dat.i$tree.s), 0.1),
#                         X2=mean(dat.i$crops.s),
#                         X3=mean(dat.i$bare.s),
#                         X4=mean(dat.i$water.s),
#                         X5=mean(dat.i$water2.s))
  
#  newdat.crops <- expand.grid(X1=mean(dat.i$tree.s),
#                         X2=seq(min(dat.i$crops.s), max(dat.i$crops.s), 0.1),
#                         X3=mean(dat.i$bare.s),
#                         X4=mean(dat.i$water.s),
#                         X5=mean(dat.i$water2.s))
  
#  newdat.bare <- expand.grid(X1=mean(dat.i$tree.s),
#                         X2=mean(dat.i$bare.s),
#                         X3=seq(min(dat.i$bare.s), max(dat.i$bare.s), 0.1),
#                         X4=mean(dat.i$water.s),
#                         X5=mean(dat.i$water2.s))
  
#  newdat.water <- expand.grid(X1=mean(dat.i$tree.s),
#                             X2=mean(dat.i$bare.s),
#                             X3=mean(dat.i$water.s),
#                             X4=seq(min(dat.i$water.s), max(dat.i$water.s), 0.1),
#                             X5=X4^2)
  
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
  
  if(scale.i=="land"){
    betas.list[[i]] <- data.frame(beta1=outM$sims.list$mu.beta1,
                                  beta2=outM$sims.list$mu.beta2,
                                  beta3=outM$sims.list$mu.beta3,
                                  beta4=outM$sims.list$mu.beta4,
                                  beta5=outM$sims.list$mu.beta5,
                                  scale=scale.i,
                                  season=season.i)
  }
  else{
    betas.list[[i]] <- data.frame(beta1=outM$sims.list$mu.beta1,
                                  beta2=outM$sims.list$mu.beta2,
                                  beta3=outM$sims.list$mu.beta3,
                                  beta4=outM$sims.list$mu.beta4,
                                  beta5=NA,
                                  scale=scale.i,
                                  season=season.i)
  }

  
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
betas <- rbindlist(betas.list) %>% 
  pivot_longer(beta1:beta5, names_to="beta", values_to="value")
fit <- rbindlist(fit.list)

#11. Inspect traceplots
outM <- readRDS("Models/DiscreteChoice_hr_FallMig.RDS")
summary(outM)
jagsUI::traceplot(outM, parameters = c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5"))
jagsUI::traceplot(outM, parameters = c("sigma.beta1", "sigma.beta2", "sigma.beta3", "sigma.beta4", "sigma.beta5"))
#non-convergence: landscape winter, landscape fall (REs only), landscape breed (REs only), hr winter, hr fall (REs only)

#12. Density plots----
betas$cov <- case_when(betas$beta=="beta1" ~ "tree",
                       betas$beta=="beta2" ~ "crop",
                       betas$beta=="beta3" ~ "bare",
                       betas$beta=="beta4" ~ "water",
                       betas$beta=="beta5" ~ "water2")

betas$scale <- factor(betas$scale, levels=c("pt", "hr", "land"))

ggplot(betas %>% 
         dplyr::filter(cov %in% c("crop", "tree", "water"))) +
  geom_density(aes(x=value, colour=season)) +
  geom_vline(aes(xintercept=0)) +
  facet_grid(scale ~ cov, scales="free")

#13. Beta overlap----
overlap <- summary %>% 
  dplyr::filter(val %in% c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4", "mu.beta5"))
View(overlap)

save.image("CONIRoosting_WorkSpace.Rdata")
load("CONIRoosting_WorkSpace.Rdata")
