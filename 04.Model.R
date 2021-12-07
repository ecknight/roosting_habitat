library(tidyverse)
library(jagsUI)
library(corrplot)
library(usdm)

options(scipen=99999)

#TO DO: DECIDE WHETHER TO SCALE COVS####

#PART A: WRANGLING####

dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID) %>% 
  dplyr::select(ptID, Year, Season, Winter)

#1. Wrangle point level----
pt.raw <- read.csv("Data/Copernicus_pt.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(tree.coverfraction))

pt.n <- table(pt.raw$ptID) %>% 
  data.frame() %>% 
  rename(ptID=Var1) %>% 
  dplyr::filter(Freq < 21)
table(pt.n$Freq)

pt <- pt.raw %>% 
  dplyr::filter(!ptID %in% pt.n$ptID) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig")


#Split into seasons, add BirdID field
pt.breed <- pt %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
                dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

pt.winter <- pt %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

pt.fall <- pt %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

pt.spring <- pt %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

#2. Wrangle home range level data----
hr <- read.csv("Data/Copernicus_hr.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(bare.coverfraction)) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig")

hr.breed <- hr %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

hr.winter <- hr %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

hr.fall <- hr %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

hr.spring <- hr %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))


#3. Wrangle landscape level data----
land <- read.csv("Data/Copernicus_land.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(bare.coverfraction)) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig")

land.breed <- land %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

land.winter <- land %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

land.fall <- land %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

land.spring <- land %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used) %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

#PART B: MODELLING####

#4. Model HR level----
#4a. Visualize----
ggplot(hr) +
  geom_hex(aes(x=tree, y=used)) +
  geom_smooth(aes(x=tree, y=used)) +
  facet_wrap(~Season)

ggplot(hr) +
  geom_hex(aes(x=shrub, y=used)) +
  geom_smooth(aes(x=shrub, y=used)) +
  facet_wrap(~Season)

ggplot(hr) +
  geom_hex(aes(x=grass, y=used)) +
  geom_smooth(aes(x=grass, y=used)) +
  facet_wrap(~Season)

ggplot(hr) +
  geom_hex(aes(x=bare, y=used)) +
  geom_smooth(aes(x=bare, y=used)) +
  facet_wrap(~Season)

ggplot(hr) +
  geom_hex(aes(x=crops, y=used)) +
  geom_smooth(aes(x=crops, y=used)) +
  facet_wrap(~Season)

ggplot(hr) +
  geom_hex(aes(x=water, y=used)) +
  geom_smooth(aes(x=water, y=used)) +
  facet_wrap(~Season)

#4b. VIF----
covs.vif <- hr.breed %>% 
  dplyr::select(tree.s, crops.s, bare.s, water.s, grass.s, shrub.s)

M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
vif(covs.vif %>% dplyr::select(-grass.s, -shrub.s))
#should probably take out grass and shrub

#4c. Variables for model----

#number of choice sets
nsets <- length(unique(hr.breed$ptID))
  
#a vector identifying how many alternatives (including the chosen one) are available for each choice set.
nchoices <- rep(21, nsets)
  
#a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
y <- matrix(hr.breed$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
#a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
hr.breed.bird <- hr.breed %>% 
  dplyr::select(BirdID, ptID) %>% 
  unique()
bird <- hr.breed.bird$BirdID
  
#a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
X1 <- matrix(hr.breed$tree.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
X2 <- matrix(hr.breed$crops.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
X3 <- matrix(hr.breed$bare.s, nrow=nsets, ncol=nchoices, byrow=TRUE)
X4 <- matrix(hr.breed$water.s, nrow=nsets, ncol=nchoices, byrow=TRUE)

#max value in the bird vector above (i.e., an integer representing how many unique birds there are)
nbirds <- max(bird)

#4d. Model specification for mixed effects conditional logistic model----

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
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j]    
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

#4e. Specifications for JAGS----
#Specify data for JAGS
win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, nbirds=nbirds)

#Specify initial values
inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))

#Specify parameters to track
params = c("mu.beta1", "sigma.beta1", "beta1",  
           "mu.beta2", "sigma.beta2", "beta2",           
           "mu.beta3", "sigma.beta3", "beta3",           
           "mu.beta4", "sigma.beta4", "beta4", 
           "p", "fit.data", "fit.sim", "bpv")

#Number of chains, iterations, burnin, and thinning 
nc=3; ni=200000; nb=100000; nt=50; na=1000 

#4f. Run the model in JAGS----
outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

#4d. Check model performance----
summary(outM)

#Rhat
outM[["Rhat"]]
#beta 2,3,4 not great

#Traceplots
jagsUI::traceplot(outM, parameters = c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4"))
jagsUI::traceplot(outM, parameters = c("sigma.beta1", "sigma.beta2", "sigma.beta3", "sigma.beta4"))
jagsUI::traceplot(outM, parameters = c("fit.data", "fit.sim", "bpv"))
#same same, esp for sigma

#Density plots
jagsUI::densityplot(outM, parameters = c("mu.beta1", "mu.beta2", "mu.beta3", "mu.beta4"))
jagsUI::densityplot(outM, parameters = c("sigma.beta1", "sigma.beta2", "sigma.beta3", "sigma.beta4"))
jagsUI::densityplot(outM, parameters = c("fit.data", "fit.sim", "bpv"))
#same same, esp for sigma

#Bayesian p value
jagsUI::whiskerplot(outM, parameters = c("fit.data", "fit.sim", "bpv"))
summary(outM[["sims.list"]][["bpv"]])
