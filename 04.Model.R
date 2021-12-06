library(tidyverse)
library(Rchoice)
library(jagsUI)

#TO DO: FIGURE OUT FOR CATEGORICAL VAR####
#TO DO: LOOP THROUGH 4 SEASONS####

#PART A: WRANGLING####

dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ID = paste0(PinpointID,"-", row)) %>% 
  arrange(ID) %>% 
  dplyr::select(ID, Year, Season, Winter)

#1. Wrangle point level----
pt <- read.csv("Data/Worldcover_point.csv")
  dplyr::select(-.geo) %>% 
  dplyr::filter(!is.na(mean)) %>% 
  mutate(class = case_when(mean==10 ~ "tree",
                           mean==20 ~ "shrub",
                           mean==30 ~ "grass",
                           mean==40 ~ "crop",
                           mean==50 ~ "urban",
                           mean==60 ~ "bare",
                           mean==70 ~ "snow",
                           mean==80 ~ "water",
                           mean==90 ~ "wetland",
                           mean==95 ~ "mangrove",
                           mean==100 ~ "moss"),
         used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n))

#Select 20 choices per choice set
pt.n <- data.frame(table(pt$ID, pt$used)) %>% 
  rename(ID = Var1, used = Var2) %>% 
  dplyr::filter(used==0, Freq==0)

pt.1 <- pt %>% 
  dplyr::filter(used==1,
                !ID %in% pt.n$ID)

set.seed(1234)
pt.0 <- pt %>% 
  dplyr::filter(used==0) %>% 
  group_by(ID) %>% 
  sample_n(20) %>% 
  ungroup()

pt.20 <- rbind(pt.1, pt.0) 

#Split into seasons, add BirdID field
pt.20.breed <- pt.20 %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt.20 %>% 
                dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ID, used)

pt.20.winter <- pt.20 %>% 
  dplyr::filter(Season=="Winter")

pt.20.fall <- pt.20 %>% 
  dplyr::filter(Season=="FallMig")

pt.20.spring <- pt.20 %>% 
  dplyr::filter(Season=="SpringMig")

#2. Wrangle home range level data----
hr <- read.csv("Data/Copernicus_hr.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(bare.coverfraction)) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n))

hr.breed <- hr %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt.20 %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ID, used)

#Landscape level data
land <- read.csv("Data/Copernicus_land.csv") %>% 
  dplyr::select(-.geo)


#PART B: MODELLING####

#4. Model point level----

#Visualize
ggplot(pt) +
  geom_jitter(aes(x=class, y=used, colour=class)) +
  facet_wrap(~Season)

#Load package


#Variables for model

#number of choice sets
nsets <- length(unique(pt.20.breed$ID))
  
#a vector identifying how many alternatives (including the chosen one) are available for each choice set.
nchoices <- rep(21, nsets)
  
#a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
y <- matrix(pt.20.breed$used, nrow=nsets, ncol=nchoices, byrow=TRUE)
  
#a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
pt.20.breed.bird <- pt.20.breed %>% 
  dplyr::select(BirdID, ID) %>% 
  unique()
bird <- pt.20.breed.bird$BirdID
  
#a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
X1 <- matrix(pt.20.breed$mean, nrow=nsets, ncol=nchoices, byrow=TRUE)

#max value in the bird vector above (i.e., an integer representing how many unique birds there are)
nbirds <- max(bird)

#Model specification for mixed effects conditional logistic model

sink("Mixed model.txt")
cat("model{    
#Priors
mu.beta1 ~ dnorm(0, 0.01)
tau.beta1 ~ dgamma(0.1, 0.1)
sigma.beta1 <- 1/sqrt(tau.beta1)

for(b in 1:nbirds){    
beta1[b] ~ dnorm(mu.beta1, tau.beta1)    
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

#Specify data for JAGS
win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, nbirds=nbirds)

#Specify initial values
inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))

#Specify parameters to track
params = c("mu.beta1", "sigma.beta1", "beta1",        
           "p", "fit.data", "fit.sim", "bpv")

#Number of chains, iterations, burnin, and thinning 
nc=3; ni=200000; nb=100000; nt=50; na=1000 

#Run the model in JAGS
outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)


#5. Model - Home range level----

#Load package
library(jagsUI)

#Variables for model

#number of choice sets
nsets <- length(unique(pt.20.breed$ID))

#a vector identifying how many alternatives (including the chosen one) are available for each choice set.
nchoices <- rep(21, nsets)

#a sets-by-alternatives matrix of 0 (available) and 1(used) values.  So if there are 100 choice sets, each with 20 possible choices, this is a 100-by-20 matrix
y <- matrix(pt.20.breed$used, nrow=nsets, ncol=nchoices, byrow=TRUE)

#a vector that is the same length as the number of choice sets, specifying a numeric bird id.  The smallest number in this vector must be 1, and the largest must be equivalent to the number of unique birds.  
pt.20.breed.bird <- pt.20.breed %>% 
  dplyr::select(BirdID, ID) %>% 
  unique()
bird <- pt.20.breed.bird$BirdID

#a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
X1 <- matrix(pt.20.breed$mean, nrow=nsets, ncol=nchoices, byrow=TRUE)

#max value in the bird vector above (i.e., an integer representing how many unique birds there are)
nbirds <- max(bird)

#Model specification for mixed effects conditional logistic model

sink("Mixed model.txt")
cat("model{    
#Priors
mu.beta1 ~ dnorm(0, 0.01)
tau.beta1 ~ dgamma(0.1, 0.1)
sigma.beta1 <- 1/sqrt(tau.beta1)

for(b in 1:nbirds){    
beta1[b] ~ dnorm(mu.beta1, tau.beta1)    
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

#Specify data for JAGS
win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, nbirds=nbirds)

#Specify initial values
inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1))

#Specify parameters to track
params = c("mu.beta1", "sigma.beta1", "beta1",        
           "p", "fit.data", "fit.sim", "bpv")

#Number of chains, iterations, burnin, and thinning 
nc=3; ni=200000; nb=100000; nt=50; na=1000 

#Run the model in JAGS
outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)



#SAVE FOR LATER WITH MORE COVS####

#a matrix with the same dimensions as y that specifies the first explanatory variable for each set-by-choice combination
X1 <- matrix(pt.20.breed$class, nrow=nsets, ncol=nchoices)
#X2 <- same as X1 (but with a different explanatory variable)
#X3 <- same as X1 (but with a different explanatory variable)
#X4 <- same as X1 (but with a different explanatory variable)
#X5 <- same as X1 (but with a different explanatory variable)
#X6 <- same as X1 (but with a different explanatory variable)

#max value in the bird vector above (i.e., an integer representing how many unique birds there are)
nbirds <- max(bird)

#Model specification for mixed effects conditional logistic model

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
    log(phi[i,j]) <- beta1[bird[i]]*X1[i,j] + beta2[bird[i]]*X2[i,j]+ beta3[bird[i]]*X3[i,j]+ beta4[bird[i]]*X4[i,j]+ beta5[bird[i]]*X5[i,j]+ beta6[bird[i]]*X6[i,j]    
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

#Specify data for JAGS
win.data = list(y=y, nsets=nsets,  nchoices=nchoices, bird=bird, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5, X6=X6, nbirds=nbirds)

#Specify initial values
inits = function()list(mu.beta1=rnorm(1), tau.beta1=runif(1),mu.beta2=rnorm(1), tau.beta2=runif(1),mu.beta3=rnorm(1), tau.beta3=runif(1),mu.beta4=rnorm(1), tau.beta4=runif(1),mu.beta5=rnorm(1), tau.beta5=runif(1),mu.beta6=rnorm(1), tau.beta6=runif(1))

#Specify parameters to track
params = c("mu.beta1", "sigma.beta1", "beta1",        
           "mu.beta2", "sigma.beta2", "beta2",           
           "mu.beta3", "sigma.beta3", "beta3",           
           "mu.beta4", "sigma.beta4", "beta4",           
           "mu.beta5", "sigma.beta5", "beta5",           
           "mu.beta6", "sigma.beta6", "beta6",           
           "p", "fit.data", "fit.sim", "bpv")

#Number of chains, iterations, burnin, and thinning 
nc=3; ni=200000; nb=100000; nt=50; na=1000 

#Run the model in JAGS
outM = jags(win.data, inits, params, "Mixed model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)