# ---
# title: CONI roost habitat selection - prediction testing
# author: Elly Knight
# created: November 9, 2021
# updated: June 11, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(lme4) # mixed effects models
library(MuMIn) # dredging

#1. Get the betas----
betas <- read.csv("Results/betas.csv")

#2. P1. Strength between scales for roost variables

# "First, we predicted that selection for variables that represent potential roost camouflage (tree cover, enhanced vegetation index) would be stronger at the local scale than the landscape scale, given that selection is expected to be driven by costs and benefits of the habitat itself that is used for a particular behaviour (Hutto 1985)."

#2a. Visualize ----
ggplot(betas |> dplyr::filter(cov %in% c("tree", "evi"))) + 
  geom_violin(aes(x=season, y=value, fill=scale)) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~cov)

#2b. Filter -----
betas.tree <- betas |> 
  dplyr::filter(cov=="tree") |> 
  mutate(scale = factor(scale, levels=c("pt", "hr")))

betas.evi <- betas |> 
  dplyr::filter(cov=="evi") |> 
  mutate(scale = factor(scale, levels=c("pt", "hr")))

#2c. Test ----
lm.tree <- lmer(value ~ scale + (1|season), data=betas.tree, na.action="na.fail")
dredge(lm.tree)
summary(lm.tree)

lm.evi <- lmer(value ~ scale + (1|season), data=betas.evi, na.action="na.fail")
dredge(lm.evi)
summary(lm.evi)

#3. P2 Strength between seasons at local scale for roost variables

#"We also predicted that the strength of selection for variables that represent potential roost camouflage would not differ across seasons at the local scale, given the importance of roosting for survival across the annual cycle (Chernetsov 2006)."

#3a. Visualize ----
ggplot(betas |> dplyr::filter(cov %in% c("tree", "evi"), scale=="pt")) + 
  geom_violin(aes(x=season, y=value, fill=scale)) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~cov)

#3b. Filter -----
betas.tree.pt <- betas |> 
  dplyr::filter(cov=="tree", scale=="pt")
mutate(season = factor(season, levels=c("Breed", "FallMig", "Winter", "SpringMig")))

betas.evi.pt <- betas |> 
  dplyr::filter(cov=="evi", scale=="pt") |> 
  mutate(season = factor(season, levels=c("Breed", "FallMig", "Winter", "SpringMig")))

#3c. Test ----
lm.tree.pt <- lm(value ~ season, data=betas.tree.pt, na.action="na.fail")
dredge(lm.tree)
summary(lm.tree)

lm.evi.pt <- lm(value ~ season, data=betas.evi.pt, na.action="na.fail")
dredge(lm.evi)
summary(lm.evi)

#4. P3 Overall strength at landscape scale ----

#"Third, we predicted that the overall strength of habitat selection at the landscape scale would be stronger during the stationary seasons than migration seasons, given the shift to more generalist habitat selection during migration and (Petit 2000, Chernetsov 2006, Stanley et al. 2020)."

#4a. Visualize ----
ggplot(betas |> dplyr::filter(scale=="hr") |> 
         mutate(group = ifelse(season %in% c("Breed", "Winter"), "Stationary", "Migration"))) + 
  geom_violin(aes(x=season, y=value, fill=group)) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~cov)

#4b. Filter ----
betas.hr <- betas |> 
  dplyr::filter(scale=="hr") |> 
  mutate(group = ifelse(season %in% c("Breed", "Winter"), "Stationary", "Migration"),
         group = factor(group, levels=c("Stationary", "Migration")),
         value = abs(value))

#4c. Test ---
lm.hr <- lmer(value ~ group + (1|cov), data=betas.hr, na.action="na.fail")
dredge(lm.hr)
summary(lm.hr)
