# ---
# title: CONI roost habitat selection - prediction testing
# author: Elly Knight
# created: November 9, 2021
# updated: June 11, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(brms) # bayesian modelling for dummies

#1. Get the betas----
betas <- read.csv("Results/betas.csv")

set.seed(1)
betas.overall <- betas |> 
  mutate(group = ifelse(season %in% c("Breed", "Winter"), "Stationary", "Migration"),
         group = factor(group, levels=c("Stationary", "Migration")),
         value = abs(value))  |> 
  group_by(season, cov, scale) |> 
  sample_n(1000) |> 
  ungroup()

#2. P1 Overall strength at landscape scale ----

lm.strength <- brm(value ~ group + (1|cov) + (1|scale), data=betas.overall,
             chains = 3, cores = 3, thin=20, control = list(adapt_delta = 0.9),
             warmup = 1000, iter = 5000)
lm.strength

save(lm.strength, file="Results/Predictions/P1SeasonStrength.Rdata")

#3. P2 Variation between stationary seasons & migration seasons ----

betas.stationary <- betas.overall |> 
  dplyr::filter(group=="Stationary")

betas.migration <- betas.overall |> 
  dplyr::filter(group=="Migration")

lm.stationary <- brm(value ~ season + (1|cov) + (1|scale), data=betas.stationary,
                     chains = 3, cores = 3, thin=20, control = list(adapt_delta = 0.9),
                     warmup = 1000, iter = 5000)
lm.stationary

lm.migration <- brm(value ~ season + (1|cov) + (1|scale), data=betas.migration,
                    chains = 3, cores = 3, thin=20, control = list(adapt_delta = 0.9),
                    warmup = 1000, iter = 5000)
lm.migration

save(lm.stationary, lm.migration, file="Results/Predictions/P2SeasonVariation.Rdata")

#4. P3 Strength between scales -----

lm.stationary.scale <- brm(value ~ scale + (1|cov) + (1|season), data=betas.stationary,
                           chains = 3, cores = 3, thin=20, control = list(adapt_delta = 0.9),
                           warmup = 1000, iter = 5000)
lm.stationary.scale

lm.migration.scale <- brm(value ~ scale + (1|cov) + (1|season), data=betas.migration,
                          chains = 3, cores = 3, thin=20, control = list(adapt_delta = 0.9),
                          warmup = 1000, iter = 5000)
lm.migration.scale

save(lm.stationary.scale, lm.migration.scale, file="Results/Predictions/P3ScaleVariation.Rdata")

#5. P4 Correlation between scales ----

#5a. Get data ----
betas.ci <- read.csv("Results/beta_summary.csv")

#5b. Wrangle ----
betas.mn <- betas.ci |> 
  dplyr::select(scale, season, cov, mean) |> 
  unique() |> 
  pivot_wider(names_from=scale, values_from=mean)

#5c. Visualize ----
ggplot(betas.mn) +
  geom_abline(aes(intercept=0, slope=1), linetype="dashed") +
  geom_point(aes(x=pt, y=hr, colour=cov, pch=season), size=3)

#5d. Test ----
cor(betas.mn$pt, betas.mn$hr)
