### Strawberry Guava Project
### Plant community characteristics
### This analysis supports the LMM results shown in the main text

### Manuscript authors: Matthew A. McCary, Zo S. Fenosoa, 
### Julieanne Montaquila, Eric F. Wuesthoff, Emile Rajeriarison,
### Ella Matsuda, Amy E. Dunham

# 8 paired plots in Talatakely and Vohiparara (Sahamalaotra) trail networks of RNP
# Plots 10 x 2 meters and at least 50 m apart from pair, 150 m apart from other pairs
# Sampled June - August, 2024

# Packages
library(tidyverse)
library(lme4)
library(lmerTest)    
library(performance) 
library(DHARMa)      
library(vegan)
library(stringr)
library(tibble)
library(cowplot)
library(rlang)
library(emmeans)

#=====import data======
#relative pathname
mada_canopy <- file.path(".", "Data","guava_canopy.csv")
mada_seedlings <- file.path(".", "Data", "guava_seedlings.csv")
mada_tree_seedlings <- file.path(".", "Data","guava_tree_seedlings.csv")

##Load data
canopy<-read_csv(mada_canopy)
seedlings<-read_csv(mada_seedlings)
tree_seedlings<-read_csv(mada_tree_seedlings)

#*************************************************************
#Canopy characteristics
#*************************************************************

#Import data
# Keep only the one replicate for canopy (Replicate == 1)
cano_open <- canopy %>%
  filter(Replicate == 1) %>%
  mutate(
    Treatment = factor(Treatment),
    Pair      = factor(Pair)
  )

canopy <- canopy %>%
  mutate(
    Treatment = factor(Treatment),
    Pair      = factor(Pair),
    Replicate = factor(Replicate)
  )

#Relevel treatments
cano_open <- cano_open %>%
  mutate(Treatment = relevel(factor(Treatment), ref = "NG"))

canopy <- canopy %>%
  mutate(Treatment = relevel(factor(Treatment), ref = "NG"))

# Quick summaries
cano_open %>% summarise(
  n = n(),
  min_open = min(Canopy_openness, na.rm = TRUE),
  max_open = max(Canopy_openness, na.rm = TRUE)
)

canopy %>% summarise(
  n = n(),
  min_robel = min(Robel_pole, na.rm = TRUE),
  max_robel = max(Robel_pole, na.rm = TRUE)
)

# ==============================================
# Canopy openness
#    Note: Openness is bounded (0–100). 
## ==============================================

# GLMM with Gamma (positive, continuous); inverse link is OK,
# but log link is often easier to interpret.
# Only use if openness behaves like continuous > 0 and heteroskedastic.
mod_open_gamma <- glmer(
  Canopy_openness ~ Treatment + (1 | Pair),
  data   = cano_open,
  family = Gamma(link = "log")
)
summary(mod_open_gamma)

# check overdispersion
simulationOutput_open <- simulateResiduals(mod_open_gamma)
plot(simulationOutput_open)

# ==============================================
# Robel pole
#    Often count-like or discrete non-negative.GLMM (Poisson or NB) with DHARMa.
# ==============================================

# GLMM (Poisson). Check for overdispersion; if present, switch to NB (glmmTMB).
mod_robel_pois <- glmer(
  Robel_pole ~ Treatment + (1 | Pair),
  data   = canopy,
  family = poisson(link = "log")
)
summary(mod_robel_pois)

# check overdispersion
simulationOutput_robel <- simulateResiduals(mod_robel_pois)
plot(simulationOutput_robel)

# ==============================================
# Simple descriptive stats (helpful context)
# ==============================================

# Means by Treatment
canopy %>%
  group_by(Treatment) %>%
  summarise(
    mean_robel = mean(Robel_pole, na.rm = TRUE),
    sd_robel   = sd(Robel_pole, na.rm = TRUE),
    n          = dplyr::n(),
    se_robel   = sd_robel / sqrt(n)
  )

cano_open %>%
  group_by(Treatment) %>%
  summarise(
    mean_open = mean(Canopy_openness, na.rm = TRUE),
    sd_open   = sd(Canopy_openness, na.rm = TRUE),
    n         = dplyr::n(),
    se_open   = sd_open / sqrt(n)
  )

#**********************************************
#Seedling data
#**********************************************

# combine subplot data
seedlings_sum <- rowsum(x = seedlings[c(1:64),c(9:33)], group = seedlings$plotnumber)

seedlings_sum <- seedlings_sum %>%
  mutate(Pair = case_when(grepl("1", rownames(seedlings_sum)) ~ "1",
                          grepl("2", rownames(seedlings_sum)) ~ "2",
                          grepl("3", rownames(seedlings_sum)) ~ "3",
                          grepl("4", rownames(seedlings_sum)) ~ "4",
                          grepl("5", rownames(seedlings_sum)) ~ "5",
                          grepl("6", rownames(seedlings_sum)) ~ "6",
                          grepl("7", rownames(seedlings_sum)) ~ "7",
                          grepl("8", rownames(seedlings_sum)) ~ "8")) %>%
  mutate(Treatment.1 = case_when(grepl("N", rownames(seedlings_sum)) ~ "a.NG"))

seedlings_sum <- seedlings_sum %>%
  mutate(Treatment = if_else(is.na(Treatment.1) == TRUE, "G", "a.NG")) 

seedlings_ren <- renyi(seedlings_sum[1:23], scale = c(0,1,2,Inf), hill = TRUE)

colnames(seedlings_ren) <- (c("X0", "X1","X2","Inf.")) 

seedlings_ren <- merge(seedlings_ren, seedlings_sum, by = 0) %>%
  mutate(seedlings_ren, plotnumber = Row.names) 

# Fit negative binomial GLMM
#total stems
total_glmm <- glmer.nb(Totalstems ~ Treatment + (1|Pair),
                       data = seedlings_ren)
summary(total_glmm)

# Check overdispersion
check_overdispersion(total_glmm)

######Species richness
rich_glmm <- glmer(X0 ~ Treatment + (1|Pair),
                   data = seedlings_ren,
                   family = poisson(link = "log"),
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e5)))

summary(rich_glmm)

# Check overdispersion
check_overdispersion(rich_glmm)

####Hill-Shannon diversity of tree seedlings
shan_glmer<-glmer(X1 ~ Treatment + (1|Pair), data = seedlings_ren,
                  family = Gamma(link = "inverse"))
summary(shan_glmer)

#=======================================================
# Tree seedling data
#=======================================================

#Organize data
tree_sum <- rowsum(x = tree_seedlings[c(1:64),c(9:27)], group = seedlings$plotnumber)

tree_sum <- tree_sum %>%
  mutate(Pair = case_when(grepl("1", rownames(tree_sum)) ~ "1",
                          grepl("2", rownames(tree_sum)) ~ "2",
                          grepl("3", rownames(tree_sum)) ~ "3",
                          grepl("4", rownames(tree_sum)) ~ "4",
                          grepl("5", rownames(tree_sum)) ~ "5",
                          grepl("6", rownames(tree_sum)) ~ "6",
                          grepl("7", rownames(tree_sum)) ~ "7",
                          grepl("8", rownames(tree_sum)) ~ "8")) %>%
  mutate(Treatment.1 = case_when(grepl("N", rownames(tree_sum)) ~ "a.NG"))

tree_sum <- tree_sum %>%
  mutate(Treatment = if_else(is.na(Treatment.1) == TRUE, "G", "a.NG")) 

tree_ren <- renyi(tree_sum[1:16], scale = c(0,1,2,Inf), hill = TRUE)

colnames(tree_ren) <- (c("X0", "X1","X2","Inf.")) 

tree_ren <- merge(tree_ren, tree_sum, by = 0) %>%
  mutate(tree_ren, plotnumber = Row.names) 

#test GLMMs to determine invasion impacts

####Total stems
tree_total_glmm <- glmer.nb(Totalstems ~ Treatment + (1|Pair),
                       data = tree_ren)

summary(tree_total_glmm)

# Check overdispersion
check_overdispersion(tree_total_glmm)

#========================================================
#Native tree seedlings
#=======================================================

# Organize data
natives <- tree_sum[c(2:22)]

native_ren <- renyi(natives[1:15], scale = c(0,1,2,Inf), hill = TRUE)

colnames(native_ren) <- (c("X0", "X1","X2","Inf.")) 

native_ren <- merge(native_ren, natives, by = 0) %>%
  mutate(native_ren, plotnumber = Row.names) 

# Evaluate impacts of guava on native tree seedlings
####Native tree seedling density
native_stem_glmm <- glmer(totalNONguava ~ Treatment + (1|Pair), data = native_ren,
                       family = poisson(link = "log"))
summary(native_stem_glmm) 

# Check overdispersion
check_overdispersion(native_stem_glmm)

####Native tree seedling Hill-Shannon Diversity
#GLMM
native_X1_glmm <- glmer(X1 ~ Treatment + (1|Pair), data = native_ren,
                        family = Gamma(link = "inverse"))
summary(native_X1_glmm) 

