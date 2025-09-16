### Strawberry Guava Project
### Invertebrate richness and diversity analyses - pitall, malaise, sticky trap
### This code represents the GLMMs to support the LMMs shown in the main text

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

#=====import data======
#relative pathname
mada_pitfall <- file.path(".", "Data","guava_arthropod_pitfall.csv")
mada_malaise <- file.path(".", "Data", "guava_arthropod_malaise.csv")
mada_sticky <- file.path(".", "Data","guava_arthropod_sticky.csv")

#load arthropod data
pitfall <- read.csv(mada_pitfall, header = T)
pitfall[is.na(pitfall)] <-0
pitfall <- pitfall[-c(39),] #remove bad row

malaise <- read.csv(mada_malaise, header = T)
malaise[is.na(malaise)] <-0

sticky <- read.csv(mada_sticky, header=T)
sticky[is.na(sticky)] <- 0
sticky <-sticky[-c(20,33),] #remove bad rows

#*********************************************
###Pitfall Trap Data ####
#********************************************

#Combine subplot data into one line for each plot####
pitfall_sum <- rowsum(x = pitfall[c(1:38),c(11:47)], group = pitfall$Plot_number)

pitfall_sum <- pitfall_sum %>%
  mutate(Pair = case_when(grepl("1", rownames(pitfall_sum)) ~ "1",
                          grepl("2", rownames(pitfall_sum)) ~ "2",
                          grepl("3", rownames(pitfall_sum)) ~ "3",
                          grepl("4", rownames(pitfall_sum)) ~ "4",
                          grepl("5", rownames(pitfall_sum)) ~ "5",
                          grepl("6", rownames(pitfall_sum)) ~ "6",
                          grepl("7", rownames(pitfall_sum)) ~ "7",
                          grepl("8", rownames(pitfall_sum)) ~ "8")) %>%
  mutate(Treatment.1 = case_when(grepl("N", rownames(pitfall_sum)) ~ "a.NG"))

pitfall_sum <- pitfall_sum %>%
  mutate(Treatment = if_else(is.na(Treatment.1) == TRUE, "G", "a.NG")) 

#Calculate diversity metrics####
pitfall_ren <- renyi(pitfall_sum[1:36], scale = c(0,1,2,Inf), hill = TRUE)
colnames(pitfall_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

pitfall_ren <- merge(pitfall_ren, pitfall_sum, by = 0) %>%
  mutate(pitfall_ren, plotnumber = Row.names) 

###Pitfall abundance and diversity Models####
#Abundance
pitfall.abun_glmm.nb <- glmer.nb(total_abundance ~ Treatment + (1|Pair),
                            data = pitfall_ren)

summary(pitfall.abun_glmm.nb)

# Check overdispersion
check_overdispersion(pitfall.abun_glmm.nb)

#Richness
pitfall.richness_glmm <- glmer(richness ~ Treatment +(1|Pair), data = pitfall_ren,family=poisson(link="log"))
summary(pitfall.richness_glmm)

# Check overdispersion
check_overdispersion(pitfall.richness_glmm)

#Shannons Diversity
pitfall.shannons_glmm <- glmer(shannons ~ Treatment +(1|Pair), data = pitfall_ren,Gamma(link = "log"))
summary(pitfall.shannons_glmm)

# Check overdispersion
simulationOutput_pitfall.shannons_glmm <- simulateResiduals(pitfall.shannons_glmm)
plot(simulationOutput_pitfall.shannons_glmm)

#****************************************
#Malaise Trap Data####
#****************************************

#Combine subplot data into one line for each plot####
malaise_sum <- rowsum(x = malaise[c(1:16),c(9:81)], group = malaise$plot_number)

malaise_sum <- malaise_sum %>%
  mutate(Pair = case_when(grepl("1", rownames(malaise_sum)) ~ "1",
                          grepl("2", rownames(malaise_sum)) ~ "2",
                          grepl("3", rownames(malaise_sum)) ~ "3",
                          grepl("4", rownames(malaise_sum)) ~ "4",
                          grepl("5", rownames(malaise_sum)) ~ "5",
                          grepl("6", rownames(malaise_sum)) ~ "6",
                          grepl("7", rownames(malaise_sum)) ~ "7",
                          grepl("8", rownames(malaise_sum)) ~ "8")) %>%
  mutate(Treatment.1 = case_when(grepl("N", rownames(malaise_sum)) ~ "a.NG"))

malaise_sum <- malaise_sum %>%
  mutate(Treatment = if_else(is.na(Treatment.1) == TRUE, "G", "a.NG")) 

#Calculate diversity metrics####
malaise_ren <- renyi(malaise_sum[1:72], scale = c(0,1,2,Inf), hill = TRUE)
colnames(malaise_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

malaise_ren <- merge( malaise_ren, malaise_sum, by = 0) %>%
  mutate(malaise_ren, plotnumber = Row.names) 

#Abundance and diversity models####
##Malaise models####

#Abundance
malaise_abundance_glmm.nb <- glmer.nb(Total_abundance ~ Treatment+(1|Pair), data = malaise_ren)
summary(malaise_abundance_glmm.nb)

# Check overdispersion
check_overdispersion(malaise_abundance_glmm.nb)

#Richness
####Hill-Shannon diversity of tree seedlings and herbaceous vegetation (x1)
malaise_rich_glmm <- glmer(richness ~ Treatment + (1|Pair), data = malaise_ren,
                   family = Gamma(link = "inverse"))
summary(malaise_rich_glmm)

# Check overdispersion
simulationOutput_malaise_rich_glmm <- simulateResiduals(malaise_rich_glmm)
plot(simulationOutput_malaise_rich_glmm)

#Hill-Shannon
malaise_shan_glmm <- glmer(shannons ~ Treatment + (1|Pair), data = malaise_ren,
                   family = Gamma(link = "inverse"))
summary(malaise_shan_glmm)

#*****************************************************
#Sticky trap data####
#*****************************************************

#Combine subplot data into one line for each plot####
sticky_sum <- rowsum(x = sticky[c(1:31),c(12:59)], group = sticky$Plot_number)

sticky_sum <- sticky_sum %>%
  mutate(Pair = case_when(grepl("1", rownames(sticky_sum)) ~ "1",
                          grepl("2", rownames(sticky_sum)) ~ "2",
                          grepl("3", rownames(sticky_sum)) ~ "3",
                          grepl("4", rownames(sticky_sum)) ~ "4",
                          grepl("5", rownames(sticky_sum)) ~ "5",
                          grepl("6", rownames(sticky_sum)) ~ "6",
                          grepl("7", rownames(sticky_sum)) ~ "7",
                          grepl("8", rownames(sticky_sum)) ~ "8")) %>%
  mutate(Treatment.1 = case_when(grepl("N", rownames(sticky_sum)) ~ "a.NG"))

sticky_sum <- sticky_sum %>%
  mutate(Treatment = if_else(is.na(Treatment.1) == TRUE, "G", "a.NG")) 

#Family-Level Analyses####
#Calculate diversity metrics####
sticky_ren <- renyi(sticky_sum[1:46], scale = c(0,1,2,Inf), hill = TRUE)
colnames(sticky_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

sticky_ren <- merge(sticky_ren, sticky_sum, by = 0) %>%
  mutate(sticky_ren, plotnumber = Row.names) 

#Sticky Models (abundance and diversity)####

#Abundance
sticky_abundance_glmm.nb <- glmer.nb(Total ~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_abundance_glmm.nb)

# Check overdispersion
check_overdispersion(sticky_abundance_glmm.nb)

#richness
sticky_rich_glmm <- glmer(richness ~ Treatment + (1|Pair), data = sticky_ren,
                           family = Gamma(link = "inverse"))
summary(sticky_rich_glmm)

#Hill-Shannon
####Hill-Shannon diversity of tree seedlings and herbaceous vegetation (x1)
sticky_shan_glmm <- glmer(shannons ~ Treatment + (1|Pair), data = sticky_ren,
                   family = Gamma(link = "inverse"))
summary(sticky_shan_glmm)
