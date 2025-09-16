### Strawberry Guava Project
### Soil analysis

### Manuscript authors: Matthew A. McCary, Zo S. Fenosoa, 
### Julieanne Montaquila, Eric F. Wuesthoff, Emile Rajeriarison,
### Ella Matsuda, Amy E. Dunham

# 8 paired plots in Talatakely and Vohiparara (Sahamalaotra) trail networks of RNP
# Plots 10 x 2 meters and at least 50 m apart from pair, 150 m apart from other pairs
# Sampled June - August, 2024

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(vegan)
library(tidyverse)
library(grid)
library(gridExtra)

#=====import data======
#relative pathname
soil_nutrients <- file.path(".", "Data","guava_soil.csv")

##Load data
soil<-read_csv(soil_nutrients)

#********************************************************
#' Soil Analysis
#********************************************************

soil <- soil %>%
  filter(row_number() %in% c(1:16)) %>%
  mutate(Treatment.1 = case_when(grepl("n", soil$Treatment) ~ "a.NG",
                                 grepl("g", soil$Treatment) ~ "G")
  ) %>%
  mutate(Treatment = Treatment.1)

colnames(soil)<- (c("Site","Pair_number","Treatment","Ntotal",
                    "Ctotal","OrganicMatter","C.N", "Ptotal",
                    "NH4", "N03", "Treatment.1")) 

##Check for normality assumption
#C total
shapiro.test(soil$Ctotal)
hist(soil$Ctotal) 

##N Total
shapiro.test(log(soil$Ntotal))
hist(soil$Ntotal)

#C:N ratio
shapiro.test(soil$C.N)
hist(soil$C.N) 

#Total NH4
shapiro.test(soil$NH4)
hist(soil$NH4) 

#Total NO3
shapiro.test(soil$N03)
hist(soil$N03) 

#Total P_available
shapiro.test(soil$Ptotal)
hist(soil$Ptotal)

#Organic matter
shapiro.test(soil$OrganicMatter)
hist(soil$OrganicMatter)

###linear mixed-effects models

#total C
C_lmm <- lmer(log(Ctotal) ~ Treatment + (1|Pair_number), 
              data = soil)
summary(C_lmm)
anova(C_lmm, ddf = "Kenward-Roger")
plot(C_lmm)

#total N
Ntotal_lmm <- lmer(log(Ntotal) ~ Treatment + (1|Pair_number), 
                   data = soil)
summary(Ntotal_lmm)
anova(Ntotal_lmm, ddf = "Kenward-Roger")
plot(Ntotal_lmm)

#C:N ratio
C.N_lmm <- lmer(log(C.N) ~ Treatment + (1|Pair_number), 
                data = soil)
summary(C.N_lmm)
anova(C.N_lmm, ddf = "Kenward-Roger")
plot(C.N_lmm)

#Total NH4
NH4_lmm <- lmer(log(NH4) ~ Treatment + (1|Pair_number), 
                data = soil)
summary(NH4_lmm)
anova(NH4_lmm, ddf = "Kenward-Roger")
plot(NH4_lmm)

#Total NO3
N03_lmm <- lmer(log(N03) ~ Treatment + (1|Pair_number), 
                data = soil)
summary(N03_lmm)
anova(N03_lmm, ddf = "Kenward-Roger")
summary(N03_lmm)
plot(P_lmm)

#Plant available P
P_lmm <- lmer(log(Ptotal) ~ Treatment + (1|Pair_number), 
              data = soil)
summary(P_lmm)
anova(P_lmm, ddf = "Kenward-Roger")
plot(P_lmm)

#Organic matter
OM_lmm <- lmer(log(OrganicMatter) ~ Treatment + (1|Pair_number), 
               data = soil)
summary(OM_lmm)
anova(OM_lmm, ddf = "Kenward-Roger")
plot(OM_lmm)

#Boxplots to visualize data

#total C
Ctotal_plot <- ggplot(data = soil, aes(Treatment, Ctotal, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Total C"), x = "")+
  annotate(geom = "text", x = 1.5, y = 31, label = "*", cex = 8)+
  annotate(geom = "text", x = 0.6, y = 30, label = "A", cex = 6)

# Total N
Ntotal_plot <- ggplot(data = soil, aes(Treatment, Ntotal, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Total N"), x = "")+
  annotate(geom = "text", x = 1.5, y = 1.5, label = "**", cex = 8)+
  annotate(geom = "text", x = 0.6, y = 1.5, label = "B", cex = 6)

# C:N
C.N_plot <- ggplot(data = soil, aes(Treatment, C.N, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("CN ratio"), x = "")+
  annotate(geom = "text", x = 1.5, y = 23, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.6, y = 22, label = "C", cex = 6)

# NH4
NH4_plot <- ggplot(data = soil, aes(Treatment, NH4, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Total NH4"), x = "")+
  annotate(geom = "text", x = 1.5, y = 290, label = "*", cex = 8)+
  annotate(geom = "text", x = 0.6, y = 286, label = "D", cex = 6)

# N03
N03_plot <- ggplot(data = soil, aes(Treatment, N03, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Total N03"), x = "")+
  annotate(geom = "text", x = 1.5, y = 300, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.6, y = 300, label = "E", cex = 6)

# Total soil phosphorus
Ptotal_plot <- ggplot(data = soil, aes(Treatment, Ptotal, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Plant-available P"), x = "")+
  annotate(geom = "text", x = 1.5, y = 52, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.6, y = 50, label = "F", cex = 6)

# Organic Matter
OM_plot <- ggplot(data = soil, aes(Treatment, OrganicMatter, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 13),
        axis.text.y = element_text(colour= "black", face = "bold", size = 13),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 20, color = 'black', face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "lightgrey", linewidth = 0.2),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5))+
  labs(y = ("Organic Matter"), x = "")+
  annotate(geom = "text", x = 1.5, y = 55, label = "*", cex = 8)+
  annotate(geom = "text", x = 0.6, y = 55, label = "G", cex = 6)

#plot all together
grid.arrange(Ctotal_plot, Ntotal_plot, C.N_plot,NH4_plot,
             N03_plot, Ptotal_plot, OM_plot, ncol = 2)
