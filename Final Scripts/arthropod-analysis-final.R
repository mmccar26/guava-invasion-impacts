### Strawberry Guava Project
### Invertebrate abundance, richness and diversity analyses - pitfall, malaise, and sticky traps

### Manuscript authors: Matthew A. McCary, Zo S. Fenosoa, 
### Julieanne Montaquila, Eric F. Wuesthoff, Emile Rajeriarison,
### Ella Matsuda, Amy E. Dunham

# 8 paired plots in Talatakely and Vohiparara (Sahamalaotra) trail networks of RNP
# Plots 10 x 2 meters and at least 50 m apart from pair, 150 m apart from other pairs
# Sampled June - August, 2024

library(dplyr)
library(lme4)
library(ggplot2)
library(vegan)
library(tidyverse)
library(tidyr)
library(grid)
library(gridExtra)
library(lmerTest)

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

#*****************************************
###Pitfall Trap Data ####
#*****************************************

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

#Abundance and diversity analyses####
#Pitfall Test for normality####

#total abundance
shapiro.test(log(pitfall_ren$total_abundance))
hist(pitfall_ren$total_abundance)

#taxonomic richness
shapiro.test(log(pitfall_ren$richness)) 
hist(pitfall_ren$richness)

#Shannon-Hill
shapiro.test(log(pitfall_ren$shannons))
hist(pitfall_ren$shannons)

###Pitfall LMM Models####
#Abundance
pitfall.abun_lm <- lmer(log(total_abundance) ~ Treatment +(1|Pair), data = pitfall_ren)
summary(pitfall.abun_lm) 
anova(pitfall.abun_lm, F = "Kenward-Rogers")
plot(pitfall.abun_lm)

#Richness
pitfall.richness_lmm <- lmer(log(richness) ~ Treatment +(1|Pair), data = pitfall_ren)
summary(pitfall.richness_lmm)
anova(pitfall.richness_lmm, F = "Kenward-Rogers")
plot(pitfall.richness_lmm)

#Shannon-Hill Diversity
pitfall.shan_lm <- lmer(log(shannons) ~ Treatment +(1|Pair), data = pitfall_ren)
summary(pitfall.shan_lm)
anova(pitfall.shan_lm, F = "Kenward-Rogers")
plot(pitfall.shan_lm)

#Boxplots to visualize data
###Abundance
pitfall_abund_box <- ggplot(data = pitfall_ren, aes(Treatment, total_abundance, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(10,200)+
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
  labs(y = ("Total invertebrate Abundance"), x = "")+
  annotate(geom = "text", x = 1.5, y = 200, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.5, y = 190, label = "A", cex = 6)

#Richness
pitfall_richness_box <- ggplot(data = pitfall_ren, aes(Treatment, richness, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(0,30)+
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
  labs(y = ("Taxonomic richness"), x = "")+
  annotate(geom = "text", x = 1.5, y = 29, label = "*", cex = 8)+
  annotate(geom = "text", x = 0.5, y = 29, label = "B", cex = 6)

####Shannon-Hill Diversity
pitfall_shannons_box <- ggplot(data = pitfall_ren, aes(Treatment, shannons, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(5,12)+
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
  labs(y = ("Hill-Shannon diversity"), x = "")+
  annotate(geom = "text", x = 1.5, y = 12, label = "**", cex = 8)+
  annotate(geom = "text", x = 0.5, y = 11.6, label = "C", cex = 6)

#Pitfall Family Manuscript Figure####
grid.arrange(pitfall_abund_box, pitfall_richness_box, 
             pitfall_shannons_box, ncol = 2)

#Pitfall Trap Functional Groups####
pitfall.fun.groups<-pitfall_sum %>%
  mutate(Predator = Araneae_other+Carabidae+Dolichopodidae,
         Detritivore = Amphipoda+Cryptophagidae+Empididae+Entomobryomorpha+Isopoda+ Meinertellidae+ Muscidae+ Mycetophilidae+
           Phoridae+Poduromorpha+Scarabaeidae+Sciaridae+Sphaeroceridae+Symphypleona+Tupilidae,
         Herbivore=Achilidae+ Bostrichidae+Cecidomyiidae+Gryllidae+Pentatomidae+Tettrigidae,
         Omnivore = Acari+Blattidae+Forfuculidae+Formicidae+Miridae+Mymaridae,
         Parasitoid = Ceraphronidae+Diapriidae+Eupelmidae+Proctotrupidae,
         Other = Ichneumonidae)%>%
  select(Pair,Treatment,Predator,Detritivore,Herbivore,Omnivore,Parasitoid,Other)

pitfall_func_ren <- renyi(pitfall.fun.groups[3:8], scale = c(0,1,2,Inf), hill = TRUE)
colnames(pitfall_func_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

pitfall_func_ren <- merge(pitfall_func_ren, pitfall.fun.groups, by = 0) %>%
  mutate(pitfall_func_ren, plotnumber = Row.names) 

#Pitfall Test for normality (Functional Groups)####
shapiro.test(log(pitfall_func_ren$richness)) 
hist(pitfall_func_ren$richness)

shapiro.test(pitfall_func_ren$shannons) 
hist(pitfall_func_ren$shannons)

#Normality check for functional group analysis
#Predator
shapiro.test(pitfall_func_ren$Predator) 
hist(pitfall_func_ren$Predator)

#Detritivore
shapiro.test(pitfall_func_ren$Detritivore) 
hist(pitfall_func_ren$Detritivore)

#Herbivore
shapiro.test(pitfall_func_ren$Herbivore) 
hist(pitfall_func_ren$Herbivore)

#Omnivore
shapiro.test(pitfall_func_ren$Omnivore) 
hist(pitfall_func_ren$Omnivore)

#Parasitoid
shapiro.test(pitfall_func_ren$Parasitoid)
hist(pitfall_func_ren$Parasitoid)

###Pitfall LMM Models(Functional Groups)####
#detritivore
pitfall_detritivore_lmm <- lmer(Detritivore ~ Treatment+(1|Pair), data = pitfall_func_ren)
summary(pitfall_detritivore_lmm)
anova(pitfall_detritivore_lmm, dff = "Kenward-Roger")
plot(pitfall_detritivore_lmm)
fixef(pitfall_detritivore_lmm)

pitfall_herbivore_lmm <- lmer(Herbivore ~Treatment + (1|Pair), data = pitfall_func_ren)
summary(pitfall_herbivore_lmm)
anova(pitfall_herbivore_lmm, ddf = "Kenward-Roger")
plot(pitfall_herbivore_lmm)
fixef(pitfall_herbivore_lmm)

pitfall_omnivore_lmm <- lmer(Omnivore ~ Treatment+(1|Pair), data = pitfall_func_ren)
summary(pitfall_omnivore_lmm) 
anova(pitfall_omnivore_lmm, ddf = "Kenward-Roger")
plot(pitfall_omnivore_lmm)
fixef(pitfall_omnivore_lmm)

pitfall_parasitoid_lmm <- lmer(Parasitoid ~ Treatment+(1|Pair), data = pitfall_func_ren)
summary(pitfall_parasitoid_lmm)
anova(pitfall_parasitoid_lmm, ddf = "Kenward-Roger")
plot(pitfall_parasitoid_lmm)
fixef(pitfall_parasitoid_lmm)

pitfall_predator_lmm <- lmer(Predator ~ Treatment+(1|Pair), data = pitfall_func_ren)
summary(pitfall_predator_lmm) 
anova(pitfall_predator_lmm, ddf = "Kenward-Roger")
plot(pitfall_predator_lmm)
fixef(pitfall_predator_lmm)

#Permanova and NMDS analysis
# ========== CONFIG ==========
csv_path  <- file.path(".", "Data", "pitfall_ren.csv")
group_var <- "Treatment"  # <- change if needed
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "beige")
point_outline <- "black"

# Optional: enforce legend order
level_order <- c("No Guava", "Guava")  # set to NULL to keep data order
set.seed(123)  # reproducible permutations / NMDS starts
# ===========================

# ---- Load & split ----
dat <- read_csv(csv_path, show_col_types = FALSE) |>
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

stopifnot(group_var %in% names(dat))
group <- factor(pull(dat, all_of(group_var)))

keep <- !is.na(group)
dat  <- dat[keep, , drop = FALSE]
group <- droplevels(group[keep])

comm <- dat |>
  select(where(is.numeric)) |>
  as.data.frame()

if (any(comm < 0, na.rm = TRUE)) stop("Negative values in community matrix.")

# Optional factor level control (for plotting/colors)
if (!is.null(level_order)) group <- factor(group, levels = level_order)

# ---- Transform & distance ----
comm_log <- log1p(as.matrix(comm))
bray_log <- vegdist(comm_log, method = "bray")

# ---- PERMANOVA ----
permanova <- adonis2(bray_log ~ group, permutations = 999, by = "terms")
print(permanova)

# ---- Homogeneity of dispersion (PERMDISP) ----
bd <- betadisper(bray_log, group)
perm_bd <- permutest(bd, permutations = 999)
cat("\n--- Test for homogeneity of multivariate dispersion (betadisper) ---\n")
print(perm_bd)

# ---- NMDS (2D, Bray on log1p data) ----
set.seed(123)
nmds <- metaMDS(comm_log, distance = "bray", k = 2,
                trymax = 200, autotransform = FALSE)
cat(sprintf("\nNMDS stress: %.3f\n", nmds$stress))

scores_df <- as.data.frame(scores(nmds, display = "sites"))

# Match "Second" block: use an Invasion_status column and control order
scores_df$Invasion_status <- as.factor(group)

# ---- Plot: NMDS with treatment ellipses ----
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "bisque2")

p_nmds <- ggplot(scores_df,
                 aes(NMDS1, NMDS2,
                     color = Invasion_status, fill = Invasion_status)) +
  # filled ellipse (95% normal), then outline
  stat_ellipse(type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.30, linewidth = 0) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8) +
  # points
  geom_point(shape = 21, size = 3.2, stroke = 1.0, alpha = 0.9) +
  coord_equal(xlim = c(-1, 1.5), ylim = c(-1, 1)) +
  scale_x_continuous(breaks = seq(-1,1.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  # hide color legend; show fill legend as "Treatment"
  scale_color_manual(values = c("No Guava" = "black", "Guava" = "black"),
                     guide = "none") +
  scale_fill_manual(values = palette_vals, name = "Treatment") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(colour= "black", face = "bold", size = 13),
    axis.text.y = element_text(colour= "black", face = "bold", size = 13),
    axis.line = element_line(colour = "black", size = .3),
    axis.ticks.x = element_line(colour = "black", size = 1),
    axis.ticks.y = element_line(colour = "black", size = 1),
    axis.title = element_text(size=14, face="bold"),
    strip.text = element_text(size = 20, color = 'black', face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    # Legend: bottom-left inside plot
    legend.position = c(0, 0),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(color = "black", size = 1, fill = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(12, "pt"),
    axis.ticks.length = unit(3, "pt")
  ) +
  # Make legend swatches clearer/bigger
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.6, size = 4)
    )
  )

#************************************
#Malaise Trap Data####
#************************************

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

#Abundance and diversity analyses####
#Malaise Test for normality (Family)####
shapiro.test(malaise_ren$Total_abundance)
hist(malaise_ren$Total_abundance)

shapiro.test(malaise_ren$richness)
hist(malaise_ren$richness)

shapiro.test(malaise_ren$shannons)
hist(malaise_ren$shannons)

##Malaise abundance and diversity LLMs####
#Abundance
malaise.abun_lm <- lmer(log(Total_abundance) ~ Treatment+(1|Pair), data = malaise_ren)
summary(malaise.abun_lm)
anova(malaise.abun_lm, ddf = "Kenward-Roger")
plot(malaise.abun_lm)

#Richness
malaise.richness_lm <- lmer(richness ~ Treatment+(1|Pair), data = malaise_ren)
summary(malaise.richness_lm) 
anova(malaise.richness_lm, ddf = "Kenward-Roger")
plot(malaise.richness_lm)

#Hill-Shannon
malaise.shannons_lm <- lmer(shannons ~ Treatment+(1|Pair), data = malaise_ren)
summary(malaise.shannons_lm)
anova(malaise.shannons_lm, ddf = "Kenward-Roger")
plot(malaise.shannons_lm)

#Malaise abundance and diversity plots####
#boxplots to visualize the data

#Abundance
malaise_abund_box <- ggplot(data = malaise_ren, aes(Treatment, Total_abundance, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(100,500)+
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
  labs(y = ("Total invertebrate abundance"), x = "")+
  annotate(geom = "text", x = 1.5, y = 498, label = "*", cex = 8)+
  annotate(geom = "text", x = 0.5, y = 500, label = "A", cex = 6)

###Richness
malaise_richness_box <- ggplot(data = malaise_ren, aes(Treatment, richness, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(15,30)+
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
  labs(y = ("Taxonomic richness"), x = "")+
  annotate(geom = "text", x = 1.5, y = 30, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.5, y = 30, label = "B", cex = 6)

###Shannons
malaise_shannons_box <- ggplot(data = malaise_ren, aes(Treatment, shannons, fill = Treatment))+
  geom_boxplot(aes(), show.legend = FALSE, outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = (c("darkolivegreen", "beige")))+
  scale_x_discrete(label = c("No Guava", "Guava"))+
  ylim(5,20)+
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
  labs(y = ("Hill-Shannon diversity"), x = "")+
  annotate(geom = "text", x = 1.5, y = 20, label = "NS", cex = 4)+
  annotate(geom = "text", x = 0.5, y = 20, label = "C", cex = 6)

#Malaise Family Manuscript figure####
grid.arrange(malaise_abund_box, malaise_richness_box, malaise_shannons_box, 
             ncol = 2)

#Malaise Trap Functional Groups####
malaise.fun.groups<-malaise_sum %>%
  mutate(Predator = Araneae_other+Carabidae+Dolicopodidae+Histeridae+Reduviidae+Saldidae+Staphilinidae,
         Detritivore = Amphipoda+Cryptophagidae+Empididae+Entomobryomorpha+Isopoda+Meinertellidae+Mycetophilidae+
           Phoridae+Poduromorpha+Psycodidae+Sciaridae+Sphaeroceridae+Symphypleona+Tupilidae+Tephritidae+Thripidae,
         Herbivore=Achilidae+Anostostomatidae+Apidae+Bostrichidae+Caterpillars+Cecidomyiidae+Cerambycidae+Cercopidae+
         Chrysomelidae+Cicadellidae+Crambidae+Curculionidae+Cynipidae+Drosophilidae+Ectopsocidae+Elateridae+Erebidae+
         Faniidae+Geometridae+Grylidae+Kalotermitidae+Mordelidae+Muscidae+Nictuidae+Nitidulidae+Phasmantidae+
         Syrphidae+Tettrigidae,
         Omnivore = Acari+Blattidae+Forfuculidae+Formicidae+Stratiomidae,
         Parasite = Ceratopogonidae+Culusidae,
         Parasitoid = Bethylidae+Ceraphronidae+Diaprudae+Eucharitidae+Eupelmidae+Mymaridae+Proctotrupidae+Tiphiidae,
         Other = Ichneumonidae+Myrmeleontidae+Caeciliusidae)%>%
  select(Pair,Treatment,Predator,Detritivore,Herbivore,Omnivore,Parasite,Parasitoid,Other)

malaise_func_ren <- renyi(malaise.fun.groups[3:9], scale = c(0,1,2,Inf), hill = TRUE)

colnames(malaise_func_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

malaise_func_ren <- merge(malaise_func_ren, malaise.fun.groups, by = 0) %>%
  mutate(malaise_func_ren, plotnumber = Row.names) 

#Malaise Test for normality (Functional Groups)####
#richness
shapiro.test(malaise_func_ren$richness)
hist(malaise_func_ren$richness)

#Hill-Shannon
shapiro.test(malaise_func_ren$shannons) 
hist(malaise_func_ren$shannons)

###Start here for functional group abundances
shapiro.test(malaise_func_ren$Predator) 
hist(malaise_func_ren$Predator)

shapiro.test(malaise_func_ren$Detritivore) 
hist(malaise_func_ren$Detritivore)

shapiro.test(malaise_func_ren$Herbivore) 
hist(malaise_func_ren$Herbivore)

shapiro.test(malaise_func_ren$Omnivore) 
hist(malaise_func_ren$Omnivore)

shapiro.test(malaise_func_ren$Parasite) 
hist(malaise_func_ren$Parasite)

shapiro.test(malaise_func_ren$Parasitoid) 
hist(malaise_func_ren$Parasitoid)

shapiro.test(malaise_func_ren$Other) 
hist(malaise_func_ren$Other)

###Malaise Models (Functional Groups)####
#Detritivore
malaise_detritivore_lm <- lmer(log(Detritivore) ~Treatment+(1|Pair), data = malaise_func_ren)
summary(malaise_detritivore_lm)
anova(malaise_detritivore_lm, ddf = "Kenward-Roger")
plot(malaise_detritivore_lm)
fixef(malaise_detritivore_lm)

#Herbivore
malaise_herbivore_lm <- lmer(log(Herbivore) ~Treatment+(1|Pair), data = malaise_func_ren)
summary(malaise_herbivore_lm) 
anova(malaise_herbivore_lm, ddf = "Kenward-Roger")
fixef(malaise_herbivore_lm)

#Omnivore
malaise_omnivore_lm <- lmer(log(Omnivore+1) ~Treatment+(1|Pair), data = malaise_func_ren)
summary(malaise_omnivore_lm)
anova(malaise_omnivore_lm, ddf = "Kenward-Roger")
plot(malaise_omnivore_lm)
fixef(malaise_omnivore_lm)

#Parasitoid
malaise_parasitoid_lm <- lmer(log(Parasitoid+1) ~ Treatment+(1|Pair), data = malaise_func_ren)
summary(malaise_parasitoid_lm)
anova(malaise_parasitoid_lm, ddf = "Kenward-Roger")
plot(malaise_parasitoid_lm)
fixef(malaise_parasitoid_lm)

#Predator
malaise_predator_lm <- lmer(log(Predator+1) ~Treatment+(1|Pair), data = malaise_func_ren)
summary(malaise_predator_lm)
anova(malaise_predator_lm, ddf = "Kenward-Roger")
plot(malaise_predator_lm)
fixef(malaise_predator_lm)

##Permanova and NMDS analysis
# ========== CONFIG ==========
csv_path  <- file.path(".", "Data", "malaise_ren.csv")
group_var <- "Treatment"  
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "beige")
point_outline <- "black"

# Optional: enforce legend order
level_order <- c("No Guava", "Guava")  # set to NULL to keep data order
set.seed(123)  # reproducible permutations / NMDS starts
# ===========================

# ---- Load & split ----
dat <- read_csv(csv_path, show_col_types = FALSE) |>
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

stopifnot(group_var %in% names(dat))
group <- factor(pull(dat, all_of(group_var)))

keep <- !is.na(group)
dat  <- dat[keep, , drop = FALSE]
group <- droplevels(group[keep])

comm <- dat |>
  select(where(is.numeric)) |>
  as.data.frame()

if (any(comm < 0, na.rm = TRUE)) stop("Negative values in community matrix.")

# Optional factor level control (for plotting/colors)
if (!is.null(level_order)) group <- factor(group, levels = level_order)

# ---- Transform & distance ----
comm_log <- log1p(as.matrix(comm))
bray_log <- vegdist(comm_log, method = "bray")

# ---- PERMANOVA ----
permanova <- adonis2(bray_log ~ group, permutations = 999, by = "terms")
print(permanova)

# ---- Homogeneity of dispersion (PERMDISP) ----
bd <- betadisper(bray_log, group)
perm_bd <- permutest(bd, permutations = 999)
cat("\n--- Test for homogeneity of multivariate dispersion (betadisper) ---\n")
print(perm_bd)

# ---- NMDS (2D, Bray on log1p data) ----
set.seed(123)
nmds <- metaMDS(comm_log, distance = "bray", k = 2,
                trymax = 200, autotransform = FALSE)
cat(sprintf("\nNMDS stress: %.3f\n", nmds$stress))

scores_df <- as.data.frame(scores(nmds, display = "sites"))

# Match "Second" block: use an Invasion_status column and control order
scores_df$Invasion_status <- as.factor(group)

# ---- Plot: NMDS with treatment ellipses ----
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "bisque2")

m_nmds <- ggplot(scores_df,
                 aes(NMDS1, NMDS2,
                     color = Invasion_status, fill = Invasion_status)) +
  # filled ellipse (95% normal), then outline
  stat_ellipse(type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.30, linewidth = 0) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8) +
  # points
  geom_point(shape = 21, size = 3.2, stroke = 1.0, alpha = 0.9) +
  coord_equal(xlim = c(-1, 1.5), ylim = c(-1, 1)) +
  scale_x_continuous(breaks = seq(-1,1.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  # hide color legend; show fill legend as "Treatment"
  scale_color_manual(values = c("No Guava" = "black", "Guava" = "black"),
                     guide = "none") +
  scale_fill_manual(values = palette_vals, name = "Treatment") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(colour= "black", face = "bold", size = 13),
    axis.text.y = element_text(colour= "black", face = "bold", size = 13),
    axis.line = element_line(colour = "black", size = .3),
    axis.ticks.x = element_line(colour = "black", size = 1),
    axis.ticks.y = element_line(colour = "black", size = 1),
    axis.title = element_text(size=14, face="bold"),
    strip.text = element_text(size = 20, color = 'black', face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    # Legend: bottom-left inside plot
    legend.position = c(0, 0),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(color = "black", size = 1, fill = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(12, "pt"),
    axis.ticks.length = unit(3, "pt")
  ) +
  # Make legend swatches clearer/bigger
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.6, size = 4)
    )
  )

#*****************************************
#Sticky trap data####
#*****************************************

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

#Abundance and Diversity Level Analyses####
#Calculate diversity metrics####
sticky_ren <- renyi(sticky_sum[1:46], scale = c(0,1,2,Inf), hill = TRUE)
colnames(sticky_ren) <- (c("richness", "shannons","inv_simpson","weighted_div")) 

sticky_ren <- merge(sticky_ren, sticky_sum, by = 0) %>%
  mutate(sticky_ren, plotnumber = Row.names) 

#Sticky Test for normality (Abundance and Diversity)

#total abundance
shapiro.test(sticky_ren$Total)
hist(sticky_ren$Total)

#Richness
shapiro.test(sticky_ren$richness) 
hist(sticky_ren$richness)

#hill-Shannon
shapiro.test(sticky_ren$shannons)
hist(sticky_ren$shannons)

#Sticky Models (Abundance and Diversity)####

#abundance
sticky_abundance_lm <- lmer(log(Total)~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_abundance_lm) 
anova(sticky_abundance_lm, ddf = "Kenward-Roger")
plot(sticky_abundance_lm)

#richness
sticky_richness_lm <- lmer(richness~ Treatment +(1|Pair), data = sticky_ren)
summary(sticky_richness_lm)
anova(sticky_richness_lm, ddf = "Kenward-Roger")
plot(sticky_richness_lm)

#shannons
sticky_shannons_lm <- lmer(shannons~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_shannons_lm)
anova(sticky_shannons_lm, ddf = "Kenward-Roger")
plot(sticky_shannons_lm)

#Sticky Trap Functional Groups####
sticky.fun.groups<-sticky_sum %>%
  mutate(Predator = Anthocoridae+Dolichpodidae..dip.+Mesostigmata+Prostigmata+Reduviidae..Hemipt.+
         Salticidae+Staphylinidae,
         Detritivore = Amphipods+Diopsidae..Dipt..stalk.eyed.flies.+Drosophilidae..dipt.+
         Empididae+Entomobryomorpha+ Meinertellidae..thysanura.+ Muscidae..dipt.+
         Mycetophilidae+ Oribatida+Phoridae..dipt.+Poduromorpha+Sciaridae+Sphaeroceridae..dipt.+
         Symphypleona+Tephritidae+Typilidae, 
         Herbivore=Braconidae..hymenop.+cecidomyidae..dip.+Chrysomelidae..coleop.+
           Cicadellidae..orthopt.+Cixiidae..hemipt.+Cynipidae..hymenopt.asked.scott.to.ID.+Faniidae..dipt.+
           Gryllidae..orthoptera. +Kinnaridae..Himiptera. +Tetrigidae..orthoptera.,
         Omnivore = Forficulidae..dermopt.+Formicidae+Psychodidae..dipt.,
         Parasite = Culicidae..dipt..mosquitos. ,
         Parasitoid = Ceraphronidae..hymenopt.+Diapriidae..hymenopt.+Eupelmidae..hymenopt. +
           Mymaridae..hymenoptera.+Scelionidae..hymenopt.,
         Other = Asteridae +Diptera.family.unknown +Ichneumonidae..hymenopt.)%>%
  select(Pair,Treatment,Predator,Detritivore,Herbivore,Omnivore,Parasite,Parasitoid,Other)

sticky_func_sum <- merge(sticky_sum, sticky.fun.groups, by = 0)

#Sticky Test for normality (Functional Groups)####

#Predator
shapiro.test(sticky_func_sum$Predator)
hist(sticky_func_sum$Predator)

#Detritivore
shapiro.test(sticky_func_sum$Detritivore) 
hist(sticky_func_sum$Detritivore)

#Herbivore
shapiro.test(sticky_func_sum$Herbivore) 
hist(sticky_func_sum$Herbivore)

#Omnivore
shapiro.test(sticky_func_sum$Omnivore) 
hist(sticky_func_sum$Omnivore)

#Parasite
shapiro.test(sticky_func_sum$Parasite) 
hist(sticky_func_sum$Parasite)

#Parasitoid
shapiro.test(sticky_func_sum$Parasitoid) 

#Sticky Models (Functional Groups)####

#Detritivore
sticky_detritivore_lm <- lmer(log(Detritivore) ~ Treatment.x+(1|Pair.x), data = sticky_func_sum)
summary(sticky_detritivore_lm) 
anova(sticky_detritivore_lm, ddf = "Kenward-Roger")
plot(sticky_detritivore_lm)
fixef(sticky_detritivore_lm)

#Herbivore
sticky_herbivore_lmm <- lmer(log(Herbivore) ~ Treatment.x +(1|Pair.x), data = sticky_func_sum)
summary(sticky_herbivore_lmm) 
anova(sticky_herbivore_lmm, ddf = "Kenward-Roger")
plot(sticky_herbivore_lmm)
fixef(sticky_herbivore_lmm)

#Omnivore
sticky_omnivore_lmm <- lmer(log(Omnivore+1) ~ Treatment.x +(1|Pair.x), data = sticky_func_sum)
summary(sticky_omnivore_lmm)
anova(sticky_omnivore_lmm, ddf = "Kenward-Roger")
plot(sticky_omnivore_lmm)
fixef(sticky_omnivore_lmm)

#Parasitoid
sticky_parasitoid_lm <- lmer(log(Parasitoid+1) ~ Treatment.x+(1|Pair.x), data = sticky_func_sum)
summary(sticky_parasitoid_lm) 
anova(sticky_parasitoid_lm, ddf= "Kenward-Roger")
plot(sticky_parasitoid_lm)
fixef(sticky_parasitoid_lm)

#Predator
sticky_predator_lm <- lmer(log(Predator+1) ~ Treatment.x+(1|Pair.x), data = sticky_func_sum)
summary(sticky_predator_lm) 
anova(sticky_predator_lm, ddf = "Kenward-Roger")
plot(sticky_predator_lm)
fixef(sticky_predator_lm)

#Individual taxa
#Poduromorpha
sticky_pod_lmm <- lmer(log(Poduromorpha+1)~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_pod_lmm)
anova(sticky_pod_lmm, ddf = "Kenward-Roger")
plot(sticky_pod_lmm)
t.test(Poduromorpha~Treatment, data = sticky_ren)

#Symphypleona
sticky_sym_lmm <- lmer(log(Symphypleona+1)~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_sym_lmm)
anova(sticky_sym_lmm, ddf = "Kenward-Roger")
plot(sticky_sym_lmm)

#Entomobryomorpha
sticky_ent_lmm <- lmer(log(Entomobryomorpha+1)~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_ent_lmm)
anova(sticky_ent_lmm, ddf = "Kenward-Roger")
plot(sticky_sym_lmm)

#total collembola
sticky_coll_lmm <- lmer(log(Total.Collembola)~ Treatment+(1|Pair), data = sticky_ren)
summary(sticky_coll_lmm)
anova(sticky_coll_lmm, ddf = "Kenward-Roger")
plot(sticky_coll_lmm)

###Permanova and NMDS analysis
# ========== CONFIG ==========
csv_path  <- file.path(".", "Data", "sticky_ren.csv")
group_var <- "Treatment"  # <- change if needed
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "beige")
point_outline <- "black"

# Optional: enforce legend order
level_order <- c("No Guava", "Guava")  # set to NULL to keep data order
set.seed(123)  # reproducible permutations / NMDS starts
# ===========================

# ---- 1) Load & split ----
dat <- read_csv(csv_path, show_col_types = FALSE) |>
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

stopifnot(group_var %in% names(dat))
group <- factor(pull(dat, all_of(group_var)))

keep <- !is.na(group)
dat  <- dat[keep, , drop = FALSE]
group <- droplevels(group[keep])

comm <- dat |>
  select(where(is.numeric)) |>
  as.data.frame()

if (any(comm < 0, na.rm = TRUE)) stop("Negative values in community matrix.")

# Optional factor level control (for plotting/colors)
if (!is.null(level_order)) group <- factor(group, levels = level_order)

# ---- 2) Transform & distance ----
comm_log <- log1p(as.matrix(comm))
bray_log <- vegdist(comm_log, method = "bray")

# ---- 3) PERMANOVA ----
permanova <- adonis2(bray_log ~ group, permutations = 999, by = "terms")
print(permanova)

# ---- 4) Homogeneity of dispersion (PERMDISP) ----
bd <- betadisper(bray_log, group)
perm_bd <- permutest(bd, permutations = 999)
cat("\n--- Test for homogeneity of multivariate dispersion (betadisper) ---\n")
print(perm_bd)

# ---- 5) NMDS (2D, Bray on log1p data) ----
set.seed(123)
nmds <- metaMDS(comm_log, distance = "bray", k = 2,
                trymax = 200, autotransform = FALSE)
cat(sprintf("\nNMDS stress: %.3f\n", nmds$stress))

scores_df <- as.data.frame(scores(nmds, display = "sites"))

# Match "Second" block: use an Invasion_status column and control order
scores_df$Invasion_status <- as.factor(group)

# ---- 6) Plot: NMDS with treatment ellipses ----
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "bisque2")

s_nmds <- ggplot(scores_df,
                 aes(NMDS1, NMDS2,
                     color = Invasion_status, fill = Invasion_status)) +
  # filled ellipse (95% normal), then outline
  stat_ellipse(type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.30, linewidth = 0) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8) +
  # points
  geom_point(shape = 21, size = 3.2, stroke = 1.0, alpha = 0.9) +
  coord_equal(xlim = c(-1, 1.5), ylim = c(-1, 1)) +
  scale_x_continuous(breaks = seq(-1,1.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  # hide color legend; show fill legend as "Treatment"
  scale_color_manual(values = c("No Guava" = "black", "Guava" = "black"),
                     guide = "none") +
  scale_fill_manual(values = palette_vals, name = "Treatment") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(colour= "black", face = "bold", size = 13),
    axis.text.y = element_text(colour= "black", face = "bold", size = 13),
    axis.line = element_line(colour = "black", size = .3),
    axis.ticks.x = element_line(colour = "black", size = 1),
    axis.ticks.y = element_line(colour = "black", size = 1),
    axis.title = element_text(size=14, face="bold"),
    strip.text = element_text(size = 20, color = 'black', face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    # Legend: bottom-left inside plot
    legend.position = c(0, 0),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(color = "black", size = 1, fill = "white"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(12, "pt"),
    axis.ticks.length = unit(3, "pt")
  ) +
  # Make legend swatches clearer/bigger
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.6, size = 4)
    )
  )
  