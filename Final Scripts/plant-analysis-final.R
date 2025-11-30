### Strawberry Guava Project
### Plant community analyses

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
library(vegan)
library(stringr)
library(tibble)
library(cowplot)
library(rlang)

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
# 1) Canopy openness
#    Note: Openness is bounded (0–100). For LMMs,
#    log-transform can stabilize variance if values > 0.
# # ==============================================

# Check transform & normality (on the transformed variable)
if (all(cano_open$Canopy_openness > 0, na.rm = TRUE)) {
  cano_open <- cano_open %>%
    mutate(log_open = log(Canopy_openness))
  
  # Shapiro on transformed raw values (quick heuristic);
  # more important is residual normality after modeling.
  shapiro.test(cano_open$log_open)
}

# LMM on log-transformed openness
mod_open_lmm <- lmer(log_open ~ Treatment + (1 | Pair), data = cano_open)
summary(mod_open_lmm)
anova(mod_open_lmm, ddf = "Kenward-Roger")  # KR F-tests
plot(mod_open_lmm)# Diagnostics for LMM

# ==============================================
# 2) Robel pole
#   Start with LMM on log(x+1)
# ==============================================

# Transform & quick test (heuristic)
canopy <- canopy %>% mutate(log_robel = log(RP_avg + 1))
shapiro.test(canopy$log_robel)

# LMM on log-transformed Robel pole
mod_robel_lmm <- lmer(log_robel ~ Treatment + (1 |Pair), data = canopy)
summary(mod_robel_lmm)
anova(mod_robel_lmm, ddf = "Kenward-Roger")
plot(mod_robel_lmm)

#**********************************************
#Seedling data
#**********************************************

# 1) combine subplot data
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

colnames(seedlings_ren) <- (c("X0", "X1","X2","Inf.")) #X0 = richness, X1 = Shannon-Hill Diversity, X2 Inverse Simpson, Inf = Weighted diversity

seedlings_ren <- merge(seedlings_ren, seedlings_sum, by = 0) %>%
  mutate(seedlings_ren, plotnumber = Row.names) 

#2) Check Normality
shapiro.test(seedlings_ren$Totalstems) 
hist(seedlings_ren$Totalstems, breaks = 4)

shapiro.test(log(seedlings_ren$X0)) 
hist(seedlings_ren$X0, breaks = 6)

shapiro.test(log(seedlings_ren$X1)) 
hist(seedlings_ren$X1)

# 3) LMM to test effects of guava invasion
#####Total stems

abund_totalstems<- lmer(log(Totalstems) ~ Treatment + (1|Pair), data = seedlings_ren)
summary(abund_totalstems)
anova(abund_totalstems, F = "Kenward-Rogers")
plot(abund_totalstems)

######Species richness

richness_stems<- lmer(log(X0) ~ Treatment + (1|Pair), data = seedlings_ren) 
summary(richness_stems)
anova(richness_stems, F = "Kenward-Rogers")
plot(richness_stems)

####Hill-Shannon diversity

shan_lmer <- lmer(log(X1) ~ Treatment + (1|Pair), data = seedlings_ren)
summary(shan_lmer)
anova(shan_lmer, F = "Kenward-Rogers")
plot(shan_lmer)

#***************************************************
# Tree seedling data
#***************************************************

# 1) Organize data
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

#2) check for normality assumptions

#Total stems
shapiro.test(log(tree_ren$Totalstems))
hist(tree_ren$Totalstems, breaks = 4)

#3) test LMMs to determine invasion impacts

####Total stems
tree_abund_lmer <- lmer(log(Totalstems) ~ Treatment + (1|Pair), 
                        data = tree_ren)
summary(tree_abund_lmer) 
anova(tree_abund_lmer, ddf = "Kenward-Roger")
plot(tree_abund_lmer)

#========================================================
#Native tree seedlings
#=======================================================

#1) Organize data
natives <- tree_sum[c(2:22)]

native_ren <- renyi(natives[1:15], scale = c(0,1,2,Inf), hill = TRUE)

colnames(native_ren) <- (c("X0", "X1","X2","Inf.")) 

native_ren <- merge(native_ren, natives, by = 0) %>%
  mutate(native_ren, plotnumber = Row.names) 

#2) Check normality
#Native tree seedling density
shapiro.test(log(native_ren$totalNONguava))
hist(native_ren$Totalstems, breaks = 4)

#Native tree Hill-Shannon diversity
shapiro.test(log(native_ren$X1)) 
hist(native_ren$X1)

#3) Evaluate impacts of guava on native tree seedlings

####Native tree seedling density
native_abund_lmer <- lmer(log(totalNONguava) ~ Treatment + (1|Pair), 
                        data = native_ren)
summary(native_abund_lmer)
anova(native_abund_lmer, ddf = "Kenward-Roger")
plot(native_abund_lmer)

####Native tree seedling Hill-Shannon Diversity
#LMM
native_X1_lmm <- lmer(log(X1) ~ Treatment + (1|Pair), data = native_ren)
summary(native_X1_lmm)
anova(native_X1_lmm, ddf = "Kenward-Roger")
plot(native_X1_lmm)

##plot figures for manuscript
# ---- Shared aesthetics ----
base_theme <- theme_minimal(base_size = 13) +
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
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)
        )

fill_vals <- c("darkolivegreen", "beige")

# ---- Panel A: Total woody seedlings ----
woody_abund <- ggplot(tree_ren, aes(Treatment, Totalstems, fill = Treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = fill_vals) +
  scale_x_discrete(labels = c("No Guava", "Guava")) +
  coord_cartesian(ylim = c(0, 50), expand = FALSE) +
  labs(y = "Tree Seedling Density", x = NULL) +
  # significance mark centered over the two boxes
  annotate("text", x = 1.5, y = Inf, label = "**", size = 6, vjust = 1.2) +
  base_theme

# ---- Panel B: Total native tree seedlings ----
native_abund <- ggplot(native_ren, aes(Treatment, totalNONguava, fill = Treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = fill_vals) +
  scale_x_discrete(labels = c("No Guava", "Guava")) +
  coord_cartesian(ylim = c(0, 15), expand = FALSE) +
  labs(y = "Native Tree Seedling Density", x = NULL) +
  annotate("text", x = 1.5, y = Inf, label = "NS", size = 5, vjust = 1.2) +
  base_theme

# ---- Panel C: Hill-Shannon of native tree seedlings ----
native_X1 <- ggplot(native_ren, aes(Treatment, X1, fill = Treatment)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.4), pch=21, size=5, alpha=0.6, show.legend = FALSE)+ 
  scale_fill_manual(values = fill_vals) +
  scale_x_discrete(labels = c("No Guava", "Guava")) +
  coord_cartesian(ylim = c(0, 5), expand = FALSE) +
  labs(y = "Native Tree Seedling Diversity", x = NULL) +
  annotate("text", x = 1.5, y = Inf, label = "NS", size = 5, vjust = 1.2) +
  base_theme

# ---- Arrange with consistent top-left panel letters ----
# label_x ~0.02 keeps the tag near the left edge; label_y ~0.98 near the top edge
cowplot::plot_grid(
  woody_abund, native_abund, native_X1,
  labels = c("A", "B", "C"),
  ncol = 2,
  label_size = 16,
  label_fontface = "bold",
  label_x = 0.17,  # left placement
  label_y = 0.95,  # top placement
  hjust = 0, vjust = 1
)

#*****************************************
#Permanova for native tree seedlings
#*****************************************

# ---- 1) Load data ----
#relative pathname
native <- file.path(".", "Data","native_ren.csv")

#import data
df <- read.csv(native, check.names = FALSE, stringsAsFactors = FALSE)

# Invasion status factor order (controls legend & palette mapping)
df$Invasion_status <- factor(df$Invasion_status, levels = c("No Guava", "Guava"))

# Pair as a factor for blocking in PERMANOVA
if ("Pair" %in% names(df)) df$Pair <- factor(df$Pair)

# ---- 2) Build species matrix (exclude numeric metadata like Pair) ----
species_mat <- df %>%
  select(where(is.numeric)) %>%
  select(-any_of(c("Pair"))) %>%
  as.data.frame()

# Drop rows with zero total abundance (avoids ordination issues)
keep_rows <- rowSums(species_mat, na.rm = TRUE) > 0
df2 <- df[keep_rows, , drop = FALSE]
species_mat <- species_mat[keep_rows, , drop = FALSE]

# ---- 3) Bray–Curtis distance + PERMANOVA (by invasion status) ----
set.seed(123)
bray <- vegdist(species_mat, method = "bray")

# Block permutations within Pair if present
perm <- adonis2(bray ~ Invasion_status, permutations = 999, data = df2)
print(perm)

#---- Homogeneity of multivariate dispersions (diagnostic) ----
bd <- betadisper(bray, df2$Invasion_status)
cat("\nDispersion ANOVA:\n"); print(anova(bd))
cat("\nDispersion permutation test:\n"); print(permutest(bd, permutations = 999))

# ---- NMDS (2D) ----
set.seed(123)
nmds <- metaMDS(species_mat, distance = "bray", k = 2,
                trymax = 200, autotransform = FALSE)
cat(sprintf("\nNMDS stress: %.3f\n", nmds$stress))

scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$Invasion_status <- df2$Invasion_status

# ---- Plot: NMDS with treatment ellipses ----
palette_vals <- c("No Guava" = "darkolivegreen", "Guava" = "bisque2")

# Ensure factor order matches the legend
scores_df$Invasion_status <- factor(scores_df$Invasion_status,
                                    levels = c("No Guava","Guava"))

p_nmds <- ggplot(scores_df,
                 aes(NMDS1, NMDS2,
                     color = Invasion_status, fill = Invasion_status)) +
  # filled ellipse (95% normal), then outline
  stat_ellipse(type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.30, linewidth = 0) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8) +
  # points
  geom_point(shape = 21, size = 3.2, stroke = 1.0, alpha = 0.9) +
  coord_equal(xlim = c(-5.5, 5), ylim = c(-2, 2.5)) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  scale_y_continuous(breaks = seq(-2, 3, by = 1)) +
  # hide color legend; show fill legend as "Treatment"
  scale_color_manual(values = c("No Guava" = "black", "Guava" = "black"),
                     guide = "none") +
  scale_fill_manual(values = palette_vals, name = "Notations") +
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
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),  # border thickness
    # Legend: bottom-left inside plot
    legend.position = c(0, 0),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(color = "black", size = 1, fill = "white"), # match panel border
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

