# =========================================================
# GLMM-only soil analysis 
# Plots 10 x 2 m, >=50 m apart within pairs, >=150 m among pairs
# Sampled June–August, 2024
#This analysis supports the LMMs shown in the main text
# =========================================================

# Packages
library(dplyr)
library(readr)
library(lme4)          
library(broom.mixed)   
library(performance)   
library(emmeans)       
library(stringr)
library(purrr)
library(tidyr)
library(knitr)
library(officer)

# ------------ Paths ------------
data_path <- file.path(".", "Data", "guava_soil.csv")
out_dir   <- file.path(".", "Outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------ Load & prep ------------
soil <- read_csv(data_path, show_col_types = FALSE)

# Keep first 16 rows as in your script
soil <- soil %>% filter(row_number() %in% 1:16)

# Harmonize Treatment factor (Uninvaded vs Guava)
soil <- soil %>%
  mutate(
    Treatment = case_when(
      str_detect(tolower(Treatment), "n") ~ "Uninvaded",
      str_detect(tolower(Treatment), "g") ~ "Guava",
      TRUE ~ as.character(Treatment)
    ),
    Treatment = factor(Treatment, levels = c("Uninvaded", "Guava"))
  )

# ----------------- GLMM helper -----------------
# Fits Gamma(log) GLMM for a positive response with (1|Pair_number)
# Adds a tiny eps to avoid zeros (Gamma requires y>0)
fit_gamma_glmm <- function(data, response, eps = 1e-6) {
  dat <- data
  
  if (any(!is.finite(dat[[response]]))) {
    stop(paste("Non-finite values in", response))
  }
  
  # Add eps if zeros/non-positive present
  zero_note <- FALSE
  if (any(dat[[response]] <= 0, na.rm = TRUE)) {
    dat[[response]] <- dat[[response]] + eps
    zero_note <- TRUE
  }
  
  fml <- as.formula(paste(response, "~ Treatment + (1|Pair_number)"))
  mod <- glmer(fml, data = dat, family = Gamma(link = "log"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  # Tidy fixed effects (Wald z) + 95% CI
  td <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    mutate(
      estimate_RR = exp(estimate),
      conf.low_RR = exp(conf.low),
      conf.high_RR = exp(conf.high)
    )
  
  # R² (marginal/conditional)
  r2s <- performance::r2(mod)
  r2_marg <- tryCatch(as.numeric(r2s$R2_marginal), error = function(e) NA_real_)
  r2_cond <- tryCatch(as.numeric(r2s$R2_conditional), error = function(e) NA_real_)
  
  # Emmeans contrast: Guava vs Uninvaded on response scale (ratio)
  emc <- suppressMessages(emmeans::emmeans(mod, ~ Treatment, type = "response"))
  cn  <- suppressMessages(emmeans::contrast(emc, method = "revpairwise"))
  cn_summ <- as.data.frame(summary(cn, infer = c(TRUE, TRUE))) %>%
    mutate(response = response)
  
  list(
    model = mod,
    tidy_table = td,
    r2 = data.frame(R2_marginal = r2_marg, R2_conditional = r2_cond),
    emm_contrast = cn_summ,
    zero_adjusted = zero_note
  )
}

# ----------------- Variables to model -----------------
responses <- c("Ntotal","Organic_matter","P_avail","Carbon","NH4","NO3","CN")

# Fit all models
fits <- purrr::map(responses, ~ fit_gamma_glmm(soil, .x))
names(fits) <- responses

# ----------------- Assemble export tables -----------------

# 1) Fixed-effects tables (long)
fixed_tables <- purrr::imap_dfr(fits, function(obj, resp) {
  obj$tidy_table %>%
    mutate(response = resp) %>%
    select(response, term, estimate, std.error, statistic, p.value,
           conf.low, conf.high, estimate_RR, conf.low_RR, conf.high_RR)
}) %>%
  mutate(
    term = dplyr::case_when(
      term == "(Intercept)" ~ "(Intercept)",
      stringr::str_detect(term, "TreatmentGuava") ~ "Treatment: Guava vs Uninvaded",
      TRUE ~ term
    )
  )

# 2) R2 table
r2_table <- purrr::imap_dfr(fits, function(obj, resp) {
  tibble(
    response = resp,
    R2_marginal = obj$r2$R2_marginal,
    R2_conditional = obj$r2$R2_conditional,
    zero_adjustment_applied = obj$zero_adjusted
  )
})

# 3) Emmeans contrasts (Guava / Uninvaded)
emm_table <- purrr::map_dfr(fits, "emm_contrast") %>%
  dplyr::rename(
    LowerCL = dplyr::any_of(c("lower.CL", "asymp.LCL")),
    UpperCL = dplyr::any_of(c("upper.CL", "asymp.UCL"))
  ) %>%
  dplyr::select(response, contrast, ratio, SE, df, z.ratio, p.value, LowerCL, UpperCL) %>%
  dplyr::rename(
    Treatment_Contrast = contrast,
    Rate_Ratio = ratio
  )
# ----------------- Write CSVs -----------------
write_csv(fixed_tables, "GLMM_fixed_effects_long.v2.csv")

