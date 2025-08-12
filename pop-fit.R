
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rstan)
library(bayesplot)
library(ggforce)
library(gridExtra)
library(here)
library(loo)
library(purrr)
library(Cairo)
library(posterior)
library(stringr)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())

theme_set(theme_bw())
theme_update(
  text = element_text(size = 8, family = "Foundry Sterling"), # Font
  plot.title = element_text(hjust = 0),      # Centre-align title
  plot.subtitle = element_text(hjust = 0),   # Centre-align subtitle
  legend.title = element_text(colour = 'black',
                              size = 8, 
                              face = 'bold'),# Legend title
  legend.position = 'right',                 # Move legend
  legend.background = element_blank(),       # Remove legend background
  legend.box.background = element_blank(),   # Remove lengend-box background
  legend.spacing.y = unit(0.01, 'mm'),       # Make legend closer
  legend.key.height = unit(0.4, "cm"),       # Make legend closer
  # panel.grid.minor = element_blank(),      # Remove minor lines
  panel.grid.minor.x = element_blank(),      # Remove minor lines on the x axis
  axis.title.x = element_text(hjust = 1),    # Move title for x-axis
  axis.ticks = element_blank(),              # Remove axis ticks
  aspect.ratio = 1,
  axis.title.y = element_text(hjust = 0.5)   # Move title for y-axis
)

# Population models -------------------------------------------------------
# Input data for Stan
df_clean <- read.csv(here('data', 'patient_data.csv'))

df_clean |> 
 group_by(pid_ori, dayofobs_pq) |> 
 summarise(n = n()) |> 
 arrange(desc(n))

df <- df_clean |> 
  arrange(pid, dayofobs_pq) |> 
  mutate(methb = methb / 100,
         methb7 = methb7 / 100,
         hb = hb / 100)

dose_group_patient <- df |> select(pid_ori, dose_group) |> distinct()
dose_group_vector <- dose_group_patient$dose_group

## Overall data
n_patients <- df$pid |> unique() |> length()

## Patient-specific data
ndays <- df |> group_by(pid) |> summarise(ndays = n()) |> pull(ndays)
nobs <- df |> group_by(pid) |> summarise(nobs = sum(!is.na(methb))) |> pull(nobs)
patient_id <- df |> filter(!is.na(methb)) |> pull(pid)
dose_group <- df |> filter(!is.na(methb)) |> pull(dose_group)

## Time-series data
dose <- df |> pull(pqmgkgday)
ind_obs <- which(!is.na(df$methb))
methb <- df$methb[!is.na(df$methb)]
sex <- df |> group_by(pid) |> summarise(sex = unique(sex)) |> pull(sex)
age_cat <- df |> group_by(pid) |> summarise(age_cat = unique(age_cat)) |> pull(age_cat)
site <- df |> group_by(pid) |> summarise(site = unique(studysite)) |> pull(site) |> as.factor() |> as.integer()
n_sites <- length(unique(df$studysite))
n_age_cats <- length(unique(df$age_cat))

# Make a data list for Stan
data_list <- list(
  n_patients = n_patients,
  ndays = ndays,
  nobs = nobs,
  methb = methb,
  dose = dose,
  dose_group = dose_group_vector,
  ind_obs = ind_obs,
  patient_id = patient_id,
  sex = sex,
  age_cat = age_cat,
  site = site,
  n_sites = n_sites,
  n_age_cats = n_age_cats
)

str(data_list)  # Check structure
sum(ndays) == length(dose)  # Must be TRUE
sum(nobs) == length(methb)  # Must be TRUE

# Stan data
saveRDS(data_list, file = 'data/model_data.rds')

# Fit Stan model ----------------------------------------------------------

n_patients <- data_list$n_patients
model <- stan_model(here('scripts',
                         'pop_2level_base.stan')) # Stan model to fit

# Initial parameter values
init_fun = function(chain_n = 4) {
  my_inits = list()
  for(cc in 1:chain_n){
    my_inits[[cc]] = list(logit_total_rate = -4,
                          logit_mu_alpha = 3,
                          sigma = 0.01, # base: 0.01
                          log_alpha_p = rep(0, n_patients),
                          sigma_alpha = 0.1
                          # site_re_raw = rep(0, n_sites),
                          # sigma_site = 0.1
                          # beta_age = rep(0, n_age_cats),
                          # beta_sex = 0
                          )
  }
  names(my_inits)= paste('chain',1:chain_n, sep = '_')
  my_inits
}

iter <- 2000
fit_multilevel <- sampling(model,
                           data = data_list,
                           iter = iter,
                           chains = 4,
                           warmup = iter/2,
                           seed = 13,
                           init = init_fun(4),
                           thin = 2,
                           control = list(adapt_delta = 0.85))

saveRDS(fit_multilevel, file = here('local-output', 'fit_base.rds'))


# Diagnostics and posterior -----------------------------------------------

fit_multilevel <- read_rds(here('local-output', 'fit_base.rds'))

parsm = c('logit_total_rate',
          'logit_mu_alpha',
          'sigma_alpha',
          'sigma'
)

parsm_ = c('logit_total_rate',
           'logit_mu_alpha',
           'sigma_alpha',
           'sigma'
)

print(fit_multilevel, pars = parsm, digits = 3)
stan_dens(fit_multilevel, pars = parsm, separate_chains = TRUE)
traceplot(fit_multilevel, pars = parsm)
mcmc_neff(neff_ratio(fit_multilevel, pars = parsm), size = 3) + yaxis_text(hjust = 0)
mcmc_acf(fit_multilevel, pars = parsm_, lags = 10)
mcmc_pairs(fit_multilevel,
           pars = parsm_,
           off_diag_args = list(size = 1.5))

posterior <- extract(fit_multilevel) 
mu_scalar <- plogis(median(posterior$logit_mu_alpha))
log_alpha_p_medians <- apply(posterior$log_alpha_p, 2, median)
alpha_medians <- mu_scalar * exp(log_alpha_p_medians)

set.seed(13)
prior_draws <- data.frame(
  logit_total_rate = rnorm(1000, mean = -4, sd = 1),
  logit_mu_alpha   = rnorm(1000, mean = 4, sd = 1),
  sigma_alpha      = truncnorm::rtruncnorm(1000, a = 0, mean = 0.1, sd = 0.05),
  sigma            = truncnorm::rtruncnorm(1000, a = 0, mean = 0.3, sd = 0.1)
)

posterior_draws <- as_draws_df(fit_multilevel, variable = parsm)
posterior_long <- posterior_draws |> 
  select(logit_total_rate, logit_mu_alpha,
         sigma_alpha,
         sigma) |> 
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |> 
  mutate(source = "Posterior")
prior_long <- prior_draws |> 
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |> 
  mutate(source = "Prior")

combined <- bind_rows(posterior_long, prior_long)

xlims <- list(
  logit_total_rate = c(-2.8, -1.8),
  logit_mu_alpha   = c(-4, 9),
  sigma_alpha      = c(0, 5),
  sigma            = c(0, 0.02)
)
filtered_combined <- combined |> 
  rowwise() |>
  filter(
    value >= xlims[[parameter]][1],
    value <= xlims[[parameter]][2]
  ) |> 
  ungroup()

# Plot
ggplot(filtered_combined, aes(x = value, fill = source, colour = source)) +
  # geom_histogram(alpha = 0.4, position = 'dodge') +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("Prior" = "grey70", "Posterior" = "aquamarine3")) +
  scale_colour_manual(values = c("Prior" = "grey70", "Posterior" = "aquamarine3")) +
  labs(x = '\nParameter value',
       y = 'Density\n',
       fill = '',
       colour = '')

#  Model comparison --------------------------------------------------------

# log_lik_intercept <- extract_log_lik(fit_multilevel, merge_chains = FALSE)
# log_lik_extended <- extract_log_lik(fit_multilevel_ext, merge_chains = FALSE)
# # log_lik_site <- extract_log_lik(fit_multilevel, merge_chains = FALSE)
# # # # log_lik_covariates <- extract_log_lik(fit_multilevel, merge_chains = FALSE)
# # # log_lik_intercept_multip <- extract_log_lik(fit_multilevel, merge_chains = FALSE)
# # # log_like_site_multip <- extract_log_lik(fit_multilevel, merge_chains = FALSE)
# # # #
# loo_intercept <- loo(log_lik_intercept); loo_intercept
# loo_extended <- loo(log_lik_extended); loo_extended
# # loo_site <- loo(log_lik_site); loo_site
# # # # loo_covariates <- loo(log_lik_covariates); loo_covariates
# # # loo_intercept_multip <- loo(log_lik_intercept_multip); loo_intercept_multip
# # # loo_site_multip <- loo(log_like_site_multip); loo_site_multip
# # # #
# # # # loo_compare(list('intercept' = loo_intercept,
# # # #                  'site' = loo_site,
# # # #                  'covariates' = loo_covariates
# # # #                  ))
# # #
# # loo_compare(list('intercept_multip' = loo_intercept_multip,
# #                  'site_multip'= loo_site_multip,
# #                  'intercept_base' = loo_intercept,
# #                  'site' = loo_site
# # ))
# 
# loo_compare(list('intercept_base' = loo_intercept,
#                  'extended' = loo_extended
# ))

# Posterior predictive checks ---------------------------------------------

# Individual
posterior_samples <- as_draws_df(fit_multilevel)

pred_summary <- posterior_samples %>%
  summarise_draws(~ quantile(.x, probs = c(0.025, 0.5, 0.975)))

# Extract day index from `variable`
pred_summary_temp <- pred_summary %>%
  filter(grepl("pred_methb", variable)) %>%
  mutate(index = as.numeric(gsub("pred_methb\\[|\\]", "", variable)))  # Extract index

patient_mapping <- tibble(
  index = 1:sum(ndays),
  patient_id = rep(1:n_patients, times = ndays),
  day = unlist(lapply(ndays, seq_len))  # Generates day sequences per patient
)

obs_data <- tibble(
  index = ind_obs,  # Observed indices in pred_methb
  methb = methb
) %>%
  left_join(patient_mapping, by = "index")

pred_summary_fin <- pred_summary_temp %>%
  mutate(index = as.numeric(str_extract(variable, "\\d+"))) %>%
  left_join(patient_mapping, by = "index") |>
  filter(day <= 30) |> 
  select(patient_id, day, `2.5%`, `50%`, `97.5%`) |> 
  rename(lower = `2.5%`,
         median = `50%`,
         upper = `97.5%`)

plot_data <- pred_summary_fin %>%
  left_join(obs_data, by = c("patient_id", "day")) |> 
  mutate(day = day - 1)

patient_id_factor <- factor(df$pid_ori)
patient_index <- as.integer(patient_id_factor)

# Original patient IDs
patient_lookup <- tibble(
  patient_index = seq_along(levels(patient_id_factor)),
  patient_id = levels(patient_id_factor)
)

plot_data <- pred_summary_fin %>%
  left_join(obs_data, by = c("patient_id", "day")) |> 
  mutate(day = day - 1) |> 
  left_join(patient_lookup, by = c("patient_id" = "patient_index")) |> 
  rename(pid = patient_id.y)

# Count number of unique patients
n_patients <- plot_data |> distinct(patient_id) |> nrow()
n_per_page <- 20
n_pages <- ceiling(n_patients / n_per_page)

# Create the PDF
CairoPDF("methb_plots_base.pdf", width = 8, height = 10)

for (i in 30:31) {
  p <- plot_data |> 
    ggplot() +
    geom_point(aes(day, methb * 100), size = 0.2) +
    geom_ribbon(aes(x = day, ymin = lower * 100, ymax = upper * 100),
                alpha = 0.5, fill = 'aquamarine3') +
    geom_line(aes(day, median * 100), linewidth = 0.3, colour = 'aquamarine3') +
    facet_wrap_paginate(~pid, ncol = 4, nrow = 5, page = i) +
    scale_y_continuous(breaks = seq(0, 50, by = 5),
                       labels = seq(0, 50, by = 5),
                       limits = c(0, 35)) +
    scale_x_continuous(
      limits = c(0, 28),
      labels = seq(0, 28, by = 7),
      breaks = seq(0, 28, by = 7),
      expand = expansion(mult = c(0.05, 0.05))) +
    theme(strip.text = element_text(size = 5),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          legend.key.width = unit(0.2, 'cm')) +
    scale_colour_discrete(na.translate = FALSE) +
    labs(x = '\nDays elapsed since primaquine administration at baseline', 
         y = 'Methaemoglobin\n')
  
  print(p)  # Add this page to the PDF
}

dev.off()
  
# Residuals ---------------------------------------------------------------

# By regimen duration
resid_mat <- extract(fit_multilevel, 'residuals')$residuals
median_resid <- apply(resid_mat, 2, median)
residual_df <- df[!is.na(df$methb), ] |> 
  mutate(median_residual = median_resid)
width_df <- residual_df |>
  count(pqdur_exp, dayofobs_pq) |>
  mutate(box_width = sqrt(n) / max(sqrt(n)))
residual_df2 <- residual_df |>
  left_join(width_df, by = c("pqdur_exp", "dayofobs_pq")) |> 
  mutate(sid = str_sub(pid_ori, 1, 5),
         sid = case_when(sid == 'BGGDZ' ~ 'Pasaribu 2013',
                         sid == 'BSAJK' ~ 'Llanos-Cuentas 2019',
                         sid == 'EOGVV' ~ 'Taylor 2019',
                         sid == 'FUCCA' ~ 'Nelwan 2015',
                         sid == 'KTMTV' ~ 'Lacerda 2019',
                         sid == 'PRSXC' ~ 'Llanos-Cuentas 2014',
                         sid == 'UJIUX' ~ 'Chu 2019',
                         sid == 'RCTFJ' ~ 'Sutanto 2013',
                         is.na(sid) ~ NA_character_,
                         TRUE ~ 'Check me!') |> 
           factor(levels = c('Pasaribu 2013',
                             'Sutanto 2013',
                             'Llanos-Cuentas 2014',
                             'Nelwan 2015',
                             'Chu 2019',
                             'Lacerda 2019',
                             'Llanos-Cuentas 2019',
                             'Taylor 2019'),
                  labels = c('Pasaribu 2013',
                             'Sutanto 2013',
                             'Llanos-Cuentas 2014',
                             'Nelwan 2015',
                             'Chu 2019',
                             'Lacerda 2019',
                             'Llanos-Cuentas 2019',
                             'Taylor 2019')))

ggplot(data = residual_df2) +
  geom_boxplot(aes(dayofobs_pq, median_residual,
                   group = dayofobs_pq),
               outlier.size = 0.5, outlier.alpha = 0.2,
               size = 0.25, varwidth = T) +
  facet_wrap(~pqdur_exp) +
  theme(strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 28, by = 7)) +
  labs(y = 'Residual\n',
       x = '\nDays elapsed since primaquine initiation at baseline')

ggplot(data = residual_df2) +
  geom_boxplot(aes(dayofobs_pq, median_residual * 100,
                   group = dayofobs_pq),
               outlier.size = 0.5, outlier.alpha = 0.2,
               size = 0.25, varwidth = T) +
  facet_grid(cols = vars(sid), rows = vars(pqdur_exp)) +
  theme(strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 28, by = 7),
                     limits = c(-0.5, 28)) +
  scale_y_continuous(limits = c(-0.1, 0.1) * 100) +
  labs(y = 'Residual\n',
       x = '\nDays elapsed since primaquine initiation at baseline')

# Residual vs dose
residual_df2 |> 
  select(sid, dayofobs_pq, pqmgkgday, median_residual, pqdur_exp) |> 
  filter(dayofobs_pq %in% c(6, 7, 13, 14)) |> 
ggplot() +
  geom_point(aes(x = pqmgkgday,
                 y = median_residual * 100),
             alpha = 0.1, size = 1) +
  labs(x = '\nPrimaquine daily dose (mg/kg)',
       y = 'Residual\n') +
  facet_grid(rows = vars(dayofobs_pq),
             cols = vars(pqdur_exp)) +
  # geom_smooth(aes(x = pqmgkgday, y = median_residual),
  #             method = "gam", formula = y ~ s(x),
  #             se = T, colour = "firebrick") +
  theme(
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  ) +
  scale_x_continuous(labels = seq(0, 1.5, by = 0.25),
                     breaks = seq(0, 1.5, by = 0.25)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    name = "Residual\n",  # Left y-axis label
    sec.axis = sec_axis(~., name = "Day of measurement", labels = NULL)  # Right y-axis label only
  )

posterior_plot <- plot_data |> 
  mutate(pid_ori = pid,
         sid = str_sub(pid_ori, 1, 5),
         sid = case_when(sid == 'BGGDZ' ~ 'Pasaribu 2013',
                         sid == 'BSAJK' ~ 'Llanos-Cuentas 2019',
                         sid == 'EOGVV' ~ 'Taylor 2019',
                         sid == 'FUCCA' ~ 'Nelwan 2015',
                         sid == 'KTMTV' ~ 'Lacerda 2019',
                         sid == 'PRSXC' ~ 'Llanos-Cuentas 2014',
                         sid == 'UJIUX' ~ 'Chu 2019',
                         sid == 'RCTFJ' ~ 'Sutanto 2013',
                         is.na(sid) ~ NA_character_,
                         TRUE ~ 'Check me!') |> 
           factor(levels = c('Pasaribu 2013',
                             'Sutanto 2013',
                             'Llanos-Cuentas 2014',
                             'Nelwan 2015',
                             'Chu 2019',
                             'Lacerda 2019',
                             'Llanos-Cuentas 2019',
                             'Taylor 2019'),
                  labels = c('Pasaribu 2013',
                             'Sutanto 2013',
                             'Llanos-Cuentas 2014',
                             'Nelwan 2015',
                             'Chu 2019',
                             'Lacerda 2019',
                             'Llanos-Cuentas 2019',
                             'Taylor 2019'))) |> 
  left_join(df |> 
              select(pqdur_exp, pid_ori) |> 
              distinct(),
            by = "pid_ori")

ggplot(data = posterior_plot) +
  geom_point(aes(x = day, y = methb * 100),
             alpha = 0.05, size = 0.5) +
  # geom_line(aes(day, median, group = patient_id),
  #           colour = 'aquamarine3', alpha = 0.075) +
  facet_grid(cols = vars(sid), rows = vars(pqdur_exp)) +
  theme(strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 28, by = 7),
                     limits = c(0, 28)) +
  labs(y = 'Methaemoglobin (%)\n',
       x = '\nDays elapsed since primaquine initiation at baseline')

# Alpha vs dose
# Duration
df |>
  arrange(pid, dayofobs_pq) |>
  group_by(pid) |>
  slice(1) |>
  ungroup() |> 
  mutate(median_alpha = alpha_medians,
         median_log_alpha_p = log_alpha_p_medians) |> 
ggplot() +
  geom_point(aes(x = pqmgkgday, y = log_alpha_p_medians,
                 # colour = pqdur_exp
                 ),
             alpha = 0.2, size = 2) +
  # labs(x = '\nPrimaquine daily dose (mg/kg)',
  #      y = 'Expected dose effect\n',
  #      colour = 'Primaquine duration') +
  scale_x_continuous(labels = seq(0, 1.5, by = 0.25),
                     breaks = seq(0, 1.5, by = 0.25)) +
  # scale_y_log10(
  #   breaks = c(0.5, 1, 2, 4, 8),  # adjust as needed
  #   labels = scales::label_number()
  # ) +
  # geom_smooth(aes(x = pqmgkgday, y = log_alpha_p_medians),
  #             method = "gam", formula = y ~ s(x),
  #             se = T, colour = "gray60", linewidth = 0.25) +
  theme(legend.position.inside = c(0.75, 0.85))


# Simulation-based visual predictive checks, by dose ----------------------

dose_cat_df <- df |> 
  filter(dayofobs_pq == 0) |> 
  select(pid, pqmgkgday) |> 
  mutate(dose_cat = case_when(
    pqmgkgday < 0.375 ~ "Low",
    pqmgkgday < 0.75 ~ "Intermediate",
    TRUE ~ "High"
  )) |> 
  mutate(dose_cat = factor(dose_cat, level = c("Low",
                                               "Intermediate",
                                               "High")))

posterior_slice <- as_draws_df(fit_multilevel) |> 
  slice_sample(n = 500) |> 
  mutate(.draw = row_number())

simulated_draws <- posterior_slice |> 
  select(.draw, starts_with("pred_methb"), sigma) |> 
  pivot_longer(cols = starts_with("pred_methb"), 
               names_to = "param",
               values_to = "pred_methb") |> 
  mutate(obs_index = as.integer(str_extract(param, "\\d+"))) |> 
  left_join(df |> mutate(obs_index = row_number()),
            by = "obs_index") |>
  left_join(dose_cat_df, by = "pid") |> 
  mutate(methb_sim = rnorm(n(), mean = pred_methb, sd = sigma)) |> 
  rename(pqmgkgday = pqmgkgday.x) |> 
  select(-pqmgkgday.y)

sim_percentiles <- simulated_draws |> 
  group_by(.draw, dayofobs_pq, dose_cat) |>
  mutate(methb_sim = if_else(methb_sim < 0, 0, methb_sim)) |> 
  summarise(p5 = quantile(methb_sim, 0.05),
            p50 = quantile(methb_sim, 0.5),
            p95 = quantile(methb_sim, 0.95),
            .groups = "drop")

ppc_summary <- sim_percentiles |> 
  pivot_longer(cols = c(p5, p50, p95),
               names_to = 'percentile',
               values_to = 'value') |> 
  group_by(dayofobs_pq, dose_cat, percentile) |> 
  summarise(
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    median = median(value),
    .groups = "drop"
  )

n_to_summarise <- 30
obs_summary <- df %>%
  left_join(dose_cat_df, by = "pid") %>%
  group_by(dayofobs_pq, dose_cat) %>%
  summarise(
    n_obs = sum(!is.na(methb)),
    obs_p5 = if (n_obs >= n_to_summarise) quantile(methb, 0.05, na.rm = TRUE) else NA_real_,
    obs_p50 = if (n_obs >= n_to_summarise) quantile(methb, 0.5, na.rm = TRUE) else NA_real_,
    obs_p95 = if (n_obs >= n_to_summarise) quantile(methb, 0.95, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  select(-n_obs) |> 
  drop_na(obs_p5)

obs_points <- df %>%
  left_join(dose_cat_df, by = "pid")

ggplot() +
  geom_point(data = obs_points,
             aes(x = dayofobs_pq,
                 y = methb * 100),
             size = 0.7, alpha = 0.015, colour = "black",
             shape = 19) +
  geom_ribbon(data = ppc_summary,
              aes(x = dayofobs_pq,
                  ymin = lower * 100,
                  ymax = upper * 100,
                  group = percentile),
              alpha = 0.7, fill = 'aquamarine3') +
  facet_grid(cols = vars(dose_cat)) +
  # geom_line(data = ppc_summary,
  #           aes(x = dayofobs_pq, y = median, group = percentile),
  #           linewidth = 0.25, colour = 'grey50') +
  geom_line(data = obs_summary, aes(x = dayofobs_pq,
                                    y = obs_p50 * 100),
            linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(data = obs_summary, aes(x = dayofobs_pq,
                                    y = obs_p5 * 100),
            linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(data = obs_summary, aes(x = dayofobs_pq,
                                    y = obs_p95 * 100),
            linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  labs(
    x = "\nDays elapsed since primaquine initiation",
    y = "Methaemoglobin (%)\n"
  ) +
  theme(strip.text = element_text(face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1) * 5, "mm")) +
  # scale_y_continuous(limits = c(0, 0.225)) +
  scale_x_continuous(limits = c(0, 28),
                     expand = c(0, 0),
                     breaks = seq(0, 28, by = 7)) +
  scale_y_continuous(limits = c(0, 25)) +
  theme(legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(1, "lines"),
        text = element_text(family = "Foundry Sterling", size = 10))

# ...
# Prepare the dataset
df <- df |>
  group_by(pid) |>
  mutate(pqmgkgday_day0 = pqmgkgday[dayofobs_pq == 0][1]) |>
  ungroup() |> 
  mutate(obs_index = row_number())

# Extract posterior draws
posterior_slice <- as_draws_df(fit_multilevel) |>
  slice_sample(n = 500) |>
  mutate(.draw = row_number())

# Get simulated predictions on day 7 only
sim_day7 <- posterior_slice |>
  select(.draw, starts_with("pred_methb"), sigma) |>
  pivot_longer(cols = starts_with("pred_methb"),
               names_to = "param", values_to = "pred_methb") |>
  mutate(obs_index = as.integer(str_extract(param, "\\d+"))) |>
  left_join(df, by = "obs_index") |>
  filter(!is.na(methb7), dayofobs_pq == 7) |>  # restrict to day 7 observations
  mutate(methb_sim = rnorm(n(), mean = pred_methb, sd = sigma))

# Smooth by pqmgkgday rounded to nearest 0.05
sim_day7 <- sim_day7 |>
  mutate(pqmgkgday_day0_smooth = round(pqmgkgday_day0 / 0.225) * 0.225)

# Summarise posterior predictive simulations
ppc_day7_summary <- sim_day7 |>
  group_by(.draw, pqmgkgday_day0_smooth) |>
  summarise(p50 = quantile(methb_sim, 0.5), .groups = "drop") |>
  group_by(pqmgkgday_day0_smooth) |>
  summarise(
    lower = quantile(p50, 0.025),
    upper = quantile(p50, 0.975),
    median = median(p50),
    .groups = "drop"
  )

# Summarise observed data
obs_day7_summary <- df |>
  filter(!is.na(methb7), dayofobs_pq == 7) |>
  mutate(pqmgkgday_day0_smooth = round(pqmgkgday_day0 / 0.225) * 0.225) |>
  group_by(pqmgkgday_day0_smooth) |>
  summarise(
    n_obs = n(),
    obs_p50 = if (n_obs >= 25) quantile(methb7, 0.5) else NA_real_,
    .groups = "drop"
  ) |>
  drop_na(obs_p50)

# Plot
ggplot() +
  geom_point(data = df |> filter(!is.na(methb7), dayofobs_pq == 7),
             aes(x = pqmgkgday_day0,
                 y = methb7 * 100),
             size = 1, alpha = 0.1, colour = "black", shape = 19) +
  geom_ribbon(data = ppc_day7_summary,
              aes(x = pqmgkgday_day0_smooth,
                  ymin = lower * 100,
                  ymax = upper * 100),
              alpha = 0.6, fill = 'aquamarine3') +
  # geom_line(data = ppc_day7_summary,
  #           aes(x = pqmgkgday_day0_smooth, y = median),
  #           colour = "aquamarine4", linewidth = 0.7) +
  geom_line(data = obs_day7_summary,
            aes(x = pqmgkgday_day0_smooth, y = obs_p50 * 100),
            linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  labs(
    x = "\nPrimaquine daily dose (mg/kg)",
    y = "Day 7 methaemoglobin (%)\n"
  ) +
  scale_x_continuous(limits = c(0.1, 1.4),
                     breaks = seq(0.25, 1.5, by = 0.25)) +
  scale_y_continuous(limits = c(0, 25)) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 10, family = "Foundry Sterling")
  )
  # scale_y_continuous(limits = c(0, 0.225)) 

# Surrogate metrics -------------------------------------------------------

# Cmax
predicted_peak_per_patient <- pred_summary_fin |> 
  group_by(patient_id) |> 
  slice_max(order_by = median, n = 1, with_ties = F) |> 
  select(patient_id, cmax = median)

# C3
predicted_c3_per_patient <- pred_summary_fin |> 
  filter(day == 3) |> 
  select(patient_id, c3 = median)

# C7
predicted_c7_per_patient <- pred_summary_fin |> 
  filter(day == 7) |> 
  select(patient_id, c7 = median)

# AUC
## 0 to 14
predicted_auc014_per_patient <- pred_summary_fin %>%
  filter(day >= 0, day <= 14) %>%
  arrange(patient_id, day) %>%
  group_by(patient_id) %>%
  summarise(
    auc014 = sum(diff(day) * (head(median, -1) + tail(median, -1)) / 2),
    .groups = "drop"
  )

## 0 to 28
predicted_auc028_per_patient <- pred_summary_fin %>%
  filter(day >= 0, day <= 28) %>%
  arrange(patient_id, day) %>%
  group_by(patient_id) %>%
  summarise(
    auc028 = sum(diff(day) * (head(median, -1) + tail(median, -1)) / 2),
    .groups = "drop"
  )

pid_ori_to_join <- df_clean |>
  distinct(pid, .keep_all = TRUE) |>
  select(pid, pid_ori) |> 
  rename(patient_id = pid)

summary_metrics <- predicted_peak_per_patient |> 
  left_join(predicted_c7_per_patient, by = c('patient_id')) |> 
  left_join(predicted_c3_per_patient, by = c('patient_id')) |> 
  left_join(predicted_auc014_per_patient, by = c('patient_id')) |> 
  left_join(predicted_auc028_per_patient, by = c('patient_id')) |> 
  left_join(pid_ori_to_join, by = c('patient_id')) |> 
  ungroup() |> 
  select(-patient_id) |> 
  select(pid_ori, everything())

write_rds(summary_metrics, file = here('data', 'summary_metrics.rds'))











