
library(tidyverse)
library(here)
library(coxme)
library(lme4)
library(broom)
library(performance)
library(broom.mixed)
library(rms)
library(brms)
library(Surrogate)
rstan::rstan_options (auto_write = TRUE)
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

main_methb7 <- read_rds(here('data', 'main_methb.rds')) |> 
  select(sid_author,
         studysite,
         pid, 
         pqday_cat,
         pqdur_exp,
         dlast120_pq,
         outcome7to120_pq,
         age,
         sex,
         log_pvdens,
         pqmgkgday,
         pqmgkgtot,
         trt1,
         methb7) |> 
  mutate(methb7 = methb7 / 100,
         methb7 = if_else(methb7 == 0, 0.01, methb7))
summary_metrics <- read_rds(here('data', 'summary_metrics.rds')) |> 
  rename(pid = pid_ori)

surrogacy_data <- main_methb7 |> # 1770 to 1747 due to NAs in parasite density
  left_join(summary_metrics, by = c('pid'))

surrogacy_data_7d <- surrogacy_data |> 
  filter(pqdur_exp == '7 days')
surrogacy_data_14d <- surrogacy_data |> 
  filter(pqdur_exp == '14 days')

id_4m <- read_rds(here('data', 'id_4m.rds'))
surrogacy_4m <- filter(surrogacy_data, pid %in% id_4m)

# Plots -------------------------------------------------------------------

surrogacy_plot |> 
  ggplot() +
  geom_point(aes(y = methb7 * 100, x = c7 * 100),
             alpha = 0.15, size = 1.75) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dotted', colour = 'firebrick3', linewidth = 0.75) +
  geom_smooth(aes(x = c7 * 100, y = methb7 * 100),
              method = "loess", span = 1, linewidth = 0.5,
              colour = "aquamarine3", size = 1) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 0.25) * 100) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.25) * 100) +
  labs(x = '\nPredicted day 7 methaemoglobin (%)',
       y = 'Observed day 7 methaemoglobin (%)\n')

bland_altman_df <- surrogacy_plot |> 
  mutate(average = 100 * (methb7 + c7) / 2,
         difference = 100 * (methb7 - c7))
mean_diff <- mean(bland_altman_df$difference, na.rm = T)
sd_diff <- sd(bland_altman_df$difference, na.rm = T)
loa_upper <- mean_diff + 1.96 * sd_diff
loa_lower <- mean_diff - 1.96 * sd_diff

ggplot(bland_altman_df, aes(x = average, y = difference)) +
  geom_point(alpha = 0.15, size = 1.75) +
  geom_hline(yintercept = mean_diff, colour = 'aquamarine3',
             linetype = 'dashed', linewidth = 1) +
  geom_hline(yintercept = loa_upper, colour = 'firebrick3',
             linetype = 'dotted', linewidth = 0.75) +
  geom_hline(yintercept = loa_lower, colour = 'firebrick3',
             linetype = 'dotted', linewidth = 0.75) +
  labs(x = '\nAverage of observed and predicted day 7 methaemoglobin (%)',
       y = 'Observed - predicted day 7 methaemoglobin (%)\n') +
  scale_y_continuous(breaks = seq(-10, 10, by = 2))

surrogacy_plot <- surrogacy_4m |> 
  mutate(outcome7to120_pq = if_else(outcome7to120_pq == 0,
                                    'No recurrence',
                                    'Recurrence'),
         pqday_cat = case_when(pqday_cat == 'Low <0.375' ~ 'Low\n[0, 0.375)',
                               pqday_cat == 'Intermediate >=0.375 and <0.75' ~ 'Intermediate\n[0.375, 0.75)',
                               pqday_cat == 'High >=0.75' ~ 'High\n[0.75, ∞)',
                               is.na(pqday_cat) ~ NA_character_,
                               TRUE ~ NA_character_),
         pqday_cat = factor(pqday_cat,
                            levels = c('Low\n[0, 0.375)',
                                       'Intermediate\n[0.375, 0.75)',
                                       'High\n[0.75, ∞)'))) |> 
  select(pqday_cat, auc028, cmax, c7, methb7, outcome7to120_pq)

surr_long <- surrogacy_plot |> 
  pivot_longer(
    cols = c(methb7, c7, cmax, auc028),
    names_to = 'metric',
    values_to = 'value'
    ) |> 
  mutate(metric = case_when(
    metric == 'methb7'  ~ 'C[day-7]*plain(", empirical")',
    metric == 'c7'      ~ 'C[day-7]*plain(", predicted")',
    metric == 'cmax'    ~ 'C[max]*plain(", predicted")',
    metric == 'auc028'  ~ 'AUC[0-28]*plain(", predicted")'
  ),
         metric = factor(metric,
                         levels = c('C[day-7]*plain(", empirical")',
                                    'C[day-7]*plain(", predicted")',
                                    'C[max]*plain(", predicted")',
                                    'AUC[0-28]*plain(", predicted")')) |> 
           as.character(),
         metric = parse_factor(metric))

surr_long |> 
  filter(metric != 'AUC[0-28]*plain(", predicted")') |> 
  ggplot() +
  geom_boxplot(aes(pqday_cat, # or specify as: reorder(pqday_cat, methb7)
                   value,
                   colour = factor(outcome7to120_pq)),
               varwidth = TRUE, linewidth = 0.5, width = 1.25, 
               position = position_dodge2(reverse = F),
               outlier.size = 0.2, outlier.alpha = 0.3) +
  scale_y_continuous(trans = scales::log_trans(),
                     breaks = c(0.025, 0.25, 0.5, 1, 2, 4, 8, 16, 32) / 100,
                     labels = function(x) {
                       pct <- scales::label_number(accuracy = 0.1, scale = 100)(x)
                       sub("\\.0$", "", pct)
                     },
                     limits = c(0.002, 0.4)
                     ) +
  facet_wrap(~metric, scales = 'free_y', nrow = 1, labeller = label_parsed) +
  theme(legend.position = 'top',
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(family = "Foundry Sterling", size = 10),
        aspect.ratio = 1,
        # panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent panel
        # plot.background = element_rect(fill = "transparent", colour = NA),   # Transparent plot area
        # legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend
        # legend.key = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_colour_manual(values = c('#918580', '#C82F46')) +
  labs(x = '\nDaily primaquine dose (mg/kg)',
       y = 'Methaemoglobin (%)\n',
       colour = '')

surr_long |> 
  filter(metric == 'AUC[0-28]*plain(", predicted")') |> 
  ggplot() +
  geom_boxplot(aes(pqday_cat, # or specify as: reorder(pqday_cat, methb7)
                   value * 100,
                   colour = factor(outcome7to120_pq)),
               varwidth = TRUE, linewidth = 0.5, width = 1.25, 
               position = position_dodge2(reverse = F),
               outlier.size = 0.2, outlier.alpha = 0.3) +
  scale_y_continuous(trans = scales::log_trans(),
                     breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
                     labels = function(x) {
                       pct <- scales::label_number(accuracy = 1)(x)
                       sub("\\.0$", "", pct)
                     },
                     limits = c(8, 512)
  ) +
  facet_wrap(~metric, scales = 'free_y', nrow = 1, labeller = label_parsed) +
  theme(legend.position = 'top',
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(family = "Foundry Sterling", size = 10),
        aspect.ratio = 1,
        # panel.background = element_rect(fill = "transparent", colour = NA),  # Transparent panel
        # plot.background = element_rect(fill = "transparent", colour = NA),   # Transparent plot area
        # legend.background = element_rect(fill = "transparent", colour = NA), # Transparent legend
        # legend.key = element_rect(fill = "transparent", colour = NA)
  ) +
  scale_colour_manual(values = c('#918580', '#C82F46')) +
  labs(x = '\nDaily primaquine dose (mg/kg)',
       y = 'Methaemoglobin (%•days)\n',
       colour = '')

# surrogacy_na_id <- surrogacy_data |> 
#   filter(is.na(sid_author)) |> 
#   select(pid) |> 
#   as.vector()
# write_rds(surrogacy_na_id, file = here('data', 'na_id.rds'))


# Logistic regression -----------------------------------------------------
mod_base_ <- glmer(outcome7to120_pq ~ log2(methb7) + 
                    age + sex + log_pvdens + pqmgkgday * pqdur_exp +
                    (1 | studysite),
                  data = surrogacy_4m,
                  family = binomial)
mod_base <- glmer(outcome7to120_pq ~ log2(methb7) + 
                  age + sex + log_pvdens + pqmgkgday + 
                  (1 | studysite),
                data = surrogacy_4m,
                family = binomial)
anova(mod_base_, mod_base)
# mod_base_x <- glmer(outcome7to120_pq ~ log2(methb7) + 
#                     age + sex + log_pvdens + pqmgkgday +
#                     (1 + log2(methb7) | studysite),
#                     data = surrogacy_4m,
#                     family = binomial)
broom.mixed::tidy(mod_base, exponentiate = T, conf.int = T) |> 
  select(term, estimate, conf.low, conf.high, p.value) |> 
  filter(term == 'log2(methb7)')
# anova(mod_base_x, mod_base, test = "Chisq")

mod_c7_ <- glmer(outcome7to120_pq ~ log2(c7) + 
                  age + sex + log_pvdens + pqmgkgday * pqdur_exp +
                  (1 | studysite),
                data = surrogacy_4m,
                family = binomial)
mod_c7 <- glmer(outcome7to120_pq ~ log2(c7) + 
                age + sex + log_pvdens + pqmgkgday +
                (1 | studysite),
                data = surrogacy_4m,
                family = binomial)
mod_c7x <- glmer(outcome7to120_pq ~ log2(c7) + 
                       age + sex + log_pvdens + #pqmgkgday +
                       (1 | studysite),
                     data = surrogacy_4m,
                     family = binomial)
anova(mod_c7x, mod_c7, test = "Chisq")
anova(mod_c7_, mod_c7, test = "Chisq")
broom.mixed::tidy(mod_c7x, exponentiate = T, conf.int = T) |> 
  select(term, estimate, conf.low, conf.high, p.value) |> 
  filter(term == 'log2(c7)')

mod_cmax_ <- glmer(outcome7to120_pq ~ log2(cmax) + 
                    age + sex + log_pvdens + pqmgkgday * pqdur_exp +
                    (1 | studysite),
                  data = surrogacy_4m,
                  family = binomial)
mod_cmax <- glmer(outcome7to120_pq ~ log2(cmax) + 
                  age + sex + log_pvdens + pqmgkgday +
                  (1 | studysite),
                data = surrogacy_4m,
                family = binomial)
broom.mixed::tidy(mod_cmax, exponentiate = T, conf.int = T) |> 
  select(term, estimate, conf.low, conf.high, p.value) |> 
  filter(term == 'log2(cmax)')

mod_auc014_ <- glmer(outcome7to120_pq ~ log2(auc014) + 
                      age + sex + log_pvdens + pqmgkgday * pqdur_exp +
                      (1 | studysite),
                    data = surrogacy_4m,
                    family = binomial)
mod_auc014 <- glmer(outcome7to120_pq ~ log2(auc014) + 
                    age + sex + log_pvdens + pqmgkgday +
                    (1 | studysite),
                  data = surrogacy_4m,
                  family = binomial)
broom.mixed::tidy(mod_auc014, exponentiate = T, conf.int = T) |> 
  select(term, estimate, conf.low, conf.high, p.value) |> 
  filter(term == 'log2(auc014)')

mod_auc028_ <- glmer(outcome7to120_pq ~ log2(auc028) + 
                     age + sex + log_pvdens + pqmgkgday * pqdur_exp +
                     (1 | studysite),
                     data = surrogacy_4m,
                     family = binomial)
mod_auc028 <- glmer(outcome7to120_pq ~ log2(auc028) + 
                      age + sex + log_pvdens + pqmgkgday +
                      (1 | studysite),
                    data = surrogacy_4m,
                    family = binomial)
mod_auc028x <- glmer(outcome7to120_pq ~ log2(auc028) + 
                      age + sex + log_pvdens + #pqmgkgday +
                      (1 | studysite),
                    data = surrogacy_4m,
                    family = binomial)
mod_auc028xy <- glmer(outcome7to120_pq ~ rcs(log2(auc028), 3) + 
                       age + sex + log_pvdens + #pqmgkgday +
                       (1 | studysite),
                     data = surrogacy_4m,
                     family = binomial)
anova(mod_auc028xy, mod_auc028x, test = "Chisq")
broom.mixed::tidy(mod_auc028xy, exponentiate = T, conf.int = T) |> 
  select(term, estimate, conf.low, conf.high, p.value) |> 
  filter(term == 'log2(auc028)')

AIC(mod_base, mod_c7, mod_cmax, mod_auc014, mod_auc028)
compare_performance(mod_base_, mod_c7_, mod_cmax_, mod_auc014_, mod_auc028_)
compare_performance(mod_base, mod_c7, mod_cmax, mod_auc014, mod_auc028)


# Cox regression ----------------------------------------------------------

coxme_base <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                    log2(methb7) + 
                    age + sex + log_pvdens + pqmgkgday * pqdur_exp + strata(trt1) +
                    (1 | studysite),
                    data = surrogacy_data,
                    ties = 'breslow')
(ph_coxme_base <- cox.zph(coxme_base))
summary(coxme_base)
confint(coxme_base) |> exp()

cox_base_brms <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(methb7) +
      age + sex + log_pvdens + pqmgkgday * pqdur_exp +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_base_brms,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2methb7')

cox_base_brms_7d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(methb7) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_7d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_base_brms_7d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2methb7')

cox_base_brms_14d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(methb7) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_14d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_base_brms_14d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2methb7')

coxme_c7 <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                    log2(c7) + 
                    age + sex + log_pvdens + pqmgkgday * pqdur_exp + strata(trt1) +
                    (1 | studysite),
                  data = surrogacy_data,
                  ties = 'breslow')
# anova(coxme_c7, coxme_c7_)
(ph_coxme_c7 <- cox.zph(coxme_c7))
summary(coxme_c7)
confint(coxme_c7) |> exp()

cox_c7_brms <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(c7) +
      age + sex + log_pvdens + pqmgkgday * pqdur_exp +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_c7_brms,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2c7')

cox_c7_brms_7d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(c7) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_7d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_c7_brms_7d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2c7')

cox_c7_brms_14d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(c7) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_14d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_c7_brms_14d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2c7')

coxme_cmax <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                      log2(cmax) + 
                      age + sex + log_pvdens + pqmgkgday * pqdur_exp + strata(trt1) +
                      (1 | studysite),
                    data = surrogacy_data,
                    ties = 'breslow')
(ph_coxme_cmax <- cox.zph(coxme_cmax))
summary(coxme_cmax)
confint(coxme_cmax) |> exp()

cox_cmax_brms <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(cmax) +
      age + sex + log_pvdens + pqmgkgday * pqdur_exp +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_cmax_brms,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2cmax')

cox_cmax_brms_7d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(cmax) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_7d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_cmax_brms_7d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2cmax')

cox_cmax_brms_14d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(cmax) +
      age + sex + log_pvdens + pqmgkgday +
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_14d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_cmax_brms_14d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2cmax')

coxme_auc014 <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                      log2(auc014) + 
                      age + sex + log_pvdens + pqmgkgday + strata(trt1) +
                      (1 | studysite),
                    data = surrogacy_data,
                    ties = 'breslow')
(ph_coxme_auc014 <- cox.zph(coxme_auc014))
summary(coxme_auc014)
confint(coxme_auc014) |> exp()

coxme_auc028 <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                        log2(methb7) + 
                        age + sex + log_pvdens + pqmgkgday * pqdur_exp + strata(trt1) +
                        (1 | studysite),
                      data = surrogacy_data,
                      ties = 'breslow')
coxme_auc028_ <- coxme(Surv(dlast120_pq, outcome7to120_pq) ~
                        log2(auc028) + 
                        age + sex + log_pvdens + pqmgkgday + 
                        strata(trt1) +
                        (1 | studysite),
                      data = surrogacy_data,
                      ties = 'breslow')
# anova(coxme_auc028, coxme_auc028_)
(ph_coxme_auc028 <- cox.zph(coxme_auc028))
summary(coxme_auc028)
confint(coxme_auc028) |> exp()

anova(coxme_base, coxme_c7, coxme_cmax, coxme_auc028, test = 'LRT')

compare_performance(coxme_base, coxme_c3, coxme_c7, coxme_cmax, coxme_auc014, coxme_auc028)

cox_auc028_brms <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(auc028) +
      age + sex + log_pvdens + pqmgkgday * pqdur_exp + 
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_auc028_brms,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2auc028')

cox_auc028_brms_7d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(auc028) +
      age + sex + log_pvdens + pqmgkgday + 
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_7d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_auc028_brms_7d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2auc028')

cox_auc028_brms_14d <- brm(
  formula = bf(
    dlast120_pq | cens(1 - outcome7to120_pq) ~
      log2(auc028) +
      age + sex + log_pvdens + pqmgkgday + 
      strata(trt1) +
      (1 | studysite)
  ),
  data = surrogacy_data_14d,
  family = brmsfamily("cox"),
  chains = 4, cores = 4, iter = 2000, seed = 13
)
tidy(cox_auc028_brms_14d,
     effects = "fixed", conf.int = TRUE, exponentiate = T) |> 
  select(term, estimate, conf.low, conf.high) |> 
  filter(term == 'log2auc028')

loo_base <- loo(cox_base_brms)
loo_auc028 <- loo(cox_auc028_brms)
loo_c7 <- loo(cox_c7_brms)
loo_compare(loo_base, loo_auc028, loo_c7)

kfold_base <- kfold(cox_base_brms, K = 10)
kfold_auc028 <- kfold(cox_auc028_brms, K = 10)
kfold_c7 <- kfold(cox_c7_brms, K = 10)
loo_compare(kfold_base, kfold_auc028, kfold_c7)

# Try next with pqdur exp multiplied by daily dose

# Population-level surrogacy ----------------------------------------------

nest_surd <- surrogacy_data |> 
  group_by(sid_author) |> 
  nest() |> 
  mutate(n = map_int(data, ~n_distinct(.x$pid))) |> 
  unnest(data)

n_sid <- surrogacy_data |> 
  group_by(sid_author) |> 
  nest() |> 
  mutate(n = map_int(data, ~n_distinct(.x$pid))) |> 
  select(sid_author, n)

# Pq on true endpoint
pop_cox <- nest_surd |> 
  # filter(sid_author != 'Pasaribu 2013') |> 
  group_by(sid_author) |> 
  nest() |> 
  mutate(
    model = map(data, ~ coxph(Surv(dlast120_pq, outcome7to120_pq) ~ pqmgkgtot, data = .x)),
    tidy_model = map(model, ~ tidy(.x, conf.int = TRUE))
  ) |> 
  select(sid_author, tidy_model) |> 
  unnest(tidy_model) |> 
  filter(term == 'pqmgkgtot') |> 
  select(sid_author, estimate, std.error) |> 
  rename(loghr = estimate,
         se_loghr = std.error)

pop_cox_brms <- nest_surd |>
  # filter(sid_author != 'Pasaribu 2013') |> 
  group_by(sid_author) |> 
  nest() |> 
  mutate(
    model = map(data, ~ brm(
      formula = bf(dlast120_pq | cens(1 - outcome7to120_pq) ~ pqmgkgtot + age + log_pvdens),
      data = .x,
      family = brmsfamily("cox"),
      chains = 2, iter = 2000, warmup = 1000, refresh = 0,
      silent = TRUE, seed = 123
    )),
    tidy_model = map(model, ~ tidy(.x, conf.int = TRUE, effects = "fixed"))
  ) |> 
  select(sid_author, tidy_model) |> 
  unnest(tidy_model) |> 
  filter(term == "pqmgkgtot") |> 
  select(sid_author, estimate, std.error) |> 
  rename(loghr = estimate,
         se_loghr = std.error)

# Empirical day-7 metHb
pop_lm <- nest_surd |>  
  # filter(sid_author != 'Pasaribu 2013') |> 
  group_by(sid_author) |> 
  nest() |> 
  mutate(
    model = map(data, ~ lm(log2(auc028) ~ pqmgkgtot + age + log_pvdens, data = .x)),
    tidy_model = map(model, ~ tidy(.x, conf.int = TRUE))
  ) |> 
  select(sid_author, tidy_model) |> 
  unnest(tidy_model) |> 
  filter(term == 'pqmgkgtot') |> 
  select(sid_author, estimate, std.error) |> 
  rename(md = estimate,
         se_md = std.error)

pop_lm_auc028 <- nest_surd |>  
  group_by(sid_author) |> 
  nest() |> 
  mutate(
    model = map(data, ~ brm(
      formula = log2(auc028) ~ pqmgkgtot + age + log_pvdens,
      data = .x,
      family = gaussian(),
      chains = 2, iter = 2000, warmup = 1000,
      refresh = 0, silent = TRUE, seed = 123
    )),
    tidy_model = map(model, ~ broom.mixed::tidy(.x, conf.int = TRUE))
  ) |> 
  select(sid_author, tidy_model) |> 
  unnest(tidy_model) |> 
  filter(term == "pqmgkgtot") |> 
  select(sid_author, estimate, std.error) |> 
  rename(md = estimate,
         se_md = std.error)

pop_auc028 <- pop_cox_brms |> left_join(pop_lm_auc028, by = "sid_author") |>
  ungroup() |> 
  left_join(n_sid, by = "sid_author")

pop_auc028$weight <- 1 / (pop_auc028$se_loghr^2 + pop_auc028$se_md^2)
res_auc028 <- TrialLevelMA(
  Alpha.Vector = pop_auc028$md,
  Beta.Vector  = pop_auc028$loghr,
  N.Vector     = pop_auc028$weight,
  Weighted     = TRUE  # More accurate if you have differing trial sizes
)

summary(res_auc028)
par(pty = "s"); plot(res_auc028)

# Fit a weighted linear model using trial sample sizes as weights
lm_model <- lm(loghr ~ md,
               data = pop_auc028, weights = weight)
summary(lm_model)

newdata <- data.frame(md = seq(min(pop_auc028$md - 0.04),
                               max(pop_auc028$md + 0.025), length.out = 100))

conf_df <- cbind(newdata,
                 as.data.frame(predict(lm_model,
                                       newdata = newdata, interval = "confidence")))

pred_df <- cbind(newdata,
                 as.data.frame(predict(lm_model,
                                       newdata = newdata, interval = "prediction")))


ggplot(data = pop_auc028, aes(x = md, y = loghr)) +
  
  # # Prediction interval (wider)
  # geom_ribbon(data = pred_df,
  #             aes(x = md, ymin = lwr, ymax = upr),
  #             fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  
  # Confidence interval (narrower)
  geom_point(aes(size = weight), alpha = 1, shape = 1) +
  scale_size_continuous(range = c(1, 20)) +
  # geom_ribbon(data = conf_df,
  #             aes(x = md, ymin = lwr, ymax = upr),
  #             fill = "grey50", alpha = 0.2, inherit.aes = FALSE) +
  
  # Fitted regression line
  geom_line(data = conf_df,
            aes(x = md, y = fit), linetype = 'dashed',
            colour = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  ggrepel::geom_text_repel(
    aes(label = sid_author),
    size = 1.85,
    direction = "y",
    family = "Foundry Sterling",
    nudge_y = -0.06,
    nudge_x = -0.035,
    vjust = 3,
    force = 0.5,
    max.overlaps = Inf,
    segment.colour = NA  # ← no line
  ) +
  scale_x_continuous(breaks = seq(-0.2, 0.5, by = 0.1),
                     # limits = c(-0.05, 0.45),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(-3, 1.5, by = 1),
                     limits = c(-2, 1),
                     # expand = expansion(mult = c(0, 0))
                     ) +
  
  labs(
    x = "\nMethaemoglobin mean difference",
    y = "Recurrence relative risk reduction\n"
  ) +
  theme(legend.position = 'none')
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



