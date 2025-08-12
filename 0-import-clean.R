
# Setup -------------------------------------------------------------------

# Packages
library(tidyverse)
library(haven)
library(here)
library(rms)
library(ggeffects)
library(lme4)

## Further refinement
theme_set(theme_bw())
theme_update(
  text = element_text(size = 10, family = "Foundry Sterling"), # Font
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
  aspect.ratio = 1,                          # Square
  legend.key = element_rect(fill = "transparent"),
  axis.title.y = element_text(hjust = 0.5)   # Move title for y-axis
)

# Data --------------------------------------------------------------------

methb_temp <- read_rds(file = here('data', 'methb_modelling.rds')) |> 
  filter(dayofobs_pq >= 0) |> 
  mutate(methb = if_else(methb == 0, 1, methb))  # Replace values of 0% with 1%, likely misreadings

methb_all.days <- methb_temp |>  
  filter(dayofobs_pq >= 0) |> 
  mutate(methb = if_else(methb == 0, 1, methb)) |>  # Replace values of 0% with 1%, likely misreadings
  # select(pid, dayofobs_pq) |> 
  group_by(pid) |> 
  arrange(dayofobs_pq) |> 
  complete(dayofobs_pq = full_seq(0:28, 1)) |> 
  ungroup() |> 
  filter(dayofobs_pq <= 28)

# Fill in the missing data of constant value variables
# Overall up to censoring (or conditional on some column)
methb <- methb_all.days |> 
  group_by(pid) |> 
  mutate( # Fill in constant values
    studysite = ifelse(is.na(studysite), unique(na.omit(studysite)), studysite),
    age = ifelse(is.na(age), unique(na.omit(age)), age),
    sex = ifelse(is.na(sex), unique(na.omit(sex)), sex),
    pqdur_exp = ifelse(is.na(pqdur_exp), unique(na.omit(pqdur_exp)), pqdur_exp),
    methb7 = ifelse(is.na(methb7), unique(na.omit(methb7)), methb7)
  ) |>
  mutate( # Fill in constant values up to some point
    pqmgkgday = case_when(
      pqdur_exp == "14 days" & dayofobs_pq <= 13 ~ ifelse(is.na(pqmgkgday), unique(na.omit(pqmgkgday[dayofobs_pq <= 13])), pqmgkgday),
      pqdur_exp == "14 days" & dayofobs_pq > 13 ~ 0,
      pqdur_exp == "7 days" & dayofobs_pq <= 6 ~ ifelse(is.na(pqmgkgday), unique(na.omit(pqmgkgday[dayofobs_pq <= 6])), pqmgkgday),
      pqdur_exp == "7 days" & dayofobs_pq > 6 ~ 0,
      TRUE ~ pqmgkgday
  )) |> 
  ungroup() |> 
  mutate(hb = 100 - methb,
         age_cat = case_when(age < 15 ~ 1,
                             age < 30 ~ 2,
                             age < 40 ~ 3,
                             age < 50 ~ 4,
                             age < 60 ~ 5,
                             age < 100 ~ 6,
                             is.na(age) ~ NA_integer_,
                             TRUE ~ 9999999),
         pid_ori = pid,
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
  mutate(pid = as.integer(as.factor(pid))) |> 
  group_by(pid_ori, dayofobs_pq) |>
  mutate(
    methb = if (all(is.na(methb))) NA_real_ else mean(methb, na.rm = TRUE)
  ) |>
  slice_head(n = 1) |>
  ungroup()
  
# Check data (random sample of size n)
set.seed(13)
id_selected <- filter(methb, dayofobs_pq == 0 & pqmgkgday > 0) |> pull(pid)
df_all <- filter(methb, pid %in% id_selected)
n_patients <- df_all$pid |> unique() |> length()
random_id20 <- sample(1:n_patients,
                      round(n_patients * 1, digits = 0)) # 1 means all patients in data 
n_patients_selected <- (df_all$pid |> unique())[random_id20]
random_patient20 <- filter(df_all, pid %in% n_patients_selected)

dose_group <- random_patient20 |> 
  select(pid_ori, dayofobs_pq, pqmgkgday) |> 
  filter(dayofobs_pq == 0) |> 
  distinct(pid_ori, .keep_all = TRUE) |> 
  select(pid_ori, pqmgkgday) |> 
  mutate(dose_group = case_when(
    pqmgkgday < 0.375 ~ 1,
    pqmgkgday < 0.75 ~ 2,
    TRUE ~ 3
    ) |> as.factor()) |> 
  select(pid_ori, dose_group)

patient_data <- random_patient20 |> left_join(dose_group, by = 'pid_ori')

# Processed data ----------------------------------------------------------

write.csv(patient_data,
          here("data", "patient_data.csv"),
          row.names = F)
























