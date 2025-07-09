
# load libraries and necessary data ---------------------------------------

library(tidyverse)
library(tidybayes)
library(ubms)
library(patchwork)
library(CAbioacoustics)
library(CAbioacousticsextras)
library(sf)
library(ggh4x)
library(tidyterra)
library(geodata)
library(cowplot)
library(ggspatial)
library(rnaturalearth)
library(ggsflabel)
library(scales)
library(ggplot.multistats)

# load model
fit_re <-
  readr::read_rds(here::here('outputs/fit_re.rds'))

# survey years
start_year <- 2021
end_year <- 2024

# study area polygon
study_area_sf <-
  cb_get_spatial('sierra_study_area')

# USFS forest boundaries
forests_sf <-
  cb_get_spatial('usfs_boundaries')

# save figures here
dir_create(here::here('figures'))


# posterior distribution plot (all submodels combined) --------------------

# see model parameters for extracting properly below
head(get_variables(fit_re@stanfit))

# detection submodel
det_parameters_df <-
  fit_re@stanfit |>
  tidy_draws() |>
  select(1:3, matches('beta_det')) |>
  pivot_longer(
    !c('.chain', '.iteration', '.draw'),
    names_to = 'parameter',
    values_to = 'value'
  ) |>
  # don't need to show intercept
  filter(!str_detect(parameter, 'Intercept')) |>
  # rename parameters
  mutate(
    parameter = case_when(
      parameter == 'beta_det[log_survey_hours_scaled]' ~ 'ARU hours',
      parameter == 'beta_det[date_scaled]' ~ 'Date',
      parameter == 'beta_det[date_scaled2]' ~ 'Date2',
      parameter == 'beta_det[year2022]' ~ 'Year (2022)',
      parameter == 'beta_det[year2023]' ~ 'Year (2023)',
      parameter == 'beta_det[year2024]' ~ 'Year (2024)'
    ),
    submodel = 'detection'
  )

# state submodel
state_parameters_df <-
  fit_re@stanfit |>
  tidy_draws() |>
  select(1:3, matches('beta_state')) |>
  pivot_longer(
    !c('.chain', '.iteration', '.draw'),
    names_to = 'parameter',
    values_to = 'value'
  ) |>
  filter(!str_detect(parameter, 'Intercept')) |>
  mutate(
    parameter = case_when(
      parameter == 'beta_state[p_high_severity_scaled]' ~ 'Prop. high severity burn',
      parameter == 'beta_state[y_coord_scaled2]' ~ 'Latitude2',
      parameter == 'beta_state[y_coord_scaled]' ~ 'Latitude',
      parameter == 'beta_state[elev_lat_resid_scaled]' ~ 'Elevation',
      parameter == 'beta_state[elev_lat_resid_scaled2]' ~ 'Elevation2',
      parameter == 'beta_state[cfo_ch_scaled]' ~ 'Canopy height'
    ),
    submodel = 'state'
  )

# extinction submodel
ext_parameters_df <-
  fit_re@stanfit |>
  tidy_draws() |>
  select(1:3, matches('beta_ext')) |>
  pivot_longer(
    !c('.chain', '.iteration', '.draw'),
    names_to = 'parameter',
    values_to = 'value'
  ) |>
  filter(!str_detect(parameter, 'Intercept')) |>
  mutate(
    parameter = case_when(
      parameter == 'beta_ext[scale(p_high_severity)]' ~ 'Prop. high severity burn',
      parameter == 'beta_ext[year2022]' ~ 'Year (2022)',
      parameter == 'beta_ext[year2023]' ~ 'Year (2023)'
    ),
    submodel = 'extinction'
  )

# colonization submodel
col_parameters_df <-
  fit_re@stanfit |>
  tidy_draws() |>
  select(1:3, matches('beta_col')) |>
  pivot_longer(
    !c('.chain', '.iteration', '.draw'),
    names_to = 'parameter',
    values_to = 'value'
  ) |>
  filter(!str_detect(parameter, 'Intercept')) |>
  mutate(
    parameter = case_when(
      parameter == 'beta_col[scale(p_high_severity)]' ~ 'Prop. high severity burn',
      parameter == 'beta_col[year2022]' ~ 'Year (2022)',
      parameter == 'beta_col[year2023]' ~ 'Year (2023)'
    ),
    submodel = 'colonization'
  )

# combine submodels
posteriors_df <-
  det_parameters_df |>
  bind_rows(state_parameters_df) |>
  bind_rows(ext_parameters_df) |>
  bind_rows(col_parameters_df)

# plot all posteriors together
posteriors_df |>
  mutate(submodel = str_to_title(submodel)) |>
  mutate(submodel = fct_relevel(submodel, 'Extinction', 'Colonization', 'State', 'Detection')) |>
  mutate(parameter = str_c(parameter, submodel)) |>
  group_by(submodel) |>
  mutate(parameter = fct_reorder(parameter, value)) |>
  ungroup() |>
  ggplot(aes(y = parameter, x = value, fill = submodel)) +
  stat_halfeye(point_interval = 'mean_qi', .width = 0.95, linewidth = 2, size = 2) +
  labs(
    x = 'Effect size',
    y = 'Parameter',
    title = 'Submodel'
  ) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = 'top',
    plot.title = element_text(vjust = -1)
  ) +
  ggokabeito::scale_fill_okabe_ito(
    order = c(1, 2, 3, 7),
    name = NULL,
    breaks = c('Detection', 'Colonization', 'State', 'Extinction'),
    labels = c(
      expression("Detection ("* italic(p) *")"),
      expression("Colonization"~(gamma)),
      expression("Occupancy"~(psi[2021])),
      expression("Extinction"~(epsilon)))
  ) +
  # more renaming to clean up
  scale_y_discrete(
    labels = c(
      "Prop. high severity burnState" = "High severity burn",
      "Prop. high severity burnColonization" = "High severity burn",
      "Prop. high severity burnExtinction" = "High severity burn",
      "Year (2022)Colonization" = "Year (2022)",
      "Year (2023)Colonization" = "Year (2023)",
      "Year (2022)Extinction" = "Year (2022)",
      "Year (2023)Extinction" = "Year (2023)",
      "Year (2022)Detection" = "Year (2022)",
      "Year (2023)Detection" = "Year (2023)",
      "Year (2024)Detection" = "Year (2024)",
      "Date2Detection" = expression("Date"^2),
      "ARU hoursDetection" = "ARU hours",
      "DateDetection" = "Date",
      "Canopy heightState" = "Canopy height",
      "Latitude2State" =  expression("Latitude"^2),
      "Elevation2State" = expression("Elevation"^2),
      "LatitudeState" = "Latitude",
      "ElevationState" = "Elevation"
    )
  )
ggsave(here::here('figures/effect_sizes.png'), width = 5, height = 8, units = 'in', dpi = 600)


# fire effect on extinction -----------------------------------------------

# Create a prediction data frame
fire_seq <-
  seq(
    min(umf@siteCovs$p_high_severity),
    max(umf@siteCovs$p_high_severity),
    length.out = 100
  )

# Scale it exactly as before
fire_scaled <- (fire_seq - mean(umf@siteCovs$p_high_severity)) / sd(umf@siteCovs$p_high_severity)

newdata <-
  tibble(
    year = as.factor(2022),
    p_high_severity_scaled = fire_scaled,
    p_high_severity = fire_seq
  )

p_ext <-
  predict(fit_re, submodel = "ext", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(p_high_severity, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(p_high_severity, Predicted), linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Proportion burned at high severity',
    y = expression(epsilon),
    title = 'c)'
  )


# fire effect on colonization ---------------------------------------------

p_col <-
  predict(fit_re, submodel = "col", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(p_high_severity, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(p_high_severity, Predicted), linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Proportion burned at high severity',
    y = expression(gamma),
    title = 'b)'
  )


# fire effect on initial occupancy ----------------------------------------

# Create a prediction data frame
fire_seq <-
  seq(
    min(umf@siteCovs$p_high_severity),
    max(umf@siteCovs$p_high_severity),
    length.out = 100
  )

# Scale it exactly as before
fire_scaled <- (fire_seq - mean(umf@siteCovs$p_high_severity)) / sd(umf@siteCovs$p_high_severity)

newdata <-
  tibble(
    elev_lat_resid_scaled = mean(umf@siteCovs$elev_lat_resid_scaled),
    elev_lat_resid_scaled2 = mean(umf@siteCovs$elev_lat_resid_scaled2),
    y_coord_scaled = mean(umf@siteCovs$y_coord_scaled),
    y_coord_scaled2 = mean(umf@siteCovs$y_coord_scaled2),
    cfo_ch_scaled = mean(umf@siteCovs$cfo_ch_scaled),
    p_high_severity_scaled = fire_scaled,
    p_high_severity = fire_seq
  )

# Predict initial occupancy (psi)
pred <-
  predict(fit_re, submodel = "state", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  as_tibble()

p_state <-
  pred |>
  ggplot() +
  geom_ribbon(aes(p_high_severity, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(p_high_severity, Predicted), linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Proportion burned at high severity',
    y = expression(psi[2021]),
    title = 'a)'
  )


# combined fire effects plot ----------------------------------------------

p_state + p_col + p_ext + plot_layout(axis_titles = "collect")
ggsave(here::here('figures/fire_effects.png'), width = 8, height = 3, units = 'in', dpi = 600)


# state effects in addition to fire ---------------------------------------

# canopy height
cfo_seq <-
  seq(
    min(umf@siteCovs$cfo_ch),
    max(umf@siteCovs$cfo_ch),
    length.out = 100
  )

# Scale it exactly as before
cfo_scaled <- (cfo_seq - mean(umf@siteCovs$cfo_ch)) / sd(umf@siteCovs$cfo_ch)

newdata <-
  tibble(
    elev_lat_resid_scaled = mean(umf@siteCovs$elev_lat_resid_scaled),
    elev_lat_resid_scaled2 = mean(umf@siteCovs$elev_lat_resid_scaled2),
    y_coord_scaled = mean(umf@siteCovs$y_coord_scaled),
    y_coord_scaled2 = mean(umf@siteCovs$y_coord_scaled2),
    cfo_ch_scaled = cfo_scaled,
    cfo_ch = cfo_seq,
    p_high_severity_scaled = mean(umf@siteCovs$p_high_severity_scaled)
  )

p_state_cfo <-
  predict(fit_re, submodel = "state", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(cfo_ch, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(cfo_ch, Predicted), linewidth = 1.5) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Canopy height (m)',
    y = expression(psi[2021]),
    title = 'a)'
  )

# elevation
ele_seq <-
  seq(
    min(umf@siteCovs$elev_lat_resid),
    max(umf@siteCovs$elev_lat_resid),
    length.out = 100
  )

# Scale it exactly as before
ele_scaled <- (ele_seq - mean(umf@siteCovs$elev_lat_resid)) / sd(umf@siteCovs$elev_lat_resid)

newdata <-
  tibble(
    elev_lat_resid_scaled = ele_scaled,
    elev_lat_resid_scaled2 = ele_scaled^2,
    y_coord_scaled = mean(umf@siteCovs$y_coord_scaled),
    y_coord_scaled2 = mean(umf@siteCovs$y_coord_scaled2),
    cfo_ch_scaled = mean(umf@siteCovs$cfo_ch_scaled),
    ele_seq = ele_seq,
    p_high_severity_scaled = mean(umf@siteCovs$p_high_severity_scaled)
  )

p_state_elev <-
  predict(fit_re, submodel = "state", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(ele_seq, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(ele_seq, Predicted), linewidth = 1.5) +
  # scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Elevation ~ latitude residual',
    y = expression(psi),
    title = 'b)'
  )

# latitude
lat_seq <-
  seq(
    min(umf@siteCovs$y_coord),
    max(umf@siteCovs$y_coord),
    length.out = 100
  )

# Scale it exactly as before
lat_scaled <- (lat_seq - mean(umf@siteCovs$y_coord)) / sd(umf@siteCovs$y_coord)

newdata <-
  tibble(
    y_coord_scaled = lat_scaled,
    y_coord_scaled2 = lat_scaled^2,
    elev_lat_resid_scaled = mean(umf@siteCovs$elev_lat_resid_scaled),
    elev_lat_resid_scaled2 = mean(umf@siteCovs$elev_lat_resid_scaled2),
    cfo_ch_scaled = mean(umf@siteCovs$cfo_ch_scaled),
    p_high_severity_scaled = mean(umf@siteCovs$p_high_severity_scaled),
    lat_seq = lat_seq
  )

p_state_lat <-
  predict(fit_re, submodel = "state", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(lat_seq, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(lat_seq, Predicted), linewidth = 1.5) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Latitude (CA Albers northing)',
    y = expression(psi),
    title = 'c)'
  )

p_state_cfo + p_state_elev + p_state_lat
ggsave(here::here('figures/state_effects.png'), width = 8, height = 3, units = 'in', dpi = 600)


# detection probability plots ---------------------------------------------

# date
date_seq <-
  seq(
    min(umf@obsCovs$date, na.rm = TRUE),
    max(umf@obsCovs$date, na.rm = TRUE),
    length.out = 100
  )

# Scale it exactly as before
date_scaled <- (date_seq - mean(umf@obsCovs$date, na.rm = TRUE)) / sd(umf@obsCovs$date, na.rm = TRUE)

newdata <-
  tibble(
    date_scaled = date_scaled,
    date_scaled2 = date_scaled^2,
    log_survey_hours_scaled = mean(umf@obsCovs$log_survey_hours_scaled, na.rm = TRUE),
    year = '2022',
    date = date_seq
  )

p_det_date <-
  predict(fit_re, submodel = "det", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  mutate(date = as.Date(date, origin = as.Date('2021-04-07'))) |>
  ggplot() +
  geom_ribbon(aes(date, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(date, Predicted), linewidth = 1.5) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Date',
    y = 'Detection probability',
    title = 'b)'
  )

# survey_hours
effort_seq <-
  seq(
    min(umf@obsCovs$survey_hours, na.rm = TRUE),
    max(umf@obsCovs$survey_hours, na.rm = TRUE),
    length.out = 100
  )

survey_seq_log <- log(effort_seq)

# Scale it exactly as before
survey_scaled <- (survey_seq_log - mean(log(umf@obsCovs$survey_hours), na.rm = TRUE)) / sd(log(umf@obsCovs$survey_hours), na.rm = TRUE)

newdata <-
  tibble(
    date_scaled = mean(umf@obsCovs$date_scaled, na.rm = TRUE),
    date_scaled2 = mean(umf@obsCovs$date_scaled2, na.rm = TRUE),
    log_survey_hours_scaled = survey_scaled,
    survey_hours = effort_seq,
    year = '2022'
  )

p_det_survey_hours <-
  predict(fit_re, submodel = "det", newdata = newdata, re.form = NA, level = c(0.95)) |>
  bind_cols(newdata) |>
  ggplot() +
  geom_ribbon(aes(survey_hours, ymin = `2.5%`, ymax = `97.5%`), fill = 'grey70', alpha = 8/10) +
  geom_line(aes(survey_hours, Predicted), linewidth = 1.5) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'ARU hours per week',
    y = expression(paste("Detection probability (", italic(p), ")")),
    title = 'a)'
  )

# year
newdata <-
  tibble(
    year = as.factor(start_year:end_year),
    date_scaled = mean(umf@obsCovs$date_scaled, na.rm = TRUE),
    date_scaled2 = mean(umf@obsCovs$date_scaled2, na.rm = TRUE),
    log_survey_hours_scaled = mean(umf@obsCovs$log_survey_hours_scaled, na.rm = TRUE),
    survey_hours = mean(umf@obsCovs$survey_hours, na.rm = TRUE)
  )

samples <- 1:nsamples(fit_re)

p_det_year <-
  # predict(fit_re, submodel = "det", newdata = newdata, re.form = NA) |>
  # bind_cols(newdata) |>
  ubms:::sim_lp(
    fit_re,
    submodel = "det",
    newdata = newdata,
    re.form = NA,
    transform = TRUE,
    samples = samples
  ) |>
  as_tibble() |>
  pivot_longer(
    everything()
  ) |>
  mutate(
    year = as.numeric(str_extract(name, "(\\d)+")) + 2020
  ) |>
  ggplot() +
  stat_halfeye(aes(x = as.factor(year), y = value), alpha = 10/10, point_interval = 'mean_qi', .width = c(0.95), point_size = 1.25, size = 2.25) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  ylim(c(0, 1)) +
  labs(
    x = 'Year',
    y = NULL,
    title = 'c)'
  )

p_det_survey_hours + p_det_date + p_det_year
ggsave(here::here('figures/det_effects.png'), width = 8, height = 3, units = 'in', dpi = 600)

predict(fit_re, submodel = "det", newdata = newdata, re.form = NA) |>
  bind_cols(newdata) |>
  rename(lower = `2.5%`, upper = `97.5%`) |>
  write_rds(here::here('outputs/detection_probabilities_df.rds'))

# for year with lowest detection probability,
min_newdata <-
  predict(fit_re, submodel = "det", newdata = newdata, re.form = NA) |>
  bind_cols(newdata) |>
  filter(Predicted == min(Predicted))

samples <- 1:nsamples(fit_re)

lp <-
  ubms:::sim_lp(
    fit_re,
    submodel = "det",
    newdata = min_newdata,
    re.form = NA,
    transform = TRUE,
    samples = samples
  ) |>
  as_tibble() |>
  transmute(
    prob_deployment = 1 - (1 - V1)^5
  ) |>
  mean_qi(prob_deployment) |>
  bind_cols(min_newdata |> select(year))

lp |>
  write_rds(here::here('outputs/deployment_probability_detection.rds'))

ubms:::sim_lp(
  fit_re,
  submodel = "det",
  newdata = min_newdata,
  re.form = NA,
  transform = TRUE,
  samples = samples
) |>
  as_tibble() |>
  transmute(
    prob_deployment = 1 - (1 - V1)^5
  ) |>
  ggplot() +
  stat_halfeye(aes(x = prob_deployment), alpha = 10/10, point_interval = 'mean_qi', .width = c(0.95), fill = 'violet', point_size = 3) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(
    x = 'Detection probability over 5-week deployment',
    y = NULL
  )
ggsave(here::here('figures/detection_over_deployment.png'), width = 8, height = 3, units = 'in', dpi = 600)


# occupancy by forest -----------------------------------------------------

proj_occ <- projected(fit_re)

# Re-arrange into array samples x primary period x nsite
proj_occ <- array(proj_occ, c(nsamples(fit_re), umf@numPrimary, numSites(umf)))

sites <- 1:dim(proj_occ)[3]

occ_list <- list()

for (i in seq_along(sites)) {

  cell_df <-
    umf@siteCovs |>
    # fire_name provided in original dataset
    select(cell_id, forest_name, fire_name) |>
    slice(sites[i])

  occ_list[[i]] <-
    proj_occ[,,sites[i]] |>
    as_tibble() |>
    mutate(posterior_sample = row_number()) |>
    pivot_longer(
      !posterior_sample,
      names_to = 'year',
      values_to = 'psi'
    ) |>
    mutate(
      year = as.numeric(str_remove(year, 'V')) + 2020
    ) |>
    bind_cols(cell_df)

}

# combine all sites
occ_df <- bind_rows(occ_list)

# occupancy across the Sierra Nevada
range_wide_occ_df <-
  occ_df |>
  group_by(posterior_sample, year) |>
  summarise(psi = mean(psi)) |>
  group_by(year) |>
  summarise(
    occ = mean(psi),
    sd = sd(psi),
    lower = quantile(psi, 0.025),
    upper = quantile(psi, 0.975)
  ) |>
  mutate(forest_name = 'bold(`All forests`)')

# occupancy by forest
forest_occ_df <-
  occ_df |>
  group_by(posterior_sample, year, forest_name) |>
  summarise(psi = mean(psi)) |>
  group_by(year, forest_name) |>
  summarise(
    occ = mean(psi),
    sd = sd(psi),
    lower = quantile(psi, 0.025),
    upper = quantile(psi, 0.975)
  ) |>
  ungroup()

# combine these for plotting
forest_occ_summary_df <-
  forest_occ_df |>
  bind_rows(range_wide_occ_df) |>
  mutate(year = as.factor(year)) |>
  mutate(forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest')))

y_max <- round(max(forest_occ_summary_df$upper), 2) + 0.03

# plot
forest_occ_summary_df |>
  ggplot() +
  geom_pointrange(aes(x = year, y = occ, ymin = lower, ymax = upper), fatten = 6/10, alpha = 8/10) +
  geom_line(aes(year, occ, group = forest_name), linetype = 2, linewidth = 4/10) +
  facet_wrap(~forest_name, ncol = 4, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    base_size = 8
  ) +
  labs(
    x = 'Year',
    y = expression(psi)
  )
ggsave(here::here('figures/occupancy_trends_by_forest.png'), width = 5, height = 3, units = 'in', dpi = 600)

# save these for supplemental table
forest_occ_summary_df |>
  arrange(forest_name, year) |>
  mutate(
    forest_name = case_when(
      forest_name == 'bold(`All forests`)' ~ 'All forests',
      TRUE ~ forest_name
    )
  ) |>
  write_csv(here::here('outputs/annual_occupancy_estimates_forest.csv'))

# CV
mean(
  forest_occ_summary_df |>
    arrange(forest_name, year) |>
    mutate(
      forest_name = case_when(
        forest_name == 'bold(`All forests`)' ~ 'All forests',
        TRUE ~ forest_name
      )
    ) |>
    summarise(
      cv = sd / occ
    ) |>
    pull(cv)
)


# occupancy by forest (tidybayes versions) --------------------------------

# occupancy across the Sierra Nevada
range_wide_occ_df <-
  occ_df |>
  group_by(posterior_sample, year) |>
  summarise(psi = mean(psi)) |>
  # group_by(year) |>
  # summarise(
  #   occ = mean(psi),
  #   lower = quantile(psi, 0.025),
  #   upper = quantile(psi, 0.975)
  # ) |>
  mutate(forest_name = 'italic("All forests")')

# occupancy by forest
forest_occ_df <-
  occ_df |>
  group_by(posterior_sample, year, forest_name) |>
  summarise(psi = mean(psi))

# combine these for plotting
forest_occ_summary_df <-
  forest_occ_df |>
  bind_rows(range_wide_occ_df) |>
  mutate(year = as.factor(year)) |>
  mutate(forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest')))

forest_occ_summary_mean_df <-
  forest_occ_summary_df|>
  group_by(year, forest_name) |>
  summarise(
    occ = mean(psi),
    lower = quantile(psi, 0.025),
    upper = quantile(psi, 0.975)
  ) |>
  mutate(year = as.factor(year)) |>
  mutate(forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest')))

# plot
forest_occ_summary_df |>
  mutate(forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest'))) |>
  ggplot() +
  stat_halfeye(aes(x = year, y = psi), alpha = 8/10, point_interval = 'mean_qi', .width = c(0.95), point_size = 6/10, size = 3/4) +
  geom_line(data = forest_occ_summary_mean_df, aes(year, occ, group = forest_name), linetype = 2, linewidth = 4/10, alpha = 7/10) +
  facet_wrap(~forest_name, ncol = 4, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, y_max)) +
  theme_light(base_size = 8) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = 'Year',
    y = expression(psi)
  )
ggsave(here::here('figures/occupancy_trends_by_forest_tidybayes.png'), width = 5, height = 3, units = 'in', dpi = 600)


# trends in occupancy by forest (lambda avg) ------------------------------

label_perc_df <-
  forest_occ_summary_df |>
  mutate(
    forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest'))
  ) |>
  pivot_wider(
    names_from = year,
    values_from = psi
  ) |>
  ungroup() |>
  rowwise() |>
  mutate(
    lambda_1 = `2022` / `2021`,
    lambda_2 = `2023` / `2022`,
    lambda_3 = `2024` / `2023`,
    lambda_avg = mean(c(lambda_1, lambda_2, lambda_3))
  ) |>
  mutate(
    forest_name = case_when(
      str_detect(forest_name, 'forests') ~ '<i>All forests</i>',
      TRUE ~ forest_name
    )
  ) |>
  ungroup() |>
  mutate(
    direction = case_when(
      lambda_avg < 1 ~ 'decreasing',
      TRUE ~ 'stable/increasing'
    )
  ) |>
  group_by(forest_name, direction) |>
  tally() |>
  mutate(
    perc = round(((n / sum(n)) * 100), 0),
    prob = round((n / sum(n)), 2)
  ) |>
  filter(perc == max(perc)) |>
  mutate(perc = str_c(perc, '%')) |>
  ungroup() |>
  select(forest_name, perc, prob)

label_x_df <-
  forest_occ_summary_df |>
  mutate(
    forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest'))
  ) |>
  pivot_wider(
    names_from = year,
    values_from = psi
  ) |>
  ungroup() |>
  rowwise() |>
  mutate(
    lambda_1 = `2022` / `2021`,
    lambda_2 = `2023` / `2022`,
    lambda_3 = `2024` / `2023`,
    lambda_avg = mean(c(lambda_1, lambda_2, lambda_3))
  ) |>
  mutate(
    forest_name = case_when(
      str_detect(forest_name, 'forests') ~ '<i>All forests</i>',
      TRUE ~ forest_name
    )
  ) |>
  ungroup() |>
  group_by(forest_name) |>
  summarise(
    lambda_avg = mean(lambda_avg)
  )

label_df <-
  label_perc_df |>
  left_join(label_x_df)

forest_occ_summary_df |>
  mutate(
    forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest'))
  ) |>
  pivot_wider(
    names_from = year,
    values_from = psi
  ) |>
  ungroup() |>
  rowwise() |>
  mutate(
    lambda_1 = `2022` / `2021`,
    lambda_2 = `2023` / `2022`,
    lambda_3 = `2024` / `2023`,
    lambda_avg = mean(c(lambda_1, lambda_2, lambda_3))
  ) |>
  mutate(
    forest_name = case_when(
      str_detect(forest_name, 'forests') ~ '<i>All forests</i>',
      TRUE ~ forest_name
    )
  ) |>
  ungroup() |>
  group_by(forest_name) |>
  summarise(
    mean = mean(lambda_avg),
    sd = sd(lambda_avg),
    .lower = quantile(lambda_avg, 0.025),
    .upper = quantile(lambda_avg, 0.975)
  ) |>
  rename(lambda_avg = mean) |>
  select(forest_name:.upper) |>
  left_join(label_df) |>
  write_rds(here::here('outputs/lambda_avg_df.rds'))

forest_occ_summary_df |>
  mutate(
    forest_name = fct_relevel(forest_name, CAbioacoustics::forests_north_south |> str_remove(' National Forest'))
  ) |>
  pivot_wider(
    names_from = year,
    values_from = psi
  ) |>
  ungroup() |>
  rowwise() |>
  mutate(
    lambda_1 = `2022` / `2021`,
    lambda_2 = `2023` / `2022`,
    lambda_3 = `2024` / `2023`,
    lambda_avg = mean(c(lambda_1, lambda_2, lambda_3))
  ) |>
  mutate(
    forest_name = case_when(
      str_detect(forest_name, 'forests') ~ '<i>All forests</i>',
      TRUE ~ forest_name
    )
  ) |>
  ungroup() |>
  mutate(forest_name = fct_relevel(forest_name, 'Lassen', 'Plumas', 'Tahoe', 'Eldorado', 'Stanislaus', 'Sierra', 'Sequoia', '<i>All forests</i>')) |>
  mutate(forest_name = fct_relevel(forest_name, rev)) |>
  ggplot() +
  geom_vline(aes(xintercept = 1), linetype = 2, color = 'grey55') +
  stat_slab(aes(y = forest_name, x = lambda_avg, fill = after_stat(x < 1))) +
  stat_pointinterval(aes(y = forest_name, x = lambda_avg), point_interval = 'mean_qi', .width = c(0.95), point_size = 2, size = 2) +
  geom_text(data = label_df, aes(y = forest_name, x = lambda_avg, label = perc), color = 'black', size = 3, nudge_y = 2/10) +
  scale_fill_manual(
    values = c('TRUE' = 'tomato3', 'FALSE' = 'gray80'),  # Set the fill colors
    breaks = rev(c(TRUE, FALSE)),
    name = 'Trend',
    labels = c('Stable or increasing', 'Decreasing')
  ) +  # Reverse the legend order
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = ggtext::element_markdown(),
    legend.position = 'top',
    plot.title = element_text(vjust = -1)
  ) +
  labs(
    y = 'National Forest',
    x = expression(bar(lambda))
  ) +
  scale_x_continuous(limits = c(0.5, 1.5))
ggsave(here::here('figures/lambda_avg.png'), width = 4, height = 6, units = 'in', dpi = 600)


# occupancy by fire -------------------------------------------------------

# get hexes that were monitored
focal_hexes_sf <-
  cb_get_spatial('sierra_hexes') |>
  filter(cell_id %in% umf@siteCovs$cell_id)

fires_for_occupancy_sf <-
  arcpullr::get_spatial_layer(
    url = 'https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/ArcGIS/rest/services/California_Historic_Fire_Perimeters/FeatureServer/2',
    where = glue::glue("YEAR_ >= {start_year - 1} AND YEAR_ < {end_year}")
  )

fires_for_occupancy_sf <-
  fires_for_occupancy_sf |>
  janitor::clean_names() |>
  st_make_valid() |>
  st_intersection(
    cb_get_spatial('sierra_study_area')
  ) |>
  select(year, fire_name)

fires_to_keep <-
  focal_hexes_sf |>
  st_intersection(fires_for_occupancy_sf) |>
  st_drop_geometry() |>
  group_by(fire_name, year) |>
  tally() |>
  ungroup() |>
  arrange(desc(n)) |>
  # keep fires that have at least 5 cells
  filter(n >= 5) |>
  arrange(year, n)

# pull out recent fires
fire_perimeters_sf <-
  arcpullr::get_spatial_layer(
    url = 'https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/ArcGIS/rest/services/California_Historic_Fire_Perimeters/FeatureServer/2',
    where = glue::glue("YEAR_ >= {start_year - 1} AND YEAR_ < {end_year}")
  ) |>
  janitor::clean_names() |>
  st_make_valid() |>
  group_by(year, fire_name) |>
  summarise() |>
  ungroup() |>
  inner_join(fires_to_keep)

fire_size_df <-
  fire_perimeters_sf |>
  st_transform(3310) %>%
  mutate(size = as.numeric(st_area(.))) |>
  select(year, fire_name, size) |>
  st_drop_geometry() |>
  mutate(fire_name = str_to_title(fire_name)) |>
  mutate(fire_name = case_when(
    fire_name == 'Knp Complex' ~ 'KNP Complex',
    TRUE ~ fire_name
  ))

# figure out which monitored cells were in a fire footprint
burned_hexes <-
  focal_hexes_sf |>
  st_intersection(fire_perimeters_sf) |>
  group_by(cell_id) |>
  add_count() |>
  filter(nn == 1) |>
  select(cell_id, fire_name, year) |>
  st_drop_geometry()

# find these burned hexes again
focal_hexes_sf <-
  cb_get_spatial('sierra_hexes') |>
  inner_join(burned_hexes)

# now get df of % burned at high severity
x <- list()

hexes <- focal_hexes_sf$cell_id

# loop through each year a cell is burned
for (i in seq_along(hexes)) {

  year <- focal_hexes_sf |> filter(cell_id == hexes[i]) |> pull(year)
  fire <- focal_hexes_sf |> filter(cell_id == hexes[i]) |> pull(fire_name)

  # read in severity raster based on fire year; path will change depending on name given when mapping drive
  fire_rast <- rast(str_c(here::here("Y:/Data-GIS/Fires/cbi_mapping_sierras/cbi_sierra_cat_rasters"), '/', 'cbi_cat_', {year}, '.tif'))

  # crop and mask it to hex
  fire_rast2 <- terra::crop(fire_rast, focal_hexes_sf |> filter(cell_id == hexes[i]) |> st_transform(st_crs(fire_rast)), mask = TRUE)

  # plot to check it
  plot(fire_rast2, main = str_c(hexes[i]))
  plot(st_geometry(focal_hexes_sf |> filter(cell_id == hexes[i]) |> st_transform(st_crs(fire_rast2))), add = TRUE)

  # extract the severity info for all raster cells within the hex
  # get tibble of cell ID and fire severity
  e <-
    terra::extract(fire_rast2, focal_hexes_sf |> filter(cell_id == hexes[i]) |> st_transform(st_crs(fire_rast2)), cellnumbers = TRUE) |>
    as_tibble()

  # save for each fire year
  x[[i]] <-
    e |>
    select(-ID) |>
    mutate(cell_id = hexes[i])

}

# now combine all year(s) for that cell
x <-
  x |>
  # year(s) are in columns ordered by cell #
  bind_rows()

# create burned categories of >= 50% or <50%
burned_cells_df <-
  x |>
  group_by(cell_id, CBI_bc) |>
  tally() |>
  mutate(prop = n / sum(n)) |>
  filter(CBI_bc == 3) |>
  ungroup() |>
  mutate(
    category = case_when(
      prop >= 0.5 & CBI_bc == 3 ~ "\u22650.5",
      prop > 0 & prop < 0.5 & CBI_bc == 3 ~ '<0.5'
    )
  )

# get list of fires and associated NFs
fires_forests_df <-
  fire_perimeters_sf |>
  st_intersection(forests_sf) |>
  distinct(year, fire_name, frst_nm) |>
  st_drop_geometry() |>
  as_tibble() |>
  mutate(fire_name = str_to_title(fire_name))

# loop through these and summarize occupancy
fires_forests <- unique(fires_forests_df$fire_name)

y <- list()

unique(occ_df$fire_name)

for (i in seq_along(fires_forests)) {

  # only keep cells for the forest a fire burned in
  forests <-
    fires_forests_df |>
    filter(fire_name == fires_forests[i]) |>
    pull(frst_nm) |>
    as.character()

  print(forests)

  # take the occupancy data, keep those from forests associated with focal fire
  # then join in % burned by that fire
  y[[i]] <-
    occ_df |>
    filter(forest_name %in% str_remove_all(forests, ' National Forest')) |>
    left_join(burned_cells_df |> select(cell_id, category)) |>
    mutate(
      category = case_when(
        is.na(category) ~ 'No high-severity fire',
        TRUE ~ category
      )
    ) |>
    filter(is.na(fire_name) | fire_name == fires_forests[i]) |>
    group_by(posterior_sample, year, category) |>
    summarise(psi = mean(psi)) |>
    group_by(year, category) |>
    summarise(
      occ = mean(psi),
      lower = quantile(psi, 0.025),
      upper = quantile(psi, 0.975)
    ) |>
    ungroup() |>
    mutate(
      category = fct_relevel(category, 'No high-severity fire', '<0.5', '\u22650.5'),
      fire = fires_forests[i]
    )

}

design <- matrix(c(1, 2, 3, NA, NA, 5, 6, 7, 8, 9, 10, NA, NA, NA, NA), 3, 5, byrow = TRUE)

fire_order <-
  fire_size_df |>
  arrange(year, desc(size)) |>
  pull(fire_name)

fire_occ_plot <-
  y |>
  bind_rows() |>
  left_join(fires_forests_df |> select(fire = fire_name, fire_year = year) |> distinct()) |>
  mutate(
    rect_min = fire_year - 0.15,
    rect_max = case_when(
      fire_year == 2020 ~ fire_year + 0.65,
      TRUE ~ fire_year + 0.15
    )
  ) |>
  mutate(
    fire = case_when(
      fire == 'Knp Complex' ~ 'KNP Complex',
      TRUE ~ fire
    )
  ) |>
  mutate(fire = fct_relevel(fire, fire_order)) |>
  ggplot() +
  geom_rect(aes(xmin = rect_min, xmax = rect_max, ymin = -Inf, ymax = Inf), fill = 'red', color = NA, alpha = 0.025) +
  # annotate("rect", xmin = 2021 - 0.15, xmax = 2021 + 0.15, ymin=-Inf, ymax=Inf, alpha=0.2, fill="red") +
  geom_line(
    aes(
      year,
      occ,
      group = category,
      color = category
    ),
    linetype = 2,
    linewidth = 4/10,
    alpha = 7/10,
    position = position_dodge2(width = 0.25),
    show.legend = FALSE
  ) +
  geom_pointrange(
    aes(
      x = year,
      y = occ,
      ymin = lower,
      ymax = upper,
      color = category
    ),
    position = position_dodge2(width = 0.25),
    fatten = 5/10,
    alpha = 10/10
  ) +
  theme_light(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 5.5,  # Top margin
      r = 5.5 * 6,  # Right margin
      b = -5,  # Bottom margin
      l = 5.5,  # Left margin
      unit = "pt"
    ),
    legend.position = c(0.83, 0.895)
  ) +
  # scale_color_manual(values = c('dodgerblue3', 'darkorange', '#DC3220'), name = 'Proportion burned at high-severity', guide = guide_legend(alpha = 1)) +
  scale_color_viridis_d(option = 'inferno', end = 0.6, direction = 1, name = 'Proportion burned at high-severity', guide = guide_legend(alpha = 1)) +
  # ggokabeito::scale_color_okabe_ito(order = c(3, 1, 6), name = 'Proportion burned at high-severity', guide = guide_legend(alpha = 1)) +
  labs(
    x = 'Year of passive acoustic monitoring',
    y = expression(psi)
  ) +
  # ylim(c(0, 1)) +
  facet_manual(~fire, design) +
  guides(colour = guide_legend(override.aes = list(size = 3/10))) +
  theme(axis.title.x = element_text(vjust = 38)) +
  coord_cartesian(xlim = c(2020.75, 2024.25))

# final plot
ggdraw(fire_occ_plot) +
  annotate(
    "rect",
    xmin = c(0.58943),
    ymin = c(0.6874),
    xmax = c(0.62),
    ymax = c(0.9259),
    fill = c("grey70")
  ) +
  annotate(
    "rect",
    xmin = c(0.58943 - 0.3535),
    ymin = c(0.6874 - 0.622),
    xmax = c(0.62 - 0.3535),
    ymax = c(0.9259 - 0.622),
    fill = c("grey70")
  ) +
  annotate(
    "rect",
    xmin = c(0.58943 + 0.353),
    ymin = c(0.6874 - (0.622 / 2)),
    xmax = c(0.62 + 0.353),
    ymax = c(0.9259 - (0.622 / 2)),
    fill = c("grey70")
  ) +
  draw_label("2020", x = 0.605, y = 0.5 + 0.31, angle = 270, size = 8, color = 'white') +
  draw_label("2021", x = 0.9575, y = 0.5, angle = 270, size = 8, color = 'white') +
  draw_label("2022", x = 0.2525, y = 0.5 - 0.31, angle = 270, size = 8, color = 'white') +
  annotate("text", x = 0.13, y = 0.85, label = "2020 fire season", size = 2, fontface = 'italic') +
  annotate("text", x = 0.17, y = (0.85 / 1.55) + 0.006, label = "Fire seasons occur after", size = 2, fontface = 'italic') +
  annotate("text", x = 0.159, y = (0.85 / 1.625) + 0.006, label = "acoustic monitoring", size = 2, fontface = 'italic') +
  annotate(
    geom = "curve", x = 0.09, y = 0.85, xend = 0.074, yend = 0.85, size = 0.5 / 2,
    curvature = .1, arrow = arrow(length = unit(0.8, "mm"))
  ) +
  annotate(
    geom = "curve", x = 0.09 + 0.023, y = 0.535 + 0.006, xend = 0.074 + 0.023, yend = 0.535 + 0.006, size = 0.5 / 2,
    curvature = .1, arrow = arrow(length = unit(0.8, "mm"))
  )
# draw_label("Year of fire ignition", x = 0.945, y = 0.821 - 0.2825, angle = 270, size = 10, color = 'black')
ggsave(here::here('figures/occupancy_by_fire.png'), width = 8, height = 4, units = 'in', dpi = 600)


# study area figure -------------------------------------------------------

focal_hexes_sf <-
  cb_get_spatial('sierra_hexes') |>
  filter(cell_id %in% umf@siteCovs$cell_id)

r_init <- elevation_30s("USA", path = here::here('data'))

# For better handling we set here the names
names(r_init) <- "elev"

ca <-
  rnaturalearth::ne_states('united states of america', returnclass = 'sf') |>
  filter(postal == 'CA')

r <- terra::crop(r_init, ca, mask = TRUE)

slope <- terrain(r, "slope", unit = "radians")
aspect <- terrain(r, "aspect", unit = "radians")
hill <- shade(slope, aspect, 30, 45)

# normalize names
names(hill) <- "shades"

# Hillshading palette
pal_greys <- hcl.colors(1000, "Grays")

# Index of color by cell
index <-
  hill %>%
  mutate(index_col = rescale(shades, to = c(1, length(pal_greys)))) %>%
  mutate(index_col = round(index_col)) %>%
  pull(index_col)

# Get cols
vector_cols <- pal_greys[index_col]

# Base hill plot
hill_plot <-
  ggplot() +
  geom_spatraster(
    data = hill,
    fill = vector_cols,
    maxcell = Inf,
    alpha = 1
  )

# Overlaying and theming
# Aware of limits of the raster
alt_limits <-
  minmax(r) %>%
  as.vector()

# Roud to lower and higher 500 integer with a min of 0
alt_limits <-
  pmax(
    c(floor(alt_limits[1] / 500), ceiling(alt_limits[2] / 500)) * 500,
    -105
  )

base_text_size <- 9

forests_sf <-
  forests_sf |>
  mutate(frst_nm = fct_relevel(frst_nm, forests_north_south))

main <-
  # hillshade
  hill_plot +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = r, maxcell = Inf) +
  scale_fill_hypso_tint_c(
    palette = "wiki-schwarzwald-cont",
    limits = alt_limits,
    alpha = 0.4,
    labels = label_comma(),
    name = 'Elevation (m)',
    guide = guide_colourbar(alpha = 7/10)
  ) +
  # Overlay the_plain
  geom_spatvector(
    data = ca,
    color = alpha("dimgrey", 0.7),
    fill = NA,
    linewidth = 0.15
  ) +
  # forests
  ggnewscale::new_scale_fill() +
  geom_spatvector(
    data = forests_sf,
    aes(
      fill = frst_nm
    ),
    linewidth = 9/10,
    color = NA,
    alpha = 6.5/10,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#F7F7F7", "#737373", "#F7F7F7", "#737373", "#F7F7F7", "#737373", "#F7F7F7")) +
  # study area boundary
  geom_spatvector(
    data = study_area_sf,
    color = alpha("violetred", 1),
    fill = NA,
    linewidth = 1.5/10
  ) +
  # cells monitored
  geom_spatvector(
    data = focal_hexes_sf,
    color = alpha("black", 1),
    fill = "black",
    linewidth = 1/10
  ) +
  theme_void() +
  theme(
    legend.position = c(0.57, 0.85),
    legend.key.size = unit(0.25,"cm"),
    legend.key.width = unit(3/10,"cm"),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 6)
  ) +
  # spatial-aware automagic scale bar
  annotation_scale(
    location = "bl",
    height = unit(0.25/2, "cm"),
    text_cex = 0.4,
    width_hint = 5/10
  ) +
  annotation_north_arrow(
    location = 'bl',
    style = north_arrow_fancy_orienteering(
      text_size = 5,
      line_width = 6/10
    ),
    which_north = "true",
    pad_y = unit(0.18, "in"),
    height = unit(0.55, "cm"),
    width = unit(0.5, "cm")
  ) +
  # forest labels
  geom_sf_label_repel(
    data = forests_sf |> mutate(frst_nm = str_remove(frst_nm, ' National Forest')) |> st_centroid(),
    aes(label = frst_nm),
    force = 100,
    nudge_x = -1.75,
    size = 1.25,
    seed = 10,
    label.padding = 0.2,
    segment.size = 0.5 / 3,
    alpha = 9/10
  )

# inset map
# ca boudning box
ca_bbox_sf <-
  st_bbox(ca) |>
  st_as_sfc() |>
  st_as_sf() |>
  st_transform("+proj=laea +y_0=0 +lon_0=-77 +lat_0=39 +ellps=WGS84 +no_defs")

# world map
world_sf <-
  rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf') |>
  st_transform("+proj=laea +y_0=0 +lon_0=-77 +lat_0=39 +ellps=WGS84 +no_defs") |>
  st_make_valid() |>
  summarise()

# inset map together
world_map <-
  ggplot() +
  geom_sf(data = world_sf, fill = 'grey70', color = NA) +
  geom_sf(data = ca, fill = '#3c3c3c', color = NA) +
  geom_sf(data = ca_bbox_sf, fill = NA, color = 'tomato', linewidth = 1/10) +
  coord_sf(ndiscr = 1000) +
  theme_light() +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = 'transparent'), #transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_line(linewidth = 0.2)
  )

# combine with main map
main_world <-
  main %>%
  ggdraw() +
  draw_plot(
    {
      world_map
    },
    x = 0.06,
    y = 0.12,
    width = 3/10,
    height = 3/10
  )

# aru inset
# example arus
example_arus <-
  tibble(
    aru = c('used', 'unused'),
    utm_zone = c(10, 10),
    utme = c(739807, 738753+200),
    utmn = c(4302051, 4301426-400)
  ) |>
  st_as_sf(coords = c('utme', 'utmn'), crs = 26910)

# buffer by 250 m
example_aru_buffers <-
  example_arus |>
  st_buffer(dist = 250)

# create cell inset map
focal_cell_map <-
  focal_hexes_sf |>
  filter(cell_id == 'C3273') |>
  ggplot() +
  geom_sf(fill = 'white', color = 'black', linewidth = 3/10) +
  geom_sf(data = example_aru_buffers, fill = 'gray90', linetype = 2) +
  geom_sf(data = example_arus, size = 8/10) +
  theme_void() +
  guides(color = 'none') +
  coord_sf(crs = "+proj=omerc +lat_0=36.934 +lonc=-90.849 +alpha=0 +k_0=.7 +datum=WGS84 +units=m +no_defs +gamma=-81")

# combine everything in one map
main_world_cell <-
  main_world %>%
  ggdraw() +
  draw_plot(
    {
      focal_cell_map
    },
    x = 6.25/10,
    y = 5.1/10,
    width = 3/10,
    height = 3/10
  )

# do some annotation
main_world_cell_annotated <-
  main_world_cell +
  annotate("text", x = .8+.05, y = .865, label = "Hexagonal\ngrid cell", size = 1.5) +
  annotate("text", x = .68+.05, y = .85, label = "ARU", size = 1.5) +
  annotate("text", x = .85, y = .45, label = "Buffer distance\nbetween ARUs", size = 1.5) +
  annotate("segment", x = .46+.04, xend = .8406, y = .66-.17, yend = .5406, colour = "black", linewidth = 2/10, linetype = 2) +
  annotate("segment", x = .46+.04, xend = .706, y = .66-.17, yend = .78, colour = "black", linewidth = 2/10, linetype = 2) +
  annotate("text", x = .8+0.0125, y = .63, label = "|", size = 4.1, angle = 270) +
  annotate("text", x = .8+0.007, y = .645, label = "500 m", size = 1.25) +
  annotate(
    'curve',
    x = .855, # Play around with the coordinates until you're satisfied
    y = .48,
    yend = .576,
    xend = .83,
    linewidth = 2/10,
    curvature = 0.25,
    arrow = arrow(length = unit(0.05, 'cm'))
  ) +
  annotate(
    'curve',
    x = .85, # Play around with the coordinates until you're satisfied
    y = .83,
    yend = .78,
    xend = .83,
    linewidth = 2/10,
    curvature = 0,
    arrow = arrow(length = unit(0.05, 'cm'))
  ) +
  annotate(
    'curve',
    x = .73, # Play around with the coordinates until you're satisfied
    y = .83,
    yend = .739,
    xend = .7545,
    linewidth = 2/10,
    curvature = 0,
    arrow = arrow(length = unit(0.05, 'cm'))
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill='white'), #transparent panel bg
    plot.background = element_rect(fill='white', color = NA)
  )
ggsave(plot = main_world_cell_annotated, here::here('figures/map_figure.png'), width = 3, height = 3, units = 'in', dpi = 600)


# posterior summary table -------------------------------------------------

# summarize and format for manuscript table
tmp <-
  summary(fit_re@stanfit, probs = c(0.025, 0.975))$summary |>
  as_tibble(rownames = 'variable') |>
  select(variable, mean, sd, `2.5%`:Rhat) |>
  filter(!str_detect(variable, 'log_lik|lp_')) %>%
  mutate(across(c(n_eff), ~format(round(.x, 0), nsmall = 0))) %>%
  mutate(across(c(Rhat), ~format(round(.x, 2), nsmall = 2))) %>%
  mutate(across(c(mean:`97.5%`), ~format(round(.x, 3), nsmall = 3))) |>
  mutate(
    submodel = str_extract(variable, 'state|col|ext|det'),
    submodel = case_when(
      submodel == 'state' ~ 'Initial occupancy',
      submodel == 'col' ~ 'Colonization',
      submodel == 'ext' ~ 'Extinction',
      TRUE ~ 'Detection'
    )
  ) |>
  select(
    Submodel = submodel,
    Parameter = variable,
    Mean = mean,
    SD = sd,
    `2.5%`,
    `97.5%`,
    Rhat,
    ESS = n_eff
  ) |>
  mutate(
    Parameter = str_remove_all(
      Parameter,
      'beta_state|beta_col|beta_ext|beta_det|b_state|b_col|b_ext|b_det|sigma_state|sigma_col|sigma_ext|sigma_det'
    )
  ) |>
  filter(!str_detect(Parameter, 'forest_name:')) |>
  mutate(
    Parameter = case_when(
      Parameter == '[(Intercept)]' ~ 'Intercept',
      str_detect(Parameter, 'year2022') ~ 'Year (2022)',
      str_detect(Parameter, 'year2023') ~ 'Year (2023)',
      str_detect(Parameter, 'year2024') ~ 'Year (2024)',
      Parameter == '[p_high_severity_scaled]' ~ 'High severity burn',
      Parameter == '[scale(p_high_severity)]' ~ 'High severity burn',
      Parameter == '[date_scaled]' ~ 'Date',
      Parameter == '[date_scaled2]' ~ 'Date2',
      Parameter == '[log_survey_hours_scaled]' ~ 'log(ARU hours)',
      Parameter == '[elev_lat_resid_scaled]' ~ 'Elevation',
      Parameter == '[elev_lat_resid_scaled2]' ~ 'Elevation2',
      Parameter == '[y_coord_scaled]' ~ 'Latitude',
      Parameter == '[y_coord_scaled2]' ~ 'Latitude2',
      Parameter == '[cfo_ch_scaled]' ~ 'Canopy height',
      Parameter == '[sigma [1|forest_name]]' ~ 'sigma',
      TRUE ~ Parameter
    )
  )

# save
tmp |>
  write_rds(here::here('outputs/posterior_summary.rds'))


# probability of direction ------------------------------------------------

# PD of parameters
prob_direction_df <-
  posteriors_df |>
  mutate(
    sign = sign(value),
    direction = case_when(
      sign == -1 ~ 'negative',
      TRUE ~ 'positive'
    )
  ) |>
  group_by(submodel, parameter, direction) |>
  tally() |>
  mutate(prob_direction = (n / sum(n))) |>
  filter(prob_direction == max(prob_direction)) |>
  select(-n) |>
  ungroup()

# add to parameter summary
parameter_stats_summary_df <-
  posteriors_df |>
  group_by(submodel, parameter) |>
  summarise(
    mean = round(mean(value), 3),
    lower = round(quantile(value, 0.025), 3),
    upper = round(quantile(value, 0.975), 3)
  ) |>
  ungroup() |>
  left_join(prob_direction_df) |>
  arrange(submodel, abs(mean))

# save
parameter_stats_summary_df |>
  write_rds(here::here('outputs/parameter_stats_summary_df.rds'))


# diagnostics of top model ------------------------------------------------

# traceplot
traceplot(fit_re, pars = c("beta_state", "sigma_state", "beta_col", "sigma_col", "beta_ext", "sigma_ext", "beta_det", "sigma_det"), ncol = 3)
ggsave(here::here('figures/traceplot.png'), width = 8.5, height = 10, units = 'in', dpi = 600)

# residual plots (95% of points should fall within shaded area)
det_residuals_plot <- plot_residuals(fit_re, submodel = "det") + ggplot2::ggtitle('b) Detection submodel')
occ_residuals_plot <- plot_residuals(fit_re, submodel = "state") + ggplot2::ggtitle('a) State submodel')

# goodness-of-fit test
# The proportion of draws for which the simulated statistic is larger than the actual statistic should be near 0.5 if the model fits well
fit_re_gof <- read_rds(here::here('outputs/gof_test.rds'))

# plot gof
ppval <- data.frame(lab=paste("P =", round(fit_re_gof@post_pred_p, 3)))

plot_theme <- function(){
  theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.background=element_blank(),
          strip.text=element_blank()
    )
}

gof_plot <-
  ggplot(fit_re_gof@samples, aes(x=.data[["obs"]], y=.data[["sim"]])) +
  geom_abline(aes(intercept=0, slope=1), linewidth=1.2, col='red') +
  geom_point(alpha=0.4) +
  geom_label(
    data=ppval,
    aes(
      x=-Inf,
      y=Inf,
      label=.data[["lab"]]
    ),
    hjust=-0.2,
    vjust=1.4,
    size=5,
    fill='white',
    label.size=0,
    label.padding=unit(0.1, "lines")
  ) +
  ggtitle(paste("Posterior predictive check:", fit_re_gof@statistic)) +
  labs(y="Simulated data", x="Observed data") +
  plot_theme() +
  theme(
    strip.background=element_rect(fill="transparent"),
    strip.text=element_text(size=14)
  ) +
  ggtitle('c) Posterior predictive check')

# all together
diagnostics_plot <- (occ_residuals_plot / det_residuals_plot) | gof_plot
ggplot2::ggsave(plot = diagnostics_plot, filename = here::here('figures/diagnostics_plot.png'), width = 12, height = 8, units = 'in', dpi = 600)


# sampling effort plot ----------------------------------------------------

cell_shp <- cb_get_spatial('sierra_hexes')
df <- read_rds(here::here('data/dynamic_encounter_effort_dfs.rds'))$encounters_long

p <-
  df |>
  select(
    cell_id,
    survey_year,
    occasion,
    survey_hours,
    arus,
    mean_date
  ) |>
  mutate(
    across(c(survey_hours, arus), \(x) na_if(x, 0)),
    mean_date = case_when(
      is.na(survey_hours) ~ NA,
      TRUE ~ mean_date
    )
  ) |>
  left_join(
    cell_shp |>
      st_drop_geometry() |>
      select(cell_id, y_coord)
  ) |>
  mutate(cell_id = fct_reorder(cell_id, y_coord)) |>
  ggplot() +
  geom_tile(aes(occasion, cell_id, fill = survey_hours, color = survey_hours)) +
  scale_fill_viridis_c(name = 'ARU hours', na.value = 'grey80', option = 'plasma') +
  scale_color_viridis_c(name = 'ARU hours', na.value = 'grey80', option = 'plasma') +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5
  ),
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5
  )
  ) +
  facet_wrap(~survey_year, ncol = 4) +
  theme_light() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'top'
  ) +
  labs(
    x = 'Survey week',
    y = 'Southern cells                                                                                               Northern cells'
  )
p

arrowdf <- tibble(survey_year = 2021)

bottom <-
  df |>
  select(
    cell_id,
    survey_year,
    occasion,
    survey_hours,
    arus,
    mean_date
  ) |>
  mutate(
    across(c(survey_hours, arus), \(x) na_if(x, 0)),
    mean_date = case_when(
      is.na(survey_hours) ~ NA,
      TRUE ~ mean_date
    )
  ) |>
  left_join(
    cell_shp |>
      st_drop_geometry() |>
      select(cell_id, y_coord)
  ) |>
  mutate(cell_id = fct_reorder(cell_id, y_coord)) |>
  distinct(cell_id) |>
  arrange(cell_id) |>
  slice(175) |>
  mutate(cell_id = as.character(cell_id)) |>
  pull()

top <-
  df |>
  select(
    cell_id,
    survey_year,
    occasion,
    survey_hours,
    arus,
    mean_date
  ) |>
  mutate(
    across(c(survey_hours, arus), \(x) na_if(x, 0)),
    mean_date = case_when(
      is.na(survey_hours) ~ NA,
      TRUE ~ mean_date
    )
  ) |>
  left_join(
    cell_shp |>
      st_drop_geometry() |>
      select(cell_id, y_coord)
  ) |>
  mutate(cell_id = fct_reorder(cell_id, y_coord)) |>
  distinct(cell_id) |>
  arrange(desc(cell_id)) |>
  slice(175) |>
  mutate(cell_id = as.character(cell_id)) |>
  pull()

effort_plot <-
  p +
  coord_cartesian(xlim = c(0, 19), clip = 'off') +
  geom_segment(
    data = arrowdf,
    aes(x = -3.3, xend = -3.3, y = top, yend = bottom),
    colour = "black",
    linewidth = 4/10,
    alpha = 0.9,
    arrow = arrow(length = unit(0.05, "npc"), ends = 'both')
  )
ggsave(plot = effort_plot, here::here('figures/effort_plot.png'), width = 6, height = 8, units = 'in')
