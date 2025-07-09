

# setup -------------------------------------------------------------------

library(tidyverse)
library(CAbioacoustics)
library(fs)
library(ubms)

# if model warnings suggest running model for more iterations change this
iter_factor <- 1

# chains, iter, etc.
n_chains <- 3
n_iter <- 1500 * iter_factor
n_warmup <- n_iter / 2
n_cores <- n_chains


# read in umf -------------------------------------------------------------

umf <-
  read_rds(here::here('data/anonymized_umf.rds'))


# run model ---------------------------------------------------------------

# full model, with random forest effects
fit_re <-
  stan_colext(
    # initial occupancy
    psiformula = ~ (1|forest_name) + # forest random effect
      elev_lat_resid_scaled + elev_lat_resid_scaled2 + # elevation
      y_coord_scaled + y_coord_scaled2 +  # latitude
      cfo_ch_scaled + # canopy height
      p_high_severity_scaled, # fire
    # extinction
    epsilonformula = ~ year + (1|forest_name) + scale(p_high_severity),
    # colonization
    gammaformula = ~ year + (1|forest_name) + scale(p_high_severity),
    # detection
    pformula = ~ year + (1|forest_name) + date_scaled + date_scaled2 + log_survey_hours_scaled,
    umf,
    chains = n_chains,
    iter = n_iter,
    warmup = n_warmup,
    seed = 123,
    cores = n_cores,
    refresh = 0
  )

# identical to this specification (same parameter estimates), but the nesting of I() and scale()
# does not allow for correct prediction plots...
# this can't be run as only scaled y_coord is provided with anonymized umf
# fit_re <-
#   stan_colext(
#     # initial occupancy
#     psiformula = ~ (1|forest_name) + # forest random effect
#       scale(elev_lat_resid) + I(scale(elev_lat_resid)^2) + # elevation
#       scale(y_coord) + I(scale(y_coord)^2) + # latitude
#       scale(cfo_ch) + # canopy height
#       scale(p_high_severity),
#     # extinction
#     epsilonformula = ~ year + (1|forest_name) + scale(p_high_severity),
#     # colonization
#     gammaformula = ~ year + (1|forest_name) + scale(p_high_severity),
#     # detection
#     pformula = ~ year + (1|forest_name) + scale(date) + I(scale(date)^2) + scale(log(survey_hours)),
#     umf,
#     chains = n_chains,
#     iter = n_iter,
#     warmup = n_warmup,
#     # for reproducibility
#     seed = 123,
#     cores = n_cores,
#     refresh = 0
#   )

# view model output
print(fit_re)

# save model
dir_create(here::here('outputs'))

fit_re |>
  write_rds(here::here('outputs/fit_re.rds'))


# model diagnostics -------------------------------------------------------

# check this; good rule of thumb of effective samples
# n_eff > 100 * number of chains (300)

# traceplots for all submodels
traceplot(fit_re, pars = c("beta_state", "beta_det", "beta_col", "beta_ext"))

# residual plots (95% of points should fall within shaded area)
plot_residuals(fit_re, submodel = "det")
plot_residuals(fit_re, submodel = "state")

# goodness-of-fit test
fit_re_gof <- gof(fit_re, quiet = TRUE)

# save to plot/report later
fit_re_gof |>
  write_rds(here::here('outputs/gof_test.rds'))
