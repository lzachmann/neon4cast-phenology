# devtools::install_github("docopt/docopt.R")
# devtools::install_git("https://gitlab.com/apis-staging/greenwave.git", ref = "neon")
library(tidyverse)
library(lubridate)
library(greenwave)

seed <- 101  # controls the data partition

n_burnin <- 7500
n_iter <- 2500
n_chains <- 2  # number of MCMC chains to run

cache <- "output"

re_prior_dist = 'flat'  # or 'flat' (the default) or 'cauchy'
reg_on_sigma = 'none'

gcc_re_sf <-
  "c(list(
    kappa = lapply(greenwave:::get_re_sf_default()[c('alpha1')], function(x) x * 40)
  ))"
gcc_priors <-
  "list(alpha1 = list(mean = 79, sd = 50),
       delta1 = list(mean = log(1/20), sd = greenwave:::get_sd_delta(20, shift = 15, log = TRUE)),
       delta2 = list(mean = log(1/20), sd = greenwave:::get_sd_delta(20, shift = 15, log = TRUE)),
       gamma1_R = list(mean = 0.35, sd = 0.1),
       gamma2 = list(mean = 0.45, sd = 0.1),
       lambda = list(mean = log(265-79), sd = 0.5))"  # plot(density(exp(rnorm(1E4, log(275-200), 0.5))))
# gcc_re_priors <- 
#   "list(theta = list(alpha1 = 5, gamma1_R = 0.03), phi = list(omega1_R = 0.01))"


opts <- gw_opts(
  seed = seed,
  # n_locs = n_locs, n_years = n_years,
  use_single_spacetime = TRUE,
  excl_params = 'c("omega1_R")',  # "omega2"
  n_burnin = n_burnin, n_iter = n_iter, n_chains = n_chains,
  resp_data = 'phenology-targets.csv.gz', cov_data = 'data/neon/field-sites-subset.csv',  # file.path(assets_path, eo_data)
  cache = cache, omit_sys_time = TRUE, 
  loc_attr = "siteID", dttm_attr = 'time',
  p_parts = c(training = 0.25, test = 0.75),
  re_sf = gcc_re_sf,
  reg_on_sigma = reg_on_sigma,
  re_prior_dist = re_prior_dist,
  # re_priors = gcc_re_priors,
  priors = gcc_priors
)
# opts <- readRDS('/datadrive/cache/synth-paper/SH/cffe6090/opts.meta')
print(opts)

# Load response and covariate data.
which_vi <- "gcc_90"

resp_data <- opts$resp_data %>%
  drop_na(gcc_90, date) %>%
  # filter(date <= ymd('2020-07-02')) %>%
  # rename(location = site) %>%
  gw_partition_vi(opts$cov_data, opts, vi = !!which_vi)

Y_in_sample <- resp_data %>% filter(partition == 'training')  # dim(Y_in_sample)

drivers <- opts$cov_data %>% 
  crossing(cal_year = seq(min(Y_in_sample$cal_year) - 1, max(Y_in_sample$cal_year) + 1)) %>% 
  left_join(
    Y_in_sample %>% 
      mutate(field_site_id = as.character(loc_id)) %>% 
      distinct(loc_id, field_site_id)
  )
effects <- list(alpha1 = c("field_mean_elevation_m",
                           "field_mean_annual_temperature_C"),
                delta1 = c("field_mean_annual_temperature_C"),
                delta2 = c("field_mean_annual_precipitation_mm"),
                gamma2 = c("field_mean_annual_temperature_C"),
                lambda = c("field_mean_elevation_m",
                           "field_mean_annual_temperature_C"))
X_in_sample <- gw_format_covariates(
  data = drivers,
  Y = Y_in_sample, 
  loc_id = 'loc_id',
  cal_year = 'cal_year',
  effects = effects,
  to_z_score = TRUE
)

obs <- gw_vis_obs(resp_data, cache_dir = opts$cache, file_ext = 'png')
obs[[1]]

draws_ckpt <- file.path(opts$cache, 'draws.ckpt')
draws_extra_ckpt <- file.path(opts$cache, 'draws-extra.ckpt')
draws <- if (file.exists(draws_ckpt)) {
  if (file.exists(draws_extra_ckpt)) readRDS(draws_extra_ckpt) else readRDS(draws_ckpt)
} else {
  gw_fit(Y_in_sample, Xs = X_in_sample, opts)
}

traces <- gw_vis_traces(draws, coefs = 'fixed', cache_dir = opts$cache, file_ext = 'png')

sf_eval <- gw_param_scaled_dists(draws, cache_dir = opts$cache, file_ext = 'png')

# Model evaluation
ppds_data <- file.path(opts$cache, 'pred', 'ppds.data')
ppds <- if (file.exists(ppds_data)) {
  readRDS(ppds_data)
} else {
  gw_ppd(draws, cache_dir = opts$cache, batch_size = 10)  # calculate once, recycled below!
}

# Y_out_of_sample <- if (length(opts$p_parts) > 1) resp_data else NULL
pred <- gw_vis_fit(ppds, opts, Y = resp_data,
                   n_yhat_draws = 100, cache_dir = opts$cache,
                   file_ext = 'png')
po <- gw_vis_po(ppds, opts, Y = resp_data, cache_dir = opts$cache,
                file_ext = 'png')

# dqs_data <- file.path(opts$cache, 'dqs', 'dqs.data')
# dqs <- if (file.exists(dqs_data)) {
#   readRDS(dqs_data)
# } else {
#   gw_dqs(draws, other_params = greenwave:::param_names(),
#          cache_dir = opts$cache, n_samples = 250, batch_size = 50)
# }

# dqs_tidy <- gw_dqs_tidy(dqs)
# dq_just_params <- dqs_tidy$fixed_plus_random %>%
#   filter(dq %in% c(greenwave:::param_names(), 'alpha2')) %>%
#   gw_vis_dqs(cache_dir = opts$cache, suffix = 'just-params', file_ext = 'png', n_cols = 1)
# focal_dqs <- c('sos', 'accum_vi_over_min', 'n_days_over_thresh')
# dq_other_errbars <- dqs_tidy$fixed_plus_random %>%
#   filter(dq %in% focal_dqs) %>%
#   # filter(!dq %in% c(greenwave:::param_names(), 'alpha2')) %>%
#   gw_vis_dqs(cache_dir = opts$cache, suffix = 'extra', file_ext = 'png', n_cols = 1)

dqs_re_data <- file.path(opts$cache, 'dqs', 'dqs-re.data')
dqs_random_effects <- if (file.exists(dqs_re_data)) {
  readRDS(dqs_re_data)
} else {
  gw_dqs(draws, which_effects = 'random', other_params = greenwave:::param_names(),
         cache_dir = opts$cache, n_samples = 250, batch_size = 50,
         suffix = 're')
}


# field_latitude field_longitude field_mean_elevation_m field_mean_annual_temperature_C field_mean_annual_precipitation_mm
candidates <- grep('mean', names(drivers), value = TRUE) # did try latitude
drivers %>% select(candidates) %>% cor
effects_new <- list(alpha1 = candidates,#grep('^.*_(JFM|AMJ)$', candidates, value = TRUE),
                    lambda = candidates,#grep('^.*_(JAS)$', candidates, value = TRUE),
                    gamma1_R = candidates,#grep('^.*_(OND)$', candidates, value = TRUE),
                    gamma2 = candidates,#grep('^.*_(AMJ|JAS)$', candidates, value = TRUE),
                    delta1 = candidates,#grep('^.*_(AMJ)$', candidates, value = TRUE),
                    delta2 = candidates#grep('^.*_(JAS)$', candidates, value = TRUE)
                    )
# effects_new <- lapply(effects_new, function(x) {
#   unique(c(x, c('elevation', 'northness')))  # 'lat', 'eastness'
# })

Xs_speculative <- gw_format_covariates(
  data = drivers,
  Y = Y_in_sample, 
  loc_id = 'loc_id',
  cal_year = 'cal_year',
  effects = effects_new,
  # end_year = 2021,
  to_z_score = TRUE
)

unexplained <- gw_vis_resid(dqs_random_effects, Xs = Xs_speculative,
                            cache_dir = opts$cache, file_ext = 'png')
