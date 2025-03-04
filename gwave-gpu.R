devtools::install_github("docopt/docopt.R")
# devtools::install(quick = TRUE, upgrade = 'never')
devtools::install_git("https://gitlab.com/apis-staging/greenwave.git", ref = "neon")
library(tidyverse)
library(lubridate)
library(greenwave)

seed <- 101  # controls the data partition

n_burnin <- 30000
n_iter <- 15000
n_chains <- 4  # number of MCMC chains to run

cache <- "/datadrive/neon"#"cache"

re_prior_dist = 'cauchy'  # or 'flat' (the default) or 'cauchy'
reg_on_sigma = 'none'

gcc_re_sf_rnd1 <-
  c(
    list(
      phi = lapply(greenwave:::get_re_sf_default()[c('alpha1')], function(x) x * 30)
    ),
    list(
      theta = lapply(greenwave:::get_re_sf_default()[c('alpha1', 'omega2')], function(x) x * 5)
    )
  )
gcc_re_sf_rnd2 <- lapply(greenwave:::set_re_sf(re_sf = gcc_re_sf_rnd1)[c('phi', 'theta')], function(z) {
  lapply(z, function(w) w * 3)
})
gcc_priors <-
  "list(alpha1 = list(mean = 115, sd = 50),
       delta1 = list(mean = log(1/30), sd = greenwave:::get_sd_delta(30, shift = 15, log = TRUE)),
       delta2 = list(mean = log(1/30), sd = greenwave:::get_sd_delta(30, shift = 15, log = TRUE)),
       gamma1_R = list(mean = 0.34, sd = 0.1),
       gamma2 = list(mean = 0.42, sd = 0.1),
       lambda = list(mean = log(235-115), sd = 0.5),
       omega2 = list(mean = -3e-04, sd = 0.1))"  # plot(density(exp(rnorm(1E4, log(275-200), 0.5))))
gcc_re_priors <-
  list(phi = list(alpha1 = 120, delta1 = 4, delta2 = 8, 
                  gamma1_R = 0.05, gamma2 = 0.03, lambda = 4, omega2 = 4e-04),
       theta = list(alpha1 = 40, delta1 = 4, delta2 = 4, 
                    gamma1_R = 0.04, gamma2 = 0.06, lambda = 8, omega2 = 0.002))

# greenwave:::set_re_priors(re_prior_dist, 
#                           re_priors = gcc_re_priors, 
#                           re_sf = greenwave:::set_re_sf(re_sf = gcc_re_sf_rnd2))



opts <- gw_opts(
  seed = seed,
  # n_locs = n_locs, n_years = n_years,
  use_single_spacetime = FALSE,
  excl_params = 'c("omega1_R")',  # "omega2"
  n_burnin = n_burnin, n_iter = n_iter, n_chains = n_chains,
  resp_data = '/home/greenwave/phenology-targets.csv.gz', cov_data = 'none',  # file.path(assets_path, eo_data)
  cache = cache, omit_sys_time = TRUE,
  loc_attr = "siteID", dttm_attr = 'time',
  p_parts = c(training = 0.80, test = 0.20),
  re_sf = gcc_re_sf_rnd2,
  reg_on_sigma = reg_on_sigma,
  re_prior_dist = re_prior_dist,
  re_priors = gcc_re_priors,
  # Lmin = 10, Lmax = 15,
  priors = gcc_priors
)
print(opts)

# Load response and covariate data.
which_vi <- "gcc_90"

resp_data <- opts$resp_data %>%
  drop_na(gcc_90, date) %>%
  # filter(date <= ymd('2020-07-02')) %>%
  # rename(location = site) %>%
  gw_partition_vi(opts$cov_data, opts, vi = !!which_vi)

Y_in_sample <- resp_data %>% filter(partition == 'training')  # dim(Y_in_sample)

obs <- gw_vis_obs(resp_data, cache_dir = opts$cache, file_ext = 'png')
# obs[[1]]

draws_ckpt <- file.path(opts$cache, 'draws.ckpt')
draws_extra_ckpt <- file.path(opts$cache, 'draws-extra.ckpt')
draws <- if (file.exists(draws_ckpt)) {
  if (file.exists(draws_extra_ckpt)) readRDS(draws_extra_ckpt) else readRDS(draws_ckpt)
} else {
  gw_fit(Y_in_sample, Xs = NULL, opts)
}

traces <- gw_vis_traces(draws, coefs = 'fixed', cache_dir = opts$cache, file_ext = 'png')


sf_eval <- gw_param_scaled_dists(draws, cache_dir = opts$cache, file_ext = 'png')
gw_param_scaled_dists(draws, cache_dir = opts$cache, scaled = FALSE, file_ext = 'png')

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

dqs_data <- file.path(opts$cache, 'dqs', 'dqs.data')
dqs <- if (file.exists(dqs_data)) {
  readRDS(dqs_data)
} else {
  gw_dqs(draws, other_params = greenwave:::param_names(), cache_dir = opts$cache, n_samples = 250,
         batch_size = 50)
  # gw_dqs(draws, other_params = greenwave:::param_names(), cache_dir = opts$cache)
}

dqs_tidy <- gw_dqs_tidy(dqs)
dq_just_params <- dqs_tidy$fixed_plus_random %>%
  filter(dq %in% c(greenwave:::param_names(), 'alpha2')) %>%
  gw_vis_dqs(cache_dir = opts$cache, suffix = 'just-params', file_ext = 'png', n_cols = 1)
# focal_dqs <- c('sos', 'accum_vi_over_min', 'n_days_over_thresh')
# dq_other_errbars <- dqs_tidy$fixed_plus_random %>%
#   filter(dq %in% focal_dqs) %>%
#   # filter(!dq %in% c(greenwave:::param_names(), 'alpha2')) %>%
#   gw_vis_dqs(cache_dir = opts$cache, suffix = 'extra', file_ext = 'png', n_cols = 1)
