## line for remote resources
devtools::install(upgrade = FALSE)

## @knitr load_packages
library(occstanhm)
library(tidyverse)
library(ggthemes)

## @knitr simulation
set.seed(123)
n_units <- 5
n_unit_revisit_mean <- 3
k_samples_mean <- 3
k_subsamples <- 6
n_spp <- 4

sim_data <-
  sim_long_form_occ_3(n_units = n_units,
                      n_unit_revisit_mean = n_unit_revisit_mean,
                      n_spp = n_spp,
                      k_samples_mean = k_samples_mean,
                      k_subsamples = k_subsamples,
                      spp_mu = rep(0.75, n_spp),
                      sigma_psi = matrix(c(1.0,  0.0,  0.0,  0.0,
                                           0.0,  1.0,  1.0, -0.7,
                                           0.0,  1.0,  1.0, -0.7,
                                           0.0, -0.7, -0.7,  1.0),
                                         nrow = n_spp, ncol = n_spp),
                      sigma_p = matrix(c(1.0,  0.0,  0.0,  0.0,
                                         0.0,  1.0,  1.0, -0.7,
                                         0.0,  1.0,  1.0, -0.7,
                                         0.0, -0.7, -0.7,  1.0),
                                       nrow = n_spp, ncol = n_spp))

## @knitr sim_format
sim_data$spp_obs <-
  sim_data$spp_obs |>
  arrange(unit, species, revisit) |>
  mutate(unit = paste("Unit", unit),
         unit_id = as.integer(unit_id),
         species_id = as.integer(factor(species))) |>
  rowid_to_column("index")

## @knitr save_sim_data
saveRDS(sim_data,
        "../extdata/occstanhm_3_sim_data.rds")

## @knitr load_sim_data
if (TRUE) {
  sim_data <-
    system.file("extdata/occstanhm_3_sim_data.rds", package = "occstanhm") |>
    readRDS()
}

## @knitr plot_one
sim_data$spp_obs |>
  filter(unit == "Unit 1") |>
  group_by(unit, species, revisit) |>
  summarize(total_y = sum(y), .groups = "drop") |>
  ggplot(aes(x = revisit, y = species, fill = total_y)) +
  geom_tile() +
  scale_fill_gradient("Samples\nwith\ndetections",
                      low = "white", high = "black") +
  theme_bw() +
  theme(strip.background = element_blank())

## @knitr plot_data
many_plot <-
  sim_data$spp_obs |>
  group_by(unit_id, species, revisit) |>
  summarize(total_y = sum(y), .groups = "drop") |>
  ggplot(aes(x = revisit, y = species, fill = total_y)) +
  geom_tile() +
  scale_fill_gradient("Samples\nwith\ndetections",
                      low = "white", high = "black") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(vars(unit_id),
             labeller = label_bquote("Unit" ~ .(unit_id))) +
  ylab("Species") +
  xlab("Revisit")
print(many_plot)
ggsave("many_plot.jpg", width = 8, height = 6)

## @knitr summarize_stan_data_revisit_species
revist_species_data <-
  sim_data$spp_obs |>
  ungroup() |>
  group_by(unit, species) |>
  mutate(min_idx = min(index) - 1) |>
  ungroup() |>
  group_by(unit, unit_id, species, species_id, revisit) |>
  summarize(z_obs = ifelse(sum(y) > 0, 1, 0),
            n_samples = n(),
            revisit_dat_start = min(index - min_idx),
            revisit_dat_stop = max(index - min_idx),
            .groups = "drop") |>
  rowid_to_column("index") |>
  ungroup()

revist_species_data |>
  print(n = 20)

## @knitr summarize_stan_data_revisit_species_3
revisit_species_at_unit_data  <-
  revist_species_data |>
  group_by(unit, species) |>
  summarize(n_revisits = n(),
            revisit_unit_start = min(index),
            revisit_unit_stop = max(index),
            .groups = "drop") |>
  ungroup()

revisit_species_at_unit_data

## @knitr summarize_stan_data_unit_species
unit_species_data <-
  sim_data$spp_obs |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")

unit_species_data |>
  print(n = 20)

## @knitr x_psi_predictor
x_psi <- model.matrix(~ species - 1, revist_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))
x_psi |> dim()
x_psi |> head()

x_psi_star <-
  sim_data$spp_obs |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))
x_psi_star |> head()
x_psi_star |> dim()

## @knitr check_psi
x_psi_hyper <-
  sim_data$spp_obs |>
  group_by(unit, species) |>
  summarize(n = n(), .groups = "drop") |>
  mutate(one = 1)

## check math for beta coefficient
beta_psi_check <- matrix(1,
                         nrow = sim_data$spp_obs |> distinct(unit) |> nrow(),
                         ncol = x_psi |> ncol())
beta_psi_check |> dim()

logit_psi_check <- matrix(1,
                          nrow = revist_species_data |> nrow(),
                          ncol = 1)

logit_psi_check |> dim()
rowSums(x_psi, beta_psi_check[revist_species_data$unit_id, ]) |> length()

## @knitr beta_priors
## prior
beta_psi_star_prior <-
  matrix(1.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())

## @knitr check_beta_priors
## check math for prior
x_psi_star_check <-
  matrix(1,
         nrow = x_psi_star |> nrow(),
         ncol = x_psi_star |> ncol())

x_psi_star_check %*% beta_psi_star_prior  |> dim()
beta_psi_check |> dim()

## @knitr check_beta_star
beta_psi_star_check <-
  matrix(1,
         nrow = 1,
         ncol = sim_data$spp_obs |> distinct(species) |> nrow())

# lhs
beta_psi_check |> dim()

# rhs elements
x_psi_star |> dim()
beta_psi_star_check |> dim()

# rhs product
x_psi_star %*%
  beta_psi_star_check |>
  dim()

## @knitr theta_coef
# theta-level
w_theta <- model.matrix(~ species - 1, sim_data$spp_obs)
colnames(w_theta) <- gsub("species", "", colnames(w_theta))
w_theta |> head()
w_theta |> dim()

## @knitr theta_coef_check
## check math for delta coefficient
alpha_theta_check <-
  matrix(1, nrow = sim_data$spp_obs |> distinct(unit) |> nrow(),
         ncol = w_theta |> ncol())

alpha_theta_check |> dim()

logit_theta_check <-
  matrix(1, nrow = sim_data$spp_obs |> nrow(), ncol = 1)
logit_theta_check |> dim()

rowSums(w_theta, alpha_theta_check[sim_data$spp_obs$unit_id, ]) |> length()

## @knitr theta_hyper
## hyper parameters
w_theta_star <-
  sim_data$spp_obs |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

w_theta_star |> dim()
w_theta_star |> head()

alpha_theta_star_prior <-
  matrix(1,
         nrow = w_theta_star |> ncol(),
         ncol = w_theta |> ncol())

## @knitr theta_check_priors
w_theta_star_check <-
  matrix(1,
         nrow = w_theta_star |> nrow(),
         ncol = w_theta_star |> ncol())

alpha_theta_star_check <-
  matrix(1,
         nrow = 1,
         ncol = sim_data$spp_obs |> distinct(species) |> nrow())

# lhs
w_theta_star_check |> dim()

# rhs elements
w_theta_star |> dim()
alpha_theta_check |> dim()

# rhs product
w_theta_star %*%
  alpha_theta_star_check |>
  dim()

## @knitr p_coef
# p-level
v_p <- model.matrix(~ species - 1, sim_data$spp_obs)
colnames(v_p) <- gsub("species", "", colnames(v_p))
v_p |> head()
v_p |> dim()

## @knitr p_coef_check
## check math for delta coefficient
delta_p_check <- matrix(1,
                        nrow = sim_data$spp_obs |> distinct(unit) |> nrow(),
                        ncol = v_p |> ncol())

delta_p_check |> dim()

logit_p_check <- matrix(1,
                        nrow = sim_data$spp_obs |> nrow(),
                        ncol = 1)
logit_p_check |> dim()

rowSums(v_p, delta_p_check[sim_data$spp_obs$unit_id, ]) |> length()

## @knitr p_hyper
## hyper parameters
v_p_star <-
  sim_data$spp_obs |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

v_p_star |> dim()
v_p_star |> head()

delta_p_star_prior <-
  matrix(1,
         nrow = v_p_star |> ncol(),
         ncol = v_p |> ncol())

## @knitr p_check_priors
v_p_star_check <-
  matrix(1,
         nrow = v_p_star |> nrow(),
         ncol = v_p_star |> ncol())

delta_p_star_check <-
  matrix(1,
         nrow = 1,
         ncol = sim_data$spp_obs |> distinct(species) |> nrow())
# lhs
delta_p_check |> dim()

# rhs elements
v_p_star |> dim()
delta_p_check |> dim()

# rhs product
v_p_star %*%
  beta_psi_star_check |>
  dim()

## @knitr save_as_stan
# first are observed data
stan_data <-
  list(y = sim_data$spp_obs |> pull(y),
       a_obs = sim_data$spp_obs |> pull(a_obs),
       k_subsamples = sim_data$spp_obs |> pull(k_subsample),
       z_obs = revist_species_data |> pull(z_obs),
       # next are the dimensions for data (including vectors)
       n_unit_species = unit_species_data |> nrow(),
       n_total_revisits = revist_species_data |> nrow(),
       n_total_samples = sim_data$spp_obs |> nrow(),
       n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
       n_samples = revist_species_data |> pull(n_samples),
       m_beta = x_psi |> ncol(),
       n_beta = sim_data$spp_obs |> distinct(unit) |> nrow(),
       m_beta_star = x_psi_star |> ncol(),
       n_beta_star = x_psi_star |> nrow(),
       m_alpha = w_theta |> ncol(),
       n_alpha = sim_data$spp_obs |> distinct(unit_id) |> nrow(),
       m_alpha_star = w_theta_star |> ncol(),
       n_alpha_star = w_theta_star |> nrow(),
       m_delta = v_p |> ncol(),
       n_delta = sim_data$spp_obs |> distinct(unit_id) |> nrow(),
       m_delta_star = v_p_star |> ncol(),
       n_delta_star = v_p_star |> nrow(),
       # indexing
       unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
       unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),
       revisit_dat_start = revist_species_data |> pull(revisit_dat_start),
       revisit_dat_stop = revist_species_data |> pull(revisit_dat_stop),
       revisit_unit_start = revisit_species_at_unit_data |>
         pull(revisit_unit_start),
       revisit_unit_stop = revisit_species_at_unit_data |>
         pull(revisit_unit_stop),
       # psi-level inputs
       x_psi = x_psi,
       x_psi_star = x_psi_star,
       # psi-level inputs
       w_theta = w_theta,
       w_theta_star = w_theta_star,
       # p-level inputs
       v_p = v_p,
       v_p_star = v_p_star,
       # psi-level priors
       beta_psi_star_prior = beta_psi_star_prior,
       beta_psi_star_sd_prior = 0.1,
       eta_psi = 1.0,
       jj_psi = revist_species_data |> dplyr::pull(unit_id),
       # theta-level priors
       alpha_theta_star_prior = alpha_theta_star_prior,
       alpha_theta_star_sd_prior = 0.1,
       eta_theta = 1.0,
       jj_theta = sim_data$spp_obs |> dplyr::pull(unit_id),
       # p-level priors
       delta_p_star_prior = delta_p_star_prior,
       delta_p_star_sd_prior = 0.1,
       eta_p = 1.0,
       jj_p = sim_data$spp_obs |> dplyr::pull(unit_id),
       # priors for sigams
       sigma_psi_prior = 2.5,
       sigma_psi_prior_sd = 0.05,
       sigma_theta_prior = 1.75,
       sigma_theta_prior_sd = 0.05,
       sigma_p_prior = 1.75,
       sigma_p_prior_sd = 0.05,
       # setting for reduce_sum
       grainsize = 1,
       ## Store data (not used in Stan)
       dat = sim_data$spp_obs,
       unit_species_data = unit_species_data)

## @knitr run_stan_model
n_chains <- 3
n_parallel_chains <- 3
n_threads_per_chain <- 5
n_refresh <- 100


# for sampling numbers
n_warmup <- 200 # 10000
n_sample <- 200 # 10000

fit <- occstanhm_3(stan_data,
                   n_chains = n_chains,
                   n_parallel_chains = n_parallel_chains,
                   n_threads_per_chain = n_threads_per_chain,
                   n_refresh = n_refresh,
                   n_warmup = n_warmup,
                   n_sample = n_sample)

## @knitr fit_summary
fit_summary <-
  fit$summary(.cores = 15)

## @knitr save_stan_outputs
saveRDS(stan_data,
        "../extdata/occstanhm_3_stan_data.rds")

saveRDS(fit,
        "../extdata/occstanhm_3_fit.rds")

saveRDS(fit_summary,
        "../extdata/occstanhm_3_fit_summary.rds",
        compress = "xz")


## @knitr load_stan_outputs
if (TRUE) {
  fit <-
    system.file("extdata/occstanhm_3_fit.rds", package = "occstanhm") |>
    readRDS()

  fit_summary <-
    system.file("extdata/occstanhm_3_fit_summary.rds", package = "occstanhm") |>
    readRDS()
}

## @knitr model_diag
fit$diagnostic_summary()

## @knitr default_fit_summary
fit_summary

## @knitr default_fit_summary_sort_1
fit_summary |>
  arrange(rhat)

## @knitr default_fit_summary_sort_1
fit_summary |>
  arrange(-rhat)

## @knitr beta_coef
## Look at psi-level model outputs
beta_detail <-
  fit_summary |>
  filter(grepl("beta_psi\\[", variable)) |>
  select(variable, q5, median, q95) |>
  mutate(l95 = plogis(q5),
         u95 = plogis(q95),
         prob = plogis(median)) |>
  select(-q5, -q95, -median) |>
  mutate(unit_id =
         as.integer(gsub("beta_psi\\[(\\d+),(\\d+)\\]", "\\1", variable)),
         species_id = as.integer(gsub("beta_psi\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(stan_data$unit_species_data,
            by = c("species_id", "unit_id")) |>
  select(variable, unit, unit_id, species, species_id, l95, prob, u95, psi_obs)

beta_detail

## @knitr beta_coef_plot
beta_plot <-
  beta_detail |>
  mutate(unit = factor(unit),
         unit = fct_reorder(unit, unit_id)) |>
  ggplot(aes(x = `unit`,
             color = species,
             y = prob, ymin = l95, ymax = u95)) +
  geom_point(position = position_dodge(width = 0.5), aes(size = psi_obs)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  scale_color_colorblind("Species") +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_blank()) +
  ylim(c(0, 1)) +
  scale_size_continuous("Propotion\nof\nrevisits\nwith\ndetections") +
  ylab(expression("estimated" ~ psi)) +
  xlab("Unit")
print(beta_plot)
ggsave("beta_plot.jpg", beta_plot, width = 6, height = 8)

## @knitr extract_alpha_coef
## Look at theta-level model estimate
dat_rough_theta <-
  stan_data$dat |>
  filter(z_obs > 0) |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(theta_rough = mean(a_obs), .groups = "drop")

dat_theta_all <-
  stan_data$dat |>
  distinct(unit, unit_id, species, species_id) |>
  full_join(dat_rough_theta, by = c("unit", "unit_id",
                                    "species",
                                    "species_id"))
theta_details <-
  fit_summary |>
  filter(grepl("alpha_theta\\[", variable)) |>
  mutate(unit_id = as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]", "\\1",
                                   variable)),
         species_id = as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(dat_theta_all, by = c("species_id", "unit_id")) |>
  mutate(prob = plogis(median),
         l95 = plogis(q5),
         u95 = plogis(q95)) |>
  select(variable, unit, unit_id, species, l95, prob, u95, theta_rough) |>
  mutate(theta_rough = ifelse(is.na(theta_rough), 0, theta_rough),
         diff = theta_rough - prob)

theta_details |>
  print(n = Inf)

## @knitr alpha_plot
theta_plot <-
  theta_details |>
  mutate(unit = factor(unit),
         unit = fct_reorder(unit, unit_id)) |>
  ggplot(aes(x = unit,
             color = species,
             y = prob, ymin = l95, ymax = u95)) +
  geom_point(position = position_dodge(width = 0.5), aes(size = theta_rough)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  scale_color_colorblind("Species") +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_blank())  +
  scale_size_continuous("Propotion\nof\npositive\nsamples") +
  ylim(c(0, 1)) +
  ylab(expression("estimated" ~ theta)) +
  xlab("Unit")
print(theta_plot)
ggsave("theta_plot.jpg", theta_plot, width = 6, height = 8)

## @knitr extract_p_coef
## Look at p-level model estimate
dat_rough_p <-
  stan_data$dat |>
  filter(z_obs > 0) |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(p_rough = mean(y / k_subsample), .groups = "drop")

dat_p_all <-
  stan_data$dat |>
  distinct(unit, unit_id, species, species_id) |>
  full_join(dat_rough_p, by = c("unit", "unit_id",
                                "species",
                                "species_id"))
p_details <-
  fit_summary |>
  filter(grepl("delta_p\\[", variable)) |>
  mutate(unit_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]", "\\1",
                                   variable)),
         species_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(dat_p_all, by = c("species_id", "unit_id")) |>
  mutate(prob = plogis(median),
         l95 = plogis(q5),
         u95 = plogis(q95)) |>
  select(variable, unit, unit_id, species, l95, prob, u95, p_rough) |>
  mutate(p_rough = ifelse(is.na(p_rough), 0, p_rough),
         diff = p_rough - prob)

p_details |>
  print(n = Inf)

## @knitr delta_plot
p_plot <-
  p_details |>
  mutate(unit = factor(unit),
         unit = fct_reorder(unit, unit_id)) |>

  ggplot(aes(x = unit,
             color = species,
             y = prob, ymin = l95, ymax = u95)) +
  geom_point(position = position_dodge(width = 0.5), aes(size = p_rough)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  scale_color_colorblind("Species") +
  coord_flip() +
  theme_bw() +
  theme(strip.background = element_blank())  +
  scale_size_continuous("Propotion\nof\npositive\ndetections") +
  ylim(c(0, 1)) +
  ylab(expression("estimated" ~ italic(p))) +
  xlab("Unit")
print(p_plot)
ggsave("p_plot.jpg", p_plot, width = 6, height = 8)

## @knitr est_corr
## Look at estimated correlations from the model
fit_summary |>
  filter(grepl("^Omega_psi\\[", variable)) |>
  select(variable, q5, median, q95)

fit_summary |>
  filter(grepl("^Omega_theta\\[", variable)) |>
  select(variable, q5, median, q95)

fit_summary |>
  filter(grepl("^Omega_p\\[", variable)) |>
  select(variable, q5, median, q95)

## @knitr plot_corr
gsub_pattern <- "^Omega_(\\w+)\\[(\\d+),(\\d+)\\]"

spp_key <-
  stan_data$unit_species_data |>
  distinct(species, species_id) |>
  mutate(species_row = species,
         species_col = species) |>
  select(-species)

corr_plot_ci <-
  fit_summary |>
  filter(grepl(gsub_pattern, variable)) |>
  mutate(species_row_id = as.integer(gsub(gsub_pattern, "\\2", variable)),
         species_col_id = as.integer(gsub(gsub_pattern, "\\3", variable)),
         level = gsub(gsub_pattern, "\\1", variable)) |>
  full_join(spp_key |> select(species_id, species_row),
            by = c("species_row_id" = "species_id")) |>
  full_join(spp_key |> select(species_id, species_col),
            by = c("species_col_id" = "species_id")) |>
  mutate(level_plot = factor(level,
                             levels = c("psi", "theta", "p"),
                             labels = c("expression(psi)",
                                        "expression(theta)",
                                        "expression(italic(p))"))) |>
  ggplot(aes(x = species_row, y = median, ymin = q5, ymax = q95)) +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  facet_grid(vars(level), vars(species_col), labeller = "label_parsed") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1.15, "lines")) +
  xlab("Species") +
  ylab("Estimated correlation")
print(corr_plot_ci)
ggsave("corr_plot_ci.jpg", corr_plot_ci, width = 6, height = 6)

## @knitr est_sigmas
fit_summary |>
  filter(grepl("^sigma", variable)) |>
  filter((q5 > 0 & q95 > 0) |
           (q5 < 0 & q95 < 0))

## @knitr logit_p
fit_summary |>
  filter(grepl("logit_p\\[", variable)) |>
  arrange(-rhat)

## @knitr logit_p
fit_summary |>
  filter(grepl("logit_psi\\[", variable)) |>
  arrange(-rhat)

## @knitr trace_plot
bayesplot::mcmc_trace(fit$draws(),
                      pars = c("lp__"))

## @knitr trace_plot_others
bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("^Omega_psi\\[[1-2]"))

bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("^Omega_theta\\[[1-2]"))

bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("^Omega_p\\[[1-2]"))



bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("sigma"))

bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("delta_p_star"))


bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("beta_psi_star"))

bayesplot::mcmc_trace(fit$draws(),
                      regex_pars = c("_p_star"))
