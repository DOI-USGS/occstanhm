## @knitr used_packages
library(tidyverse)
library(GGally)
library(R2jags)
library(occstanhm)
library(MCMCpack)
library(mclust)
library(corrplot)

## @knitr run_stan_model
n_chains <- 2
n_parallel_chains <- 2
n_threads_per_chain <- 4 
n_refresh <- 100

# for sampling numbers
n_warmup <- 400 # 2000
n_sample <- 400 # 2000

## @knitr simulate_data
n_units <- 50
n_unit_revisit_mean <- 10
n_unit_revisit_use <- 10
k_samples_mean <- 1
k_use <- 1
n_spp <- 3

## create covariance matrixes
a_sim_psi <- matrix(rbeta(n_spp^2, 0.01, 0.01) * 2 - 1, ncol = n_spp) 
sigma_use_psi <- t(a_sim_psi) %*% a_sim_psi
sigma_use_p <- matrix(0, nrow = n_spp, ncol = n_spp)
diag(sigma_use_p) <- 1

## create empty predictor matrices
x_psi_sim <-
  tibble(x_1 = rep(1, n_units), # intercept
         x_2 = rnorm(n_units, mean = 0.0, sd = 1.0),
         x_3 = rnorm(n_units, mean = 0.0, sd = 1.0))

v_p_sim <- NULL

psi_predictor_in <-
  rowSums(x_psi_sim)

sim_data <-
  sim_long_form_occ_2(n_units = n_units,
                      n_unit_revisit_mean = n_unit_revisit_mean,
                      n_spp = n_spp,
                      k_samples_mean = k_samples_mean,
                      spp_mu = rep(0.5, n_spp),
                      sigma_psi = sigma_use_psi,
                      sigma_p = sigma_use_p,
                      psi_predictor = psi_predictor_in,
                      p_predictor = NULL,
                      p_same = TRUE)

# known correlations among species
cor(sim_data$sigma_corr_psi$sigma_sim[,c(seq(1, n_spp, by = 1), n_spp)])

dat <-
  sim_data$spp_obs |>
  filter(sample <= k_use & revisit <= n_unit_revisit_use) |>
  mutate(species_id = as.integer(factor(species))) |>
  group_by(unit, unit_id, revisit, species, species_id) |>
  summarize(k_samples = n(),
            y = sum(y), .groups = "keep") |>
  mutate(z_obs = ifelse(sum(y) > 0, 1, 0)) |>
  ungroup()

clone_spp <-
  dat |>
  filter(species_id == max(species_id)) |>
  mutate(species_id = max(species_id) + 1L,
         species = paste0("Spp_", n_spp + 1))

dat <-
  dat |>
  bind_rows(clone_spp) |>
  arrange(unit_id, species_id) |>
  rowid_to_column("index")

revist_species_data <-
  dat |>
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
  ungroup() |>
  full_join(x_psi_sim |> rowid_to_column("unit"), by = "unit",
  )

revisit_species_at_unit_data  <-
  revist_species_data |>
  group_by(unit, species) |>
  summarize(n_revisits = n(),
            revisit_unit_start = min(index),
            revisit_unit_stop = max(index),
            .groups = "drop") |>
  ungroup()

unit_species_data <-
  dat |>
  group_by(unit, unit_id, species_id, species) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")

x_psi <- model.matrix(~ species - 1 + x_1 + x_2 + x_3, revist_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))

x_psi_star <-
  dat |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))

## @knitr check_psi
x_psi_hyper <-
  dat |>
  group_by(unit, species) |>
  summarize(n = n(), .groups = "drop") |>
  mutate(one = 1)

## check math for beta coefficient
beta_psi_check <- matrix(1,
                         nrow = dat |> distinct(unit) |> nrow(),
                         ncol = x_psi |> ncol())

logit_psi_check <- matrix(1,
                          nrow = revist_species_data |> nrow(),
                          ncol = 1)

## prior
beta_psi_star_prior <-
  matrix(0.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())

## @knitr check_beta_priors
## check math for prior
x_psi_star_check <-
  matrix(1,
         nrow = x_psi_star |> nrow(),
         ncol = x_psi_star |> ncol())

## @knitr check_beta_star
beta_psi_star_check <-
  matrix(1,
         nrow = 1,
         ncol = dat |> distinct(species) |> nrow())

## @knitr p_coef
# p-level
v_p <- model.matrix(~ 1, dat)
colnames(v_p) <- gsub("species", "", colnames(v_p))

## @knitr p_coef_check
## check math for delta coefficient
delta_p_check <- matrix(1,
                        nrow = dat |> distinct(unit) |> nrow(),
                        ncol = v_p |> ncol())

logit_p_check <- matrix(1,
                        nrow = dat |> nrow(),
                        ncol = 1)

## @knitr save_as_stan
stan_data <-
  # first are observed data
  list(y = dat |> pull(y),
       k_samples = dat |> pull(k_samples),
       z_obs = revist_species_data |> pull(z_obs),
       
       # next are the dimensions for data (including vectors)
       n_unit_species = unit_species_data |> nrow(),
       n_total_revisits = revist_species_data |> nrow(),
       n_total_samples = dat |> nrow(),
       
       n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
       n_samples = revist_species_data |> pull(n_samples),
       
       m_beta = x_psi |> ncol(),
       n_beta = dat |> distinct(unit) |> nrow(),
       
       m_beta_star = x_psi_star |> ncol(),
       n_beta_star = x_psi_star |> nrow(),
       
       m_delta = v_p |> ncol(),
       n_delta = 1, # dat |> distinct(unit_id) |> nrow(),
       
       # indexing
       unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
       unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),
       
       revisit_dat_start = revist_species_data |> pull(revisit_dat_start),
       revisit_dat_stop = revist_species_data |> pull(revisit_dat_stop),
       
       revisit_unit_start = revisit_species_at_unit_data |> pull(revisit_unit_start),
       revisit_unit_stop = revisit_species_at_unit_data |> pull(revisit_unit_stop),
       
       # psi-level inputs
       x_psi = x_psi,
       x_psi_star = x_psi_star,
       
       # p-level inputs
       v_p = v_p,
       
       # psi-level priors
       beta_psi_star_prior = beta_psi_star_prior,
       beta_psi_star_sd_prior = 1,
       eta_psi = 0.1,
       jj_psi = revist_species_data |> dplyr::pull(unit_id),
       
       # p-level priors
       delta_p_star_sd_prior = 0.1,
       eta_p = 1.0,
       jj_p = rep(1, nrow(dat)), # dat |> dplyr::pull(unit_id),
       
       # priors for sigams
       sigma_psi_prior = 0.5,
       sigma_psi_prior_sd = 0.05,
       sigma_p_prior = 1.4,
       sigma_p_prior_sd = 0.05,
       
       # setting for reduce_sum
       grainsize = 1,
       
       ## Store data (not used in Stan)
       dat = dat,
       unit_species_data = unit_species_data
  )

## Look at correlation
psi_cor <-
  stan_data$dat |>
  dplyr::group_by(unit, species_id) |>
  dplyr::summarize(z_obs = mean(z_obs),
                   .groups = "drop") |>
  dplyr::select(unit, species_id, z_obs) |>
  pivot_wider(names_from = unit,
              values_from = z_obs) |>
  dplyr::select(-species_id) |>
  as.matrix() |>
  t() |>
  cor()

# get inputs in ball park of data
init_fun <- function(){
  list(logit_p = rnorm(n = stan_data$n_total_samples,
                       mean = qlogis(
                         (stan_data$y/stan_data$k_samples + .01) * 0.9),
                       sd = 0.05),
       delta_p = matrix(rnorm(1), nrow = 1, ncol = 1),
       logit_psi = rnorm(n = stan_data$n_total_samples,
                         mean = qlogis((stan_data$z + 0.01) * 0.9),
                         sd = 0.05),
       beta_psi = rnorm(n = n_units * (n_spp + 3),
                        mean = 0,
                        sd = 0.5) |>
         matrix(nrow = n_units, ncol = n_spp + 3, byrow = TRUE),
       l_omega_psi = (rbeta(n = choose(stan_data$m_beta, 2) - 1, 3, 2) * 2 - 1),
       tau_unif_psi = runif(stan_data$m_beta, min = 0, max = pi / 2),
       beta_psi_star = 
         matrix(rnorm(x_psi_star |> ncol() *  x_psi |> ncol(), 0, sd = 0.1),
                nrow = x_psi_star |> ncol(),
                ncol = x_psi |> ncol()),
       sigma_psi = abs(rnorm(1, 0, 1)),
       sigma_p = abs(rnorm(1, 0, 1)),
       z_psi =
         matrix(rnorm(stan_data$m_beta * stan_data$n_beta_star, 0, sd = 1),
                nrow = stan_data$m_beta,
                ncol = stan_data$n_beta_star),
       r_2_psi =
         runif(n = stan_data$m_beta - 1L, 0, 1)
  )}

# occstanhm_fast_2(
fit <- occstanhm_fast_2(stan_data = stan_data,
                        n_warmup = n_warmup,
                        n_sample = n_sample,
                        n_chains = n_chains,
                        n_parallel_chains = n_parallel_chains,
                        init = init_fun,
                        n_threads_per_chain = n_threads_per_chain)

## @knitr fit_summary
fit_summary <-
  fit$summary(.cores = 15)

print(fit_summary)

print(fit$time())
print(fit$diagnostic_summary())

stan_posterior <-
  fit_summary |>
  filter(!grepl("star", variable)) |>
  filter(grepl("^Omega_psi\\[|beta_psi\\[|delta_p\\[", variable)) |>
  mutate(model = "stan",
         time = fit$time()$total)

## @knitr setup_tobler
### Latent variable multi-species co-occurrence model
modelFile <- system.file("./extdata/MSCoOcc_LVM.txt", package = "occstanhm")

#Number of latent variables to use
nlv <- 5

y_use <-
  stan_data$dat |>
  dplyr::select(unit, species_id, y) |>
  dplyr::group_by(unit, species_id) |>
  dplyr::summarize(y = sum(y), .groups = "drop") |>
  pivot_wider(names_from = unit,
              values_from = y) |>
  dplyr::select(-species_id) |>
  as.matrix() |>
  t()


#Specify the data
# N is number of species
# J is number of sites
# number of occasions for each site
# K<-rep(T,n.sites)
# y is simulated detection history
X_occ_in <- x_psi_sim |> as.matrix()
X_obs_in <- model.matrix(~ 1, x_psi_sim)
X_occ_in |> head()

n <- n_spp + 1
rownames(y_use) <- NULL
y <- y_use

Xocc <- X_occ_in
Vocc <- ncol(X_occ_in)

Xobs <- X_obs_in
Vobs <- ncol(X_obs_in)

k_in <- 
  stan_data$dat |>
  dplyr::group_by(unit, species) |>
  dplyr::summarize(k = sum(k_samples), .groups = "drop") |>
  dplyr::group_by(unit) |>
  dplyr::summarize(k = mean(k), .groups = "drop") |>
  dplyr::pull(k)

J <- nrow(y)

occ_data <-
  list(n = n_spp + 1,
       J = nrow(y),
       k = k_in,
       y = y_use,
       Xocc = X_occ_in,
       Vocc = ncol(X_occ_in),
       Xobs = X_obs_in,
       Vobs = ncol(X_obs_in),
       nlv = nlv)

#Specify the parameters to be monitored
occ_params <-
  c('z','u.b','v.b','mu.u.b','tau.u.b','mu.v.b','tau.v.b','LV','lv.coef')

#Specify the initial values
occ_inits = function() {
  lv.coef<-matrix(1,n,nlv)
  for(l in 1:nlv-1){
    lv.coef[l,(l+1):nlv]<-NA
  }
  list(
    u.b=matrix(rnorm(Vocc*n),c(n,Vocc)),
    # u.b=u.b,
    v.b=matrix(rnorm(Vobs*n),c(n,Vobs)),
    u=(y>0)-runif(1,0.1,0.8),
    LV=matrix(rnorm(nlv*J),J,nlv)
  )
}

## @knitr run_tobler
jags_time <-
  system.time({
    fit_raw <- jags.parallel(occ_data, occ_inits, occ_params, modelFile,
                             n.chains = 3L, n.iter= 15000L,
                             n.burnin = 10000L, n.thin = 5L)
  })

print(jags_time)

## @kntir compare_outputs
# draws_rvas_jags <- posterior::as_draws_rvars(fit_raw$BUGSoutput)
cols_keep <- grepl("LV|u.b|v.b", fit_raw$BUGSoutput$sims.matrix |> colnames())
jags_use <- posterior::rvar(fit_raw$BUGSoutput$sims.matrix[ , cols_keep])


jags_posterior <-
  posterior::summarize_draws(jags_use) |>
  arrange(-ess_bulk) |>
  mutate(model = "jags",
         time = jags_time[3])

#calculate the correlation matrix from the latent variables
lv.coef<-fit_raw$BUGSoutput$sims.list$lv.coef
cmall<-array(NA,dim=c(dim(lv.coef)[1], n_spp + 2,  n_spp + 2))
## RAE, had to use + 2 because of cloned species (originally there was + 1)

eps.res<-apply(lv.coef,c(1,2),function(x)1-sum(x^2))
for(i in 1:dim(lv.coef)[1]){ #for each mcmc sample
  cmall[i,,]<-cov2cor(tcrossprod(lv.coef[i,,]) + diag(apply(eps.res,2,mean)))
}
cmest <- apply(cmall,c(2,3),mean, na.rm = TRUE)
clower <- apply(cmall,c(2,3),quantile, 0.025, na.rm = TRUE)
cupper <- apply(cmall,c(2,3),quantile, 0.975, na.rm = TRUE)


#compare simulated and estimated correlation matrix
colnames(cmest) <- c(paste0("Spp_", seq(1, n_spp + 2)))
rownames(cmest) <- c(paste0("Spp_", seq(1, n_spp + 2)))

cmest_df <-
  cmest |>
  as.data.frame() |>
  rownames_to_column("species_row") |>
  pivot_longer(-species_row,
               names_to = "species_col",
               values_to = "correlation")

colnames(clower) <- c(paste0("Spp_", seq(1, n_spp + 2)))
rownames(clower) <- c(paste0("Spp_", seq(1, n_spp + 2)))

clower_df <-
  clower |>
  as.data.frame() |>
  rownames_to_column("species_row") |>
  pivot_longer(-species_row,
               names_to = "species_col",
               values_to = "l95")

colnames(psi_cor) <- c(paste0("Spp_", seq(1, n_spp + 1)))
rownames(psi_cor) <- c(paste0("Spp_", seq(1, n_spp + 1)))

psi_cor_df <-
  psi_cor |>
  as.data.frame() |>
  rownames_to_column("species_row") |>
  pivot_longer(-species_row,
               names_to = "species_col",
               values_to = "true")

colnames(cupper) <- c(paste0("Spp_", seq(1, n_spp + 2)))
rownames(cupper) <- c(paste0("Spp_", seq(1, n_spp + 2)))

cupper_df <-
  cupper |>
  as.data.frame() |>
  rownames_to_column("species_row") |>
  pivot_longer(-species_row,
               names_to = "species_col",
               values_to = "u95")

## extract from Stan
gsub_pattern <- "^Omega_(\\w+)\\[(\\d+),(\\d+)\\]"
spp_key <- 
  tibble(species_id = seq(1, ncol(stan_data$x_psi)),
         species_row = stan_data$x_psi |> colnames(),
         species_col = stan_data$x_psi |> colnames())

cor_stan <-
  fit_summary |>
  filter(grepl(gsub_pattern, variable)) |>
  mutate(species_row_id = as.integer(gsub(gsub_pattern, "\\2", variable)),
         species_col_id = as.integer(gsub(gsub_pattern, "\\3", variable)),
         level = gsub(gsub_pattern, "\\1", variable)) |>
  full_join(spp_key |> dplyr::select(species_id, species_row),
            by = c("species_row_id" = "species_id")) |>
  full_join(spp_key |> dplyr::select(species_id, species_col),
            by = c("species_col_id" = "species_id")) |>
  dplyr::select(species_row, species_col, correlation = mean,
                l95 = q5, u95 = q95) |>
  full_join(psi_cor_df, by = c("species_row", "species_col")) |>
  mutate(model = "stan") 

cor_jags <-
  cmest_df |>
  full_join(clower_df, by = c("species_row", "species_col")) |>
  full_join(cupper_df, by = c("species_row", "species_col")) |>
  full_join(psi_cor_df, by = c("species_row", "species_col")) |>
  mutate(model = "jags")


cor_all <-
  bind_rows(cor_jags, cor_stan) |>
  mutate(bias = true - correlation,
         coverage = ifelse(true >= l95 & true <= u95, TRUE, FALSE))

## @knitr compare_cor

cor_all |>
  filter(species_row != "Spp_5" & species_col != "Spp_5" ) |>
  mutate(cor_name = paste0(species_row, "-", species_col)) |>
  filter(grepl("Spp_", species_row) & grepl("Spp_", species_col)) |>
  ggplot(aes(x = cor_name, y = correlation,
             ymin = l95, ymax = u95, color = model)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_linerange(position = position_dodge(width = 1)) +
  coord_flip() +
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_manual(values = c("red", "blue")) +
  facet_grid(rows = vars(species_row, species_col), scales = "free") +
  geom_hline(aes(yintercept = true), size = 1.5, alpha = 0.5)



## @knitr compare_coef
#occupancy coefficients
BetaOcc <- fit_raw$BUGSoutput$sims.list$u.b

bocc1<-apply(BetaOcc,c(2,3),mean) #mean
bocc2<-apply(BetaOcc,c(2,3),quantile,0.025) #lower CI
bocc3<-apply(BetaOcc,c(2,3),quantile,0.975) #upper CI
bocc4<-((bocc2<0 & bocc3<0) | (bocc2>0 & bocc3>0)) #
bocc5<-bocc1*bocc4
colnames(bocc1)<-colnames(Xocc)
colnames(bocc5)<-colnames(Xocc)

##Occupancy probabilities for all species and sites
#occupancy probability on normal scale for each site and species
mu.psi<-Xocc %*% t(bocc1)
#occupancy probability for each site
psi<-1-pnorm(0,mu.psi,1)
# psi
colnames(psi) <- c(paste0("Spp_", seq(1, n_spp)), "Clone_Spp_9")

site_occ_jags <-
  psi |>
  as.data.frame() |>
  rownames_to_column("unit_id") |>
  pivot_longer(-unit_id,
               names_to = "species",
               values_to = "mean_jags") |>
  mutate(model = "jags")


spp_key_2 <- 
  tibble(species_id = seq(1, ncol(stan_data$x_psi)),
         species = stan_data$x_psi |> colnames())
gsub_pattern <- "beta_psi\\[(\\d+),(\\d+)\\]"

site_occ_prob <-
  stan_posterior |>
  filter(grepl("beta", variable)) |>
  mutate(unit_id = gsub(gsub_pattern, "\\1", variable),
         species_id = as.integer(gsub(gsub_pattern, "\\2", variable))) |>
  full_join(spp_key_2, by = "species_id") |>
  dplyr::select(unit_id, species, mean, q5, q95) |>
  mutate(mean_stan = plogis(mean),
         l95 = plogis(q5),
         u95 = plogis(q95)) |>
  dplyr::select(-q5, -q95) |>
  full_join(site_occ_jags, by = c("unit_id", "species"))

site_occ_prob
