#' Three-level, multispecies hierarchical occupancy model
#'
#' occstanhm_3 is a three-level, hierarchical occupancy model.
#'
#' @param stan_data Data formatted for Stan.
#' @return fit Model fit from cmdstanr.
#'
#' @export
#' @examples
#'
#' \dontrun{
#' library(occstanhm)
#' library(tidyverse)
#' library(ggthemes)
#'
#' set.seed(123)
#' n_units <- 5
#' n_unit_revisit_mean <- 3
#' k_samples_mean <- 3
#' k_subsamples <- 6
#' n_spp <- 4
#'
#' sim_data <-
#'   sim_long_form_occ_3(n_units = n_units,
#'                       n_unit_revisit_mean = n_unit_revisit_mean,
#'                       n_spp = n_spp,
#'                       k_samples_mean = k_samples_mean,
#'                       k_subsamples = k_subsamples,
#'                       spp_mu = rep(0.75, n_spp),
#'                       sigma_psi = matrix(c(1.0,  0.0,  0.0,  0.0,
#'                                            0.0,  1.0,  1.0, -0.7,
#'                                            0.0,  1.0,  1.0, -0.7,
#'                                            0.0, -0.7, -0.7,  1.0),
#'                                          nrow = n_spp, ncol = n_spp),
#'                       sigma_p = matrix(c(1.0,  0.0,  0.0,  0.0,
#'                                          0.0,  1.0,  1.0, -0.7,
#'                                          0.0,  1.0,  1.0, -0.7,
#'                                          0.0, -0.7, -0.7,  1.0),
#'                                        nrow = n_spp, ncol = n_spp))
#'
#' sim_data$spp_obs <-
#'   sim_data$spp_obs |>
#'   arrange(unit, species, revisit) |>
#'   mutate(unit = paste("Unit", unit),
#'          unit_id = as.integer(unit_id),
#'          species_id = as.integer(factor(species))) |>
#'   rowid_to_column("index")
#'
#'
#' sim_data$spp_obs |>
#'   filter(unit == "Unit 1") |>
#'   group_by(unit, species, revisit) |>
#'   summarize(total_y = sum(y), .groups = "drop") |>
#'   ggplot(aes(x = revisit, y = species, fill = total_y)) +
#'   geom_tile() +
#'   scale_fill_gradient("Samples\nwith\ndetections",
#'                       low = "white", high = "black") +
#'   theme_bw() +
#'   theme(strip.background = element_blank())
#'
#' many_plot <-
#'   sim_data$spp_obs |>
#'   group_by(unit_id, species, revisit) |>
#'   summarize(total_y = sum(y), .groups = "drop") |>
#'   ggplot(aes(x = revisit, y = species, fill = total_y)) +
#'   geom_tile() +
#'   scale_fill_gradient("Samples\nwith\ndetections",
#'                       low = "white", high = "black") +
#'   theme_bw() +
#'   theme(strip.background = element_blank(),
#'         panel.spacing = unit(2, "lines"))+
#'   facet_wrap(vars(unit_id),
#'              labeller = label_bquote("Unit"~.(unit_id))) +
#'   ylab("Species") +
#'   xlab("Revisit")
#' print(many_plot)
#' revist_species_data <-
#'   sim_data$spp_obs |>
#'   ungroup() |>
#'   group_by(unit, species) |>
#'   mutate(min_idx = min(index) - 1) |>
#'   ungroup() |>
#'   group_by(unit, unit_id, species, species_id, revisit) |>
#'   summarize(z_obs = ifelse(sum(y) > 0, 1, 0),
#'             n_samples = n(),
#'             revisit_dat_start = min(index - min_idx),
#'             revisit_dat_stop = max(index - min_idx),
#'             .groups = "drop") |>
#'   rowid_to_column("index") |>
#'   ungroup()
#' revist_species_data |>
#'   print(n = 20)
#'
#' revisit_species_at_unit_data  <-
#'   revist_species_data |>
#'   group_by(unit, species) |>
#'   summarize(n_revisits = n(),
#'             revisit_unit_start = min(index),
#'             revisit_unit_stop = max(index),
#'             .groups = "drop") |>
#'   ungroup()
#'
#' revisit_species_at_unit_data
#' unit_species_data <-
#'   sim_data$spp_obs |>
#'   group_by(unit, unit_id, species, species_id) |>
#'   summarize(psi_obs = mean(z_obs),
#'             n_samples = n(),
#'             unit_spp_dat_start = min(index),
#'             unit_spp_dat_stop = max(index),
#'             .groups = "drop")
#' x_psi <- model.matrix(~ species - 1, revist_species_data)
#' colnames(x_psi) <- gsub("species", "", colnames(x_psi))
#' x_psi |> dim()
#' x_psi |> head()
#' x_psi_star <-
#'   sim_data$spp_obs |>
#'   distinct(unit) |>
#'   model.matrix(object = ~ 1)
#' colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))
#' beta_psi_star_prior <-
#'   matrix(1.0,
#'          nrow = x_psi_star |> ncol(),
#'          ncol = x_psi |> ncol())
#' # theta-level
#' w_theta <- model.matrix(~ species - 1, sim_data$spp_obs)
#' colnames(w_theta) <- gsub("species", "", colnames(w_theta))
#' ## hyper parameters
#' w_theta_star <-
#'   sim_data$spp_obs |>
#'   distinct(unit) |>
#'   model.matrix(object = ~ 1)
#'
#' alpha_theta_star_prior <-
#'   matrix(1,
#'          nrow = w_theta_star |> ncol(),
#'          ncol = w_theta |> ncol())
#' # p-level
#' v_p <- model.matrix(~ species - 1, sim_data$spp_obs)
#' colnames(v_p) <- gsub("species", "", colnames(v_p))
#' ## check math for delta coefficient
#' delta_p_check <- matrix(1,
#'                         nrow = sim_data$spp_obs |> distinct(unit) |> nrow(),
#'                         ncol = v_p |> ncol())
#' ## hyper parameters
#' v_p_star <-
#'   sim_data$spp_obs |>
#'   distinct(unit) |>
#'   model.matrix(object = ~ 1)
#'
#' delta_p_star_prior <-
#'   matrix(1,
#'          nrow = v_p_star |> ncol(),
#'          ncol = v_p |> ncol())
#' stan_data <-
#'   # first are observed data
#'   list(y = sim_data$spp_obs |> pull(y),
#'        a_obs = sim_data$spp_obs |> pull(a_obs),
#'        k_subsamples = sim_data$spp_obs |> pull(k_subsample),
#'        z_obs = revist_species_data |> pull(z_obs),
#'
#'        # next are the dimensions for data (including vectors)
#'        n_unit_species = unit_species_data |> nrow(),
#'        n_total_revisits = revist_species_data |> nrow(),
#'        n_total_samples = sim_data$spp_obs |> nrow(),
#'
#'        n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
#'        n_samples = revist_species_data |> pull(n_samples),
#'
#'        m_beta = x_psi |> ncol(),
#'        n_beta = sim_data$spp_obs |> distinct(unit) |> nrow(),
#'
#'        m_beta_star = x_psi_star |> ncol(),
#'        n_beta_star = x_psi_star |> nrow(),
#'
#'        m_alpha = w_theta |> ncol(),
#'        n_alpha = sim_data$spp_obs |> distinct(unit_id) |> nrow(),
#'
#'        m_alpha_star = w_theta_star |> ncol(),
#'        n_alpha_star = w_theta_star |> nrow(),
#'
#'        m_delta = v_p |> ncol(),
#'        n_delta = sim_data$spp_obs |> distinct(unit_id) |> nrow(),
#'
#'        m_delta_star = v_p_star |> ncol(),
#'        n_delta_star = v_p_star |> nrow(),
#'
#'        # indexing
#'        unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
#'        unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),
#'
#'        revisit_dat_start = revist_species_data |> pull(revisit_dat_start),
#'        revisit_dat_stop = revist_species_data |> pull(revisit_dat_stop),
#'
#'        revisit_unit_start =
#'           revisit_species_at_unit_data |> pull(revisit_unit_start),
#'        revisit_unit_stop =
#'           revisit_species_at_unit_data |> pull(revisit_unit_stop),
#'
#'        # psi-level inputs
#'        x_psi = x_psi,
#'        x_psi_star = x_psi_star,
#'
#'        # psi-level inputs
#'        w_theta = w_theta,
#'        w_theta_star = w_theta_star,
#'
#'        # p-level inputs
#'        v_p = v_p,
#'        v_p_star = v_p_star,
#'
#'        # psi-level priors
#'        beta_psi_star_prior = beta_psi_star_prior,
#'        beta_psi_star_sd_prior = 0.1,
#'        eta_psi = 1.0,
#'        jj_psi = revist_species_data |> dplyr::pull(unit_id),
#'
#'        # theta-level priors
#'        alpha_theta_star_prior = alpha_theta_star_prior,
#'        alpha_theta_star_sd_prior = 0.1,
#'        eta_theta = 1.0,
#'        jj_theta = sim_data$spp_obs |> dplyr::pull(unit_id),
#'
#'        # p-level priors
#'        delta_p_star_prior = delta_p_star_prior,
#'        delta_p_star_sd_prior = 0.1,
#'        eta_p = 1.0,
#'        jj_p = sim_data$spp_obs |> dplyr::pull(unit_id),
#'
#'        # priors for sigams
#'        sigma_psi_prior = 2.5,
#'        sigma_psi_prior_sd = 0.05,
#'
#'        sigma_theta_prior = 1.75,
#'        sigma_theta_prior_sd = 0.05,
#'
#'        sigma_p_prior = 1.75,
#'        sigma_p_prior_sd = 0.05,
#'
#'        # setting for reduce_sum
#'        grainsize = 1,
#'
#'        ## Store data (not used in Stan)
#'        dat = sim_data$spp_obs,
#'        unit_species_data = unit_species_data
#'   )
#'
#' n_chains <- 3
#' n_parallel_chains <- 3
#' n_threads_per_chain <- 5
#' n_refresh <- 100
#'
#' # for sampling numbers
#' n_warmup <- 200 # 10000
#' n_sample <- 200 # 10000
#'
#' fit <- occstanhm_3(stan_data,
#'                    n_chains = n_chains,
#'                    n_parallel_chains = n_parallel_chains,
#'                    n_threads_per_chain = n_threads_per_chain,
#'                    n_refresh = n_refresh,
#'                    n_warmup = n_warmup,
#'                    n_sample = n_sample)
#' fit_summary <-
#'   fit$summary(.cores = 15)
#' fit_summary
#' }

occstanhm_3 <- function(
    stan_data,
    n_chains = 2,
    n_parallel_chains = 2,
    n_threads_per_chain = 4,
    n_refresh = 100,
    n_warmup = 200,
    n_sample = 200,
    ...) {

  stan_file <- system.file("./stan_models/occstanhm_3.stan",
                           package = "occstanhm")

  mod <- cmdstanr::cmdstan_model(stan_file,
                                 cpp_options = list(stan_threads = TRUE))
  fit <- mod$sample(stan_data,
                    chains = n_chains,
                    parallel_chains = n_parallel_chains,
                    threads_per_chain = n_threads_per_chain,
                    refresh = n_refresh,
                    iter_warmup = n_warmup,
                    iter_sampling = n_sample,
                    ...)
  return(fit)
}
