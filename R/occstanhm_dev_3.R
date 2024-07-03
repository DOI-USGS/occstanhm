#' Three-level occupancy for multiple species with correlation
#'
#' This model is the "complete" three-level model. The model has three levels
#' with correlated error structures at all three levels. The model also allows
#' for different coefficients to be estimated for each level unlike the
#' building block models that estimated the same value for each unit.
#'
#' This model also includes priors on the estimated coefficients in addition to
#' a prior.
#' @param stan_data Data formatted for Stan.
#' @return fit Model fit from cmdstanr.
#'
#' @export
#' @examples
#' # Increase the number of chains to 4
#' # Increase the number of parallel chains and threads based upon your CPUs
#' # For example, an 8 CPU machine might do 4 chains on 2 threads
#' # or 2 chains on 2 threads that are run back-to-back.
#' # Likewise, the number of warmup and samples to 100s for debugging and
#' # 1,000s or 10,000 for actual mode fits.
#' \dontrun{
#' ## Load packages used
#' library(tidyverse)
#' library(ggthemes)
#'
#' ## Set simulation details
#' # Number of units
#' n_units_sim <- 6 # max is 26 due to number of letters
#'
#' # Number of species
#' n_species_sim <- 3
#'
#' # Number of subunits (such as habitats within unit)
#' n_subunit <- 3
#'
#' # number of markers
#' n_markers <- 2
#'
#' # number of molecular samples
#' n_subsamples <- 6
#'
#' # Poisson mean for number of stochastic samples per site
#' sample_lamba <- 10
#'
#' ## Simulate data
#' set.seed(12345)
#' # create species-unit data
#' dat_unit <-
#'   expand_grid(unit = letters[seq(1, n_units_sim)],
#'               species = paste("species", seq(1, n_species_sim))) |>
#'   mutate(psi = rbeta(n(), 6, 3))
#' dat_unit
#'
#' # create subunit details
#' dat_subunit <-
#'   expand_grid(unit = letters[seq(1, n_units_sim)],
#'               subunit = paste("subunit", seq(1, n_subunit))) |>
#'   full_join(dat_unit, by = "unit", relationship = "many-to-many") |>
#'   arrange(unit, subunit, species) |>
#'   mutate(z = rbinom(n(), 1, psi),
#'          theta = rbeta(n(), 6, 3),
#'          unit_subunit = paste(unit, subunit, sep = "-"),
#'          subunit_species = paste(subunit, species, sep = "-"))
#'
#' dat_subunit
#'
#' dat_subunit |>
#'   group_by(unit, species) |>
#'   summarize(psi = mean(psi),
#'             psi_obs = mean(z))
#'
#' # create sample data
#' n_samples_per_subunit <-
#'   rpois(n = n_units_sim * n_subunit, lambda = sample_lamba) + 1L
#'
#' dat_sample <-
#'   tibble(
#'     unit_subunit = rep(dat_subunit |>
#'                         distinct(unit_subunit) |>
#'                         pull(unit_subunit),
#'                         times = n_samples_per_subunit)) |>
#'   full_join(dat_subunit, by = "unit_subunit",
#'             relationship = "many-to-many") |>
#'   mutate(a = z * rbinom(n(), 1, theta))
#'
#' dat_sample |>
#'   arrange(unit, species, subunit) |>
#'   print(n = 10)
#'
#' # summarize data to see simualted theta values
#' theta_obs <-
#'   dat_sample |>
#'   group_by(unit, subunit, subunit_species, species, theta) |>
#'   summarize(n = n(), z = mean(z), theta_obs = mean(a),
#'             .groups = "drop") |>
#'   filter(z > 0) |>
#'   select(unit, subunit_species, subunit, species, n, z, theta, theta_obs)
#'
#' theta_obs |>
#'   print(n = 10)
#'
#' # Create/simulate full data set
#' # Arranged by unit, species, and then subunit
#' # Subunit arrange is tricky part (needs to be last)
#' dat <-
#'   expand_grid(unit = letters[seq(1, n_units_sim)],
#'               species = paste("species", seq(1, n_species_sim)),
#'               plain_marker = paste("M", seq(1, n_markers))) |>
#'   mutate(marker = paste(species, plain_marker, sep = "-"),
#'          unit_id = as.integer(factor(unit))) |>
#'   select(-plain_marker) |>
#'   mutate(p = rbeta(n(), 3, 3)) |>
#'   full_join(dat_sample, by = c("unit", "species"),
#'             relationship = "many-to-many") |>
#'   select(unit_id, unit, subunit, species, marker, z, a, p) |>
#'   mutate(subunit_species = paste(subunit, species, sep = "-"),
#'          k_subsamples = n_subsamples,
#'          y = a * rbinom(n(), k_subsamples, prob = p),
#'          a_obs = ifelse(y > 0, 1, 0)) |>
#'   group_by(unit_id,unit, subunit, species) |>
#'   mutate(z_obs = ifelse(sum(a_obs) > 0, 1, 0),
#'          subunit_id = as.integer(factor(subunit))) |>
#'   ungroup() |>
#'   arrange(unit, species, subunit) |>
#'   rowid_to_column("index")
#'
#' dat |>
#'   print(n = 10)
#'
#' # look at simulated or "observed" psi values
#' psi_obs <-
#'   dat |>
#'   group_by(unit_id, unit, species) |>
#'   filter(a > 0) |>
#'   summarize(p = mean(p), p_obs = mean(y/k_subsamples),
#'             .groups = "drop")
#'
#' psi_obs |>
#'   print(n = Inf)
#'
#' dat |>
#'   filter(z_obs != z)
#'
#' dat |>
#'   filter(a_obs != a)
#'
#' # Summarize simulated data for model inputs
#' # Need to make sure group by unit first
#' # This gives us  the start and stop for each species/unit combination
#' unit_species_summary <-
#'   dat |>
#'   group_by(unit_id, unit, species) |>
#'   summarize(psi_obs = mean(z_obs),
#'             unit_species_dat_start = min(index),
#'             unit_species_dat_stop = max(index),
#'             n_revisits_per_site = n(),
#'             .groups = 'keep') |>
#'   mutate(logit_psi = qlogis(psi_obs))
#'
#' unit_species_summary |>
#'   print(n = Inf)
#'
#' # examine indexing to compare to Stan model
#' # (used for debugging and understanding)
#' test_idx <- 1
#' test_range <-
#'   seq(unit_species_summary$unit_species_dat_start[test_idx],
#'       unit_species_summary$unit_species_dat_stop[test_idx])
#'
#' dat$z[test_range] |> mean()
#'
#' # This gives the start and stop for each subunit (or revisit)
#' # This includes the z_obs for each subunit
#'
#' # summarize each subunit or revisit to a unit
#' revisit_species_summary <-
#'   dat |>
#'   group_by(unit_id, unit, species, subunit) |>
#'   summarize(n_samples_per_grouping = n(),
#'             revisit_start = min(index),
#'             revisit_stop = max(index),
#'             z_obs = mean(z_obs),
#'             .groups = "drop") |>
#'   group_by(unit_id, unit, species) |>
#'   mutate(
#'     revisit_start_in_unit_species = revisit_start - min(revisit_start) + 1L,
#'     revisit_stop_in_unit_species = revisit_stop - min(revisit_start) + 1L)
#'
#' revisit_species_summary
#'
#' # summarize for slicing index
#' unit_species_revisit_summary <-
#'   revisit_species_summary |>
#'   rowid_to_column("index") |>
#'   group_by(unit_id, unit, species) |>
#'   summarize(unit_species_revisit_start = min(index),
#'             unit_species_revisit_stop = max(index),
#'             n_revisits_per_grouping = n(),
#'             psi_obs = mean(z_obs),
#'             .groups = "drop")
#'
#'
#' unit_species_revisit_summary |>
#'   print(n = Inf)
#'
#' # look at indexes to understand and debug Stan code
#' revisit_species_summary |>
#'   filter(unit == "a" & species == "species 1") |>
#'   summarize(mean(z_obs))
#'
#' test_idx <- 1
#' test_range <-
#'   seq(unit_species_revisit_summary$unit_species_revisit_start[test_idx],
#'       unit_species_revisit_summary$unit_species_revisit_stop[test_idx])
#'
#' revisit_species_summary$z_obs[test_range] |>
#'   mean()
#'
#' ## Create parameters for model
#' # psi-level
#' x_psi <- model.matrix(~ species - 1, revisit_species_summary)
#' x_psi |> head()
#' x_psi |> dim()
#'
#' unit_species_summary_hyper <-
#'   unit_species_summary |>
#'   group_by(unit) |>
#'   summarize(n = n(),
#'             hyper_psi = mean(psi_obs),
#'             hyper_logit = mean(logit_psi))
#'
#' x_psi_star <-
#'   model.matrix(~ 1, unit_species_summary_hyper)
#'
#' beta_psi_star_prior <-
#'   matrix(0,
#'          nrow = x_psi_star |> ncol(),
#'          ncol = x_psi |> ncol())
#'
#' # theta-level
#' w_theta <- model.matrix(~ subunit_species - 1, dat)
#' colnames(w_theta) <- gsub("subunit_species", "", colnames(w_theta))
#' w_theta |> head()
#' w_theta |> dim()
#'
#' w_theta_star <-
#'   model.matrix(~ 1, unit_species_summary_hyper)
#' w_theta_star |> dim()
#'
#' alpha_theta_star_prior <-
#'   matrix(0,
#'          nrow = w_theta_star |> ncol(),
#'          ncol = w_theta |> ncol())
#'
#' # p-level
#' v_p <- model.matrix(~ marker - 1, dat)
#' v_p |> head()
#' v_p |> dim()
#'
#' v_p_star <-
#'   model.matrix(~ 1, unit_species_summary_hyper)
#'
#' delta_p_star_prior <-
#'   matrix(0,
#'          nrow = v_p_star |> ncol(),
#'          ncol = v_p |> ncol())
#'
#' ## Save data in list for Stan
#' stan_data <-
#'   # observation and latent observation data
#'   list(
#'     y = dat |> pull(y),
#'     k_subsamples = dat |> pull(k_subsamples),
#'     a_obs = dat |> pull(a_obs),
#'     z_obs = revisit_species_summary |> pull(z_obs),
#'
#'     # indexing
#'     unit_species_dat_start = unit_species_summary |>
#'       pull(unit_species_dat_start),
#'     unit_species_dat_stop = unit_species_summary |>
#'       pull(unit_species_dat_stop),
#'     revisit_start_in_unit = revisit_species_summary |>
#'       pull(revisit_start_in_unit_species),
#'     revisit_stop_in_unit = revisit_species_summary |>
#'       pull(revisit_stop_in_unit_species),
#'     unit_species_revisit_start = unit_species_revisit_summary |>
#'       pull(unit_species_revisit_start),
#'     unit_species_revisit_stop = unit_species_revisit_summary |>
#'       pull(unit_species_revisit_stop),
#'
#'     # all different dimensions, starting with "n_"
#'     n_units = dat |> distinct(unit) |> nrow(),
#'     n_groupings = unit_species_summary |> nrow(),
#'     n_total_samples = dat |> nrow(),
#'     n_revisits_per_grouping = unit_species_revisit_summary |>
#'       pull(n_revisits_per_grouping),
#'     n_samples_per_grouping = revisit_species_summary |>
#'       pull(n_samples_per_grouping),
#'     n_total_groupings = revisit_species_summary |> nrow(),
#'
#'     # psi-level inputs
#'     x_psi = x_psi,
#'     m_beta = x_psi |> ncol(),
#'     x_psi_star = x_psi_star,
#'     m_beta_star = x_psi_star |> ncol(),
#'     beta_psi_star_prior = beta_psi_star_prior,
#'     beta_psi_star_sd_prior = 1,
#'     eta_psi = 1,
#'     jj_psi = revisit_species_summary |> dplyr::pull(unit_id),
#'
#'     # theta-level inputs
#'     w_theta = w_theta,
#'     m_alpha = w_theta |> ncol(),
#'     w_theta_star = w_theta_star,
#'     m_alpha_star = w_theta_star |> ncol(),
#'     alpha_theta_star_prior = alpha_theta_star_prior,
#'     alpha_theta_star_sd_prior = 1,
#'     eta_theta = 1,
#'     jj_theta = dat |> dplyr::pull(unit_id),
#'
#'     # theta-level inputs
#'     v_p = v_p,
#'     m_delta = v_p |> ncol(),
#'     v_p_star = v_p_star,
#'     m_delta_star = v_p_star |> ncol(),
#'     delta_p_star_prior = delta_p_star_prior,
#'     delta_p_star_sd_prior = 1,
#'     eta_p = 1,
#'     jj_p = dat |> dplyr::pull(unit_id),
#'
#'        # priors
#'        sigma_psi_prior = 1.5,
#'        sigma_psi_prior_sd = 0.05,
#'        sigma_theta_prior = 1.5,
#'        sigma_theta_prior_sd = 0.05,
#'        sigma_p_prior = 2.0,
#'        sigma_p_prior_sd = 0.05,
#'
#'        # reduce_sum setting
#'        grainsize = 1
#'   )
#'
#' ## Set simulation details
#' # this includes within and among chain parallel settings
#' # (something that most Stan users may not be familiar with)
#' n_chains <- 4
#' n_parallel_chains <- 4
#' n_threads_per_chain <- 4
#' n_refresh <- 100
#'
#' # for sampling numbers
#' # will want to set to be at least 2000 each, more would be better.
#' ## Takes ~1 hr to run with 2000 for each
#' # used 600 and 300 for demonstration
#' n_warmup <- 600 # 2000
#' n_sample <- 300 # 2000
#'
#' ## Compile and run model
#' fit <- occstanhm_dev_3(stan_data,
#'                        n_chains = n_chains,
#'                        n_parallel_chains = n_parallel_chains,
#'                        n_threads_per_chain = n_threads_per_chain,
#'                        n_refresh = n_refresh,
#'                        n_warmup = n_warmup,
#'                        n_sample = n_sample)
#'
#' ## Summarize model fit
#' fit_summary <-
#'   fit$summary(.cores = 15)
#'
#' fit_summary
#'
#' # Look at psi-level estimatse from the model
#' unit_species_revisit_summary_2 <-
#'   unit_species_revisit_summary|>
#'   group_by(species, unit) |>
#'   summarize(psi_obs = mean(psi_obs),
#'             .groups = "drop")
#'
#' unit_species_revisit_summary_2 |>
#'   select(unit, species, psi_obs)
#'
#' beta_detail <-
#'   fit_summary |>
#'   filter(grepl("beta_psi\\[", variable)) |>
#'   select(variable, mean, median, q5, q95) |>
#'   mutate(mean = mean,
#'          median = median,
#'          l95 = plogis(q5),
#'          prob = plogis(median),
#'          u95 = plogis(q95)) |>
#'   select(-mean, -median, -q5, -q95) |>
#'   cbind(unit_species_revisit_summary |>
#'           arrange(species, unit) |>
#'           select(unit, species, psi_obs)) |>
#'   as_tibble() |>
#'   mutate(diff = prob - psi_obs) |>
#'   select(variable, unit, species, l95, prob, u95, psi_obs, diff)
#'
#' beta_detail
#'
#' beta_detail |>
#'   ggplot(aes(x = `unit`,
#'              color = species,
#'              y = prob, ymin = l95, ymax = u95)) +
#'   geom_point(position = position_dodge(width = 0.5)) +
#'   geom_linerange(position = position_dodge(width = 0.5)) +
#'   scale_color_colorblind() +
#'   coord_flip() +
#'   theme_bw() +
#'   geom_point(aes(y = psi_obs),
#'              position = position_dodge(width = 0.5), shape = "square",
#'              size = 3, alpha = 0.5) +
#'   theme(strip.background = element_blank()) +
#'   ylim(c(0, 1))
#'
#' ## Look at theta-level model outputs
#' dat |>
#'   mutate(theta_cal =
#'            fit_summary |>
#'            filter(grepl("logit_theta\\[", variable)) |>
#'            pull(mean)) |>
#'   group_by(unit, species) |>
#'   summarize(t_1 = mean(theta_cal)) |>
#'   mutate(plogis(t_1))
#'
#' stan_data$a_obs
#'
#' subunit_species_summary_hyper <-
#'   dat |>
#'   group_by(unit, subunit, species, subunit_species) |>
#'   summarize(n = n(), theta_rough = mean(a_obs),
#'             .groups = "drop") |>
#'   select(unit, subunit, species, subunit_species, theta_rough) |>
#'   mutate(unit_id = as.integer(factor(unit)),
#'          subunit_species_id = as.integer(factor(subunit_species)))
#'
#' subunit_species_summary_hyper
#'
#' alpha_detail <-
#'   fit_summary |>
#'   filter(grepl("alpha_theta\\[", variable)) |>
#'   select(variable, q5, median, q95) |>
#'   mutate(l95 = plogis(q5),
#'          u95 = plogis(q95),
#'          prob = plogis(median)) |>
#'   select(-q5, -q95, -median) |>
#'   mutate(
#'     unit_id = as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]" , "\\1",
#'                               variable)),
#'     subunit_species_id = as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]",
#'                                          "\\2", variable))) |>
#'   full_join(subunit_species_summary_hyper,
#'             by = c("unit_id", "subunit_species_id")) |>
#'   select(variable, unit, subunit, species, l95, prob, u95, theta_rough) |>
#'   mutate(diff = theta_rough - prob)
#'
#' alpha_detail
#'
#' ggplot(alpha_detail, aes(x = subunit,
#'                          color = species,
#'                          y = prob, ymin = l95, ymax = u95)) +
#'   geom_point(position = position_dodge(width = 0.5)) +
#'   geom_linerange(position = position_dodge(width = 0.5)) +
#'   scale_color_colorblind() +
#'   facet_wrap(vars(`unit`), nrow = 2) +
#'   coord_flip() +
#'   theme_bw() +
#'   geom_point(aes(y = theta_rough),
#'              position = position_dodge(width = 0.5), shape = "square",
#'              size = 4) +
#'   theme(strip.background = element_blank()) +
#'   ylim(c(0, 1))
#'
#' ## Look at p-level model estimate
#' dat_rough_p <-
#'   dat |>
#'   filter(a_obs > 0) |>
#'   mutate(marker_id = as.integer(factor(marker))) |>
#'   group_by(unit, unit_id, species, marker, marker_id) |>
#'   summarize(p_rough = mean(y/k_subsamples), .groups = "drop")
#'
#' dat_p_all <-
#'   dat |>
#'   mutate(marker_id = as.integer(factor(marker))) |>
#'   distinct(unit, unit_id, species, marker, marker_id) |>
#'   full_join(dat_rough_p, by = c("unit", "unit_id",
#'                                 "species",
#'                                 "marker",
#'                                 "marker_id"))
#'
#' p_details <-
#'   fit_summary |>
#'   filter(grepl("delta_p\\[", variable)) |>
#'   mutate(unit_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]" , "\\1",
#'                                    variable)),
#'          marker_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]",
#'                                      "\\2", variable))) |>
#'   full_join(dat_p_all, by = c("marker_id", "unit_id")) |>
#'   mutate(prob = plogis(median),
#'          l95 = plogis(q5),
#'          u95 = plogis(q95)) |>
#'   select(variable, unit, species, marker, l95, prob, u95, p_rough) |>
#'   mutate(diff = p_rough - prob)
#'
#' p_details |>
#'   print(n = Inf)
#'
#' p_details |>
#'   ggplot(aes(x = marker,
#'              color = species,
#'              y = prob, ymin = l95, ymax = u95)) +
#'   geom_point(position = position_dodge(width = 0.5)) +
#'   geom_linerange(position = position_dodge(width = 0.5)) +
#'   scale_color_colorblind() +
#'   facet_grid(rows = vars(`unit`)) +
#'   coord_flip() +
#'   theme_bw() +
#'   geom_point(aes(y = p_rough),
#'              position = position_dodge(width = 0.5), shape = "square",
#'              size = 4) +
#'   theme(strip.background = element_blank())  +
#'   ylim(c(0, 1))
#'
#' ## Look at estimated correlations from the model
#' fit_summary |>
#'   filter(grepl("^Omega_psi\\[", variable)) |>
#'   select(variable, q5, median, q95)
#'
#' fit_summary |>
#'   filter(grepl("^Omega_theta\\[", variable)) |>
#'   select(variable, q5, median, q95)
#'
#' fit_summary |>
#'   filter(grepl("^Omega_p\\[", variable)) |>
#'   select(variable, q5, median, q95)
#'
#' }
occstanhm_dev_3 <- function(stan_data,
                            n_chains = 2,
                            n_parallel_chains = 2,
                            n_threads_per_chain = 4,
                            n_refresh = 100,
                            n_warmup = 200,
                            n_sample = 200,
                            ...) {

  stan_file <- system.file("./stan_models/tutorial/occstanhm_dev_3.stan",
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
