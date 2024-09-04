functions {
#include ../functions/helper_functions.stan
#include ../functions/scale_cholesky_factor.stan
}
data {
  int<lower=1> n_units_species;
  int<lower=1> n_total_samples;
  int<lower=1> n_units;
  array[n_units_species] int<lower=1> n_revisits_per_unit_species;
  int<lower=1> n_total_revisits_species;
  array[n_total_revisits_species] int n_samples_per_revisit_species;
  array[n_total_revisits_species] int z_obs;
  array[n_total_samples] int a;
  array[n_units_species] int unit_species_revisit_start;
  array[n_units_species] int unit_species_revisit_stop;
  array[n_units_species] int unit_species_dat_start;
  array[n_units_species] int unit_species_dat_stop;
  array[n_total_revisits_species] int revisit_start_in_unit;
  array[n_total_revisits_species] int revisit_stop_in_unit;
  
  int<lower=1> m_beta;
  int<lower=1> m_alpha;

  // psi-level coffients
  array[n_units_species] int jj_psi;
  matrix[n_units_species, m_beta] x_psi;
  
  int<lower = 1> m_beta_star;
  matrix[n_units, m_beta_star] x_psi_star;

  // theta-level coffients
  array[n_units_species] int jj_theta;
  matrix[n_units_species, m_alpha] w_theta;
  
  int<lower = 1> m_alpha_star;
  matrix[n_units, m_alpha_star] w_theta_star;

  // priors for beta
  matrix[m_beta_star, m_beta] beta_psi_star_prior; // prior values for beta_star
  real<lower=0> beta_psi_star_sd_prior; // prior values for gamms's SD
  real<lower=0> eta_psi; // prior for Chol. matrix

  // priors for theta
  matrix[m_alpha_star, m_alpha] alpha_theta_star_prior; // prior values for beta_star
  real<lower=0> alpha_theta_star_sd_prior; // prior values for gamms's SD
  real<lower=0> eta_theta; // prior for Chol. matrix

  // priors for sigmas
  real<lower=0> sigma_psi_prior;
  real<lower=0> sigma_psi_prior_sd;

  real<lower=0> sigma_theta_prior;
  real<lower=0> sigma_theta_prior_sd;

  // reduce_sum setting
  int grainsize;
}
parameters {
  // psi-level parameters
  matrix[m_beta_star, m_beta] beta_psi_star; // group coeffs
  matrix[m_beta, n_units] z_psi;
  vector<lower=0, upper= pi()/2>[m_beta] tau_unif_psi;
  real<lower=0> sigma_psi;

  row_vector[choose(m_beta, 2) - 1]  l_omega_psi;  // do NOT init with 0 for all elements
  vector<lower = 0,upper = 1>[m_beta - 1] r_2_psi; // first element is not really a R^2 but is on (0,1)
  
  vector[n_units_species] logit_psi;

  // theta-level parameters
  matrix[m_alpha_star, m_alpha] alpha_theta_star; // group coeffs
  matrix[m_alpha, n_units] z_theta;
  vector<lower=0, upper= pi()/2>[m_alpha] tau_unif_theta;
  real<lower=0> sigma_theta;

  row_vector[choose(m_alpha, 2) - 1]  l_omega_theta;  // do NOT init with 0 for all elements
  vector<lower = 0,upper = 1>[m_alpha - 1] r_2_theta; // first element is not really a R^2 but is on (0,1)
  
 vector[n_units_species] logit_theta;
}
transformed parameters {
  // psi-level parameters
  matrix[n_units, m_beta] beta_psi;
  real<lower = 0> alpha_psi = eta_psi + (m_beta - 2) / 2.0;
  vector<lower = 0>[m_beta-1] shape_1_psi;
  vector<lower = 0>[m_beta-1] shape_2_psi;
  vector<lower=0>[m_beta] tau_psi; // prior scale
  
  matrix[m_beta, m_beta] L_Omega_psi =
    custom_cholesky(m_beta, r_2_psi, l_omega_psi);

  vector[n_units_species] logit_psi_in;

  // theta-level parameters
  matrix[n_units, m_alpha] alpha_theta;
  real<lower = 0> alpha_theta_psi = eta_theta + (m_alpha - 2) / 2.0;
  vector<lower = 0>[m_alpha-1] shape_1_theta;
  vector<lower = 0>[m_alpha-1] shape_2_theta;
  vector<lower=0>[m_alpha] tau_theta; // prior scale
  
  matrix[m_alpha, m_alpha] L_Omega_theta =
    custom_cholesky(m_alpha, r_2_theta, l_omega_theta);

  vector[n_units_species] logit_theta_in;
  //priors
  // psi-level
  shape_1_psi[1] = alpha_psi;
  shape_2_psi[1] = alpha_psi;

  for(k in 2:(m_beta-1)) {
    alpha_psi -= 0.5;
    shape_1_psi[k] = k / 2.0;
    shape_2_psi[k] = alpha_psi;
  }

  for (k in 1:m_beta) tau_psi[k] = 2.5 * tan(tau_unif_psi[k]);

  // theta-level
  shape_1_theta[1] = alpha_theta_psi;
  shape_2_theta[1] = alpha_theta_psi;

  for(k in 2:(m_alpha-1)) {
    alpha_theta_psi -= 0.5;
    shape_1_theta[k] = k / 2.0;
    shape_2_theta[k] = alpha_theta_psi;
  }

  for (k in 1:m_alpha) tau_theta[k] = 2.5 * tan(tau_unif_theta[k]);

  // modeled parameters
  // psi-level
  beta_psi = x_psi_star * beta_psi_star + (diag_pre_multiply(tau_psi, L_Omega_psi) * z_psi)';
  logit_psi_in = rows_dot_product(beta_psi[jj_psi], x_psi);

  // theta-level
  alpha_theta = w_theta_star * alpha_theta_star + (diag_pre_multiply(tau_theta, L_Omega_theta) * z_theta)';
  logit_theta_in = rows_dot_product(alpha_theta[jj_theta], w_theta);
}

model {
  vector[n_units_species] log_psi = log_inv_logit(logit_psi);
  vector[n_units_species] log1m_psi = log1m_inv_logit(logit_psi);
  
  l_omega_psi ~ std_normal();
  to_vector(z_psi) ~ std_normal();
  r_2_psi ~ beta(shape_1_psi, shape_2_psi);
  
  l_omega_theta ~ std_normal();
  to_vector(z_theta) ~ std_normal();
  r_2_theta ~ beta(shape_1_theta, shape_2_theta);

  to_vector(beta_psi_star) ~
    normal(to_vector(beta_psi_star_prior), beta_psi_star_sd_prior);
  to_vector(alpha_theta_star) ~
    normal(to_vector(alpha_theta_star_prior), alpha_theta_star_sd_prior);

  logit_psi ~ normal(logit_psi_in, sigma_psi);
  logit_theta ~ normal(logit_theta_in, sigma_theta);

  sigma_theta ~ normal(sigma_theta_prior, sigma_theta_prior_sd);
  sigma_psi ~ normal(sigma_psi_prior, sigma_psi_prior_sd);

  // loop over unit-species combination, that's a discrete sampling event
  
  target +=
    reduce_sum(unit_loop_theta, n_revisits_per_unit_species, grainsize,
               z_obs,
               a,
               log_psi,
               log1m_psi,
               logit_theta,
               unit_species_revisit_start,
               unit_species_revisit_stop,
               unit_species_dat_start,
               unit_species_dat_stop,
               revisit_start_in_unit,
               revisit_stop_in_unit,
               n_samples_per_revisit_species);
  
}


