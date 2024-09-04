functions {
#include ./functions/scale_cholesky_factor.stan
#include ./functions/help_occstanhm_2.stan
}
data {
  int<lower=1> n_unit_species; // number of unit-species combs
  int<lower=1> n_total_revisits; // number of all revisits across above
  int<lower=1> n_total_samples; // total number of y-level samples across above

  array[n_unit_species] int n_revisits; // revisits per unit
  array[n_total_revisits] int n_samples; // samples per revisit

  int<lower=1> m_beta; // number of rows
  int<lower=1> m_delta; // number of row

  int<lower=1> n_beta; // number of coefficients
  int<lower=1> n_delta; // number of coefficients

  array[n_total_samples] int<lower=0, upper=1> y; // detection during sampling event
  array[n_total_revisits] int<lower=0, upper=1> z_obs; // detection during revisit

  array[n_unit_species] int unit_spp_dat_start; // start of unit in obs-data
  array[n_unit_species] int unit_spp_dat_stop; // stop of unit in obs-data
  
  array[n_total_revisits] int revisit_dat_start; // start of revist in obs-data
  array[n_total_revisits] int revisit_dat_stop; // stop of revist in obs-data
  
  array[n_unit_species] int revisit_unit_start; // start of revist each unit
  array[n_unit_species] int revisit_unit_stop; // stop of revist each unit

  // psi-level coffients
  array[n_total_revisits] int jj_psi;
  matrix[n_total_revisits, m_beta] x_psi;

  int<lower=1> m_beta_star;
  int<lower=1> n_beta_star;

  matrix[n_beta_star, m_beta_star] x_psi_star;

  // p-level coffients
  array[n_total_samples] int jj_p;
  matrix[n_total_samples, m_delta] v_p;

  int<lower = 1> m_delta_star;
  int<lower = 1> n_delta_star;

  matrix[n_delta_star, m_delta_star] v_p_star;

  // priors for psi
  matrix[m_beta_star, m_beta] beta_psi_star_prior; // prior values for beta_star
  real<lower=0> beta_psi_star_sd_prior; // prior values for gamms's SD
  real<lower=0> eta_psi; // prior for Chol. matrix

  // priors for p
  matrix[m_delta_star, m_delta] delta_p_star_prior; // prior values for delta_star
  real<lower=0> delta_p_star_sd_prior; // prior values for gamms's SD
  real<lower=0> eta_p; // prior for Chol. matrix

  // priors for sigmas
  real<lower=0> sigma_psi_prior;
  real<lower=0> sigma_psi_prior_sd;

  real<lower=0> sigma_p_prior;
  real<lower=0> sigma_p_prior_sd;

  // reduce_sum setting
  int grainsize;
}
parameters {
  // // psi-level parameters
  matrix[m_beta_star, m_beta] beta_psi_star; // group coeffs
  real<lower=0> sigma_psi;

  matrix[m_beta, n_beta_star] z_psi;
  vector<lower=0, upper= pi()/2>[m_beta] tau_unif_psi;

  row_vector[choose(m_beta, 2) - 1]  l_omega_psi;  // do NOT init with 0 for all elements
  vector<lower = 0,upper = 1>[m_beta - 1] r_2_psi; // first element is not really a R^2 but is on (0,1)

  // p-level parameters
  matrix[m_delta_star, m_delta] delta_p_star; // group coeffs
  real<lower=0> sigma_p;

  matrix[m_delta, n_delta_star] z_p;
  vector<lower=0, upper= pi()/2>[m_delta] tau_unif_p;

  row_vector[choose(m_delta, 2) - 1]  l_omega_p;  // do NOT init with 0 for all elements
  vector<lower = 0,upper = 1>[m_delta - 1] r_2_p; // first element is not really a R^2 but is on (0,1)

  vector[n_total_revisits] logit_psi;
  vector[n_total_samples] logit_p;
}
transformed parameters {
  // psi-level parameters
  matrix[n_beta, m_beta] beta_psi;
  real<lower = 0> beta_psi_psi = eta_psi + (m_beta - 2) / 2.0;
  vector<lower = 0>[m_beta-1] shape_1_psi;
  vector<lower = 0>[m_beta-1] shape_2_psi;
  vector<lower=0>[m_beta] tau_psi; // prior scale

  matrix[m_beta, m_beta] L_Omega_psi =
    custom_cholesky(m_beta, r_2_psi, l_omega_psi);

  vector[n_total_revisits] logit_psi_in;

  // p-level parameters
  matrix[n_delta, m_delta] delta_p;
  real<lower = 0> beta_p_psi = eta_p + (m_delta - 2) / 2.0;
  vector<lower = 0>[m_delta-1] shape_1_p;
  vector<lower = 0>[m_delta-1] shape_2_p;
  vector<lower=0>[m_delta] tau_p; // prior scale

  matrix[m_delta, m_delta] L_Omega_p =
    custom_cholesky(m_delta, r_2_p, l_omega_p);
  vector[n_total_samples] logit_p_in;

  // psi-level
  shape_1_psi[1] = beta_psi_psi;
  shape_2_psi[1] = beta_psi_psi;

  for(k in 2:(m_beta-1)) {
    beta_psi_psi -= 0.5;
    shape_1_psi[k] = k / 2.0;
    shape_2_psi[k] = beta_psi_psi;
  }

  for (k in 1:m_beta) tau_psi[k] = 2.5 * tan(tau_unif_psi[k]);

  // p-level
  shape_1_p[1] = beta_p_psi;
  shape_2_p[1] = beta_p_psi;

  for(k in 2:(m_delta-1)) {
    beta_p_psi -= 0.5;
    shape_1_p[k] = k / 2.0;
    shape_2_p[k] = beta_p_psi;
  }

  for (k in 1:m_delta) tau_p[k] = 2.5 * tan(tau_unif_p[k]);

  // psi-level
  beta_psi = x_psi_star * beta_psi_star + (diag_pre_multiply(tau_psi, L_Omega_psi) * z_psi)';
  logit_psi_in = rows_dot_product(beta_psi[jj_psi], x_psi);

  // p-level
  delta_p = v_p_star * delta_p_star + (diag_pre_multiply(tau_p, L_Omega_p) * z_p)';
  logit_p_in = rows_dot_product(delta_p[jj_p], v_p);
}

model {
  vector[n_total_revisits] log_psi = log_inv_logit(logit_psi);
  vector[n_total_revisits] log1m_psi = log1m_inv_logit(logit_psi);

  l_omega_psi ~ std_normal();
  to_vector(z_psi) ~ std_normal();
  r_2_psi ~ beta(shape_1_psi, shape_2_psi);
  
  l_omega_p ~ std_normal();
  to_vector(z_p) ~ std_normal();
  r_2_p ~ beta(shape_1_p, shape_2_p);

  to_vector(beta_psi_star) ~
    normal(to_vector(beta_psi_star_prior), beta_psi_star_sd_prior);
  to_vector(delta_p_star) ~
    normal(to_vector(delta_p_star_prior), delta_p_star_sd_prior);

  sigma_psi ~ normal(sigma_psi_prior, sigma_psi_prior_sd);
  sigma_p ~ normal(sigma_p_prior, sigma_p_prior_sd);

  logit_psi ~ normal(logit_psi_in, sigma_psi);
  logit_p ~ normal(logit_p_in, sigma_p);

  // loop over unit-species combination, that's a discrete sampling event
  target +=
    reduce_sum(unit_loop_psi,
               n_revisits, grainsize,
               z_obs,
               y,
               n_samples,
               log_psi,
               log1m_psi,
               logit_p,
               unit_spp_dat_start,
               unit_spp_dat_stop,
               revisit_dat_start,
               revisit_dat_stop,
               revisit_unit_start,
               revisit_unit_stop);
}
generated quantities {
  matrix[m_beta, m_beta] Omega_psi;
  matrix[m_delta, m_delta] Omega_p;

  Omega_psi = multiply_lower_tri_self_transpose(L_Omega_psi);
  Omega_p     = multiply_lower_tri_self_transpose(L_Omega_p);
}
