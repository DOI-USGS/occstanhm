functions {
  #include ../functions/helper_functions.stan
}
data {
  int<lower=1> n_revisits;
  int<lower=1> n_units_species;
  int<lower=1> m_beta;
  int<lower=1> n_units;
  array[n_revisits] int<lower=0, upper=1> z_obs;
  array[n_units_species] int<lower=1> n_revisits_per_unit;
  array[n_units_species] int<lower=1> unit_start;
  array[n_units_species] int<lower=1> unit_stop;
  array[n_units_species] int jj_psi;
  matrix[n_units_species, m_beta] x_psi;
  int grainsize;
}
parameters {
  matrix[n_units, m_beta] beta_psi;
}
transformed parameters {
  vector[n_units_species] logit_psi;
  logit_psi = rows_dot_product(beta_psi[jj_psi], x_psi);
}
model {
  vector[n_units_species] log_psi = log_inv_logit(logit_psi);
  vector[n_units_species] log1m_psi = log1m_inv_logit(logit_psi);
  
  target +=
    reduce_sum(unit_loop, n_revisits_per_unit, grainsize,
               z_obs,
               log_psi,
               log1m_psi,
               unit_start,
               unit_stop);
}

