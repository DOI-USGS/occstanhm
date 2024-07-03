functions {
#include ../functions/helper_functions.stan
}
data {
  int<lower=1> n_units;
  int<lower=1> n_total_samples;
  array[n_units] int<lower=1> n_revisits_per_unit;
  int<lower=1> n_total_revisits;
  array[n_total_revisits] int n_samples_per_revisit;
  array[n_total_revisits] int z_obs;
  array[n_total_samples] int a;
  array[n_units] int unit_revisit_start;
  array[n_units] int unit_revisit_stop;
  array[n_units] int unit_dat_start;
  array[n_units] int unit_dat_stop;
  array[n_total_revisits] int revisit_start_in_unit;
  array[n_total_revisits] int revisit_stop_in_unit;
  int grainsize;
}
parameters {
  vector[n_units] logit_psi;
  vector[n_units] logit_theta;
}
model {
  vector[n_units] log_psi = log_inv_logit(logit_psi);
  vector[n_units] log1m_psi = log1m_inv_logit(logit_psi);

  target +=
    reduce_sum(unit_loop_theta, n_revisits_per_unit, grainsize,
               z_obs,
               a,
               log_psi,
               log1m_psi,
               logit_theta,
               unit_revisit_start,
               unit_revisit_stop,
               unit_dat_start,
               unit_dat_stop,
               revisit_start_in_unit,
               revisit_stop_in_unit,
               n_samples_per_revisit);
  
}
