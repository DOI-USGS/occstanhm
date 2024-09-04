functions {
#include ../functions/helper_functions.stan
}
data {
  int<lower=1> n_units_species;
  int<lower=1> n_total_samples;
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
  int grainsize;
}
parameters {
  vector[n_units_species] logit_psi;
  vector[n_units_species] logit_theta;
}
model {
  vector[n_units_species] log_psi = log_inv_logit(logit_psi);
  vector[n_units_species] log1m_psi = log1m_inv_logit(logit_psi);
  
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


