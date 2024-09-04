functions {
  #include ../functions/helper_functions.stan
}
data {
  int<lower=0> n_revisits;
  int<lower=0> n_units;
  array[n_revisits] int<lower=0, upper=1> z_obs;
  array[n_units] int<lower=1> n_revisits_per_unit;
  array[n_units] int<lower=1> unit_start;
  array[n_units] int<lower=1> unit_stop;
  int grainsize;
}
parameters {
  vector[n_units] logit_psi;
}
model {
  vector[n_units] log_psi = log_inv_logit(logit_psi);
  vector[n_units] log1m_psi = log1m_inv_logit(logit_psi);
  
  target +=
    reduce_sum(unit_loop, n_revisits_per_unit, grainsize,
               z_obs,
               log_psi,
               log1m_psi,
               unit_start,
	             unit_stop);
}

