functions {
#include ../functions/helper_functions.stan
}
data {
  int<lower=1> n_total_samples;
  int<lower=1> n_revisits;
  array[n_revisits] int n_samples_per_revisit;
  array[n_revisits] int z_obs;
  array[n_total_samples] int a;
  array[n_revisits] int unit_start;
  array[n_revisits] int unit_stop;
}
parameters {
  real logit_psi;
  real logit_theta;
}
model {
  real log_psi = log_inv_logit(logit_psi);
  real log1m_psi = log1m_inv_logit(logit_psi);

  for (revisit_idx in 1:n_revisits) {
    if (z_obs[revisit_idx] > 0){
      
      target +=
        detection_lpmf(
          a[unit_start[revisit_idx]:unit_stop[revisit_idx]] |
          log_psi, logit_theta,
          n_samples_per_revisit[revisit_idx]);
      } else {
      target += 
        log_sum_exp(log1m_psi,
                    detection_lpmf(
                      a[unit_start[revisit_idx]:unit_stop[revisit_idx]] |
                      log_psi, logit_theta,
                      n_samples_per_revisit[revisit_idx]));
        }
  }
}


