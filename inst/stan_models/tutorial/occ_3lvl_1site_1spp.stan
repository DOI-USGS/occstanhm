functions {
#include ../functions/helper_functions.stan
}
data {
  int<lower=1> n_total_samples;
  int<lower=1> n_revisits;
  array[n_revisits] int n_samples_per_revisit;
  array[n_revisits] int z_obs;
  array[n_total_samples] int a_obs;
  array[n_total_samples] int y;
  array[n_total_samples] int k_subsamples;
  array[n_revisits] int unit_start;
  array[n_revisits] int unit_stop;
}
parameters {
  real logit_psi;
  real logit_theta;
  real logit_p;
}
model {
  real log_psi = log_inv_logit(logit_psi);
  real log1m_psi = log1m_inv_logit(logit_psi);
  
  real log_theta = log_inv_logit(logit_theta);
  real log1m_theta = log1m_inv_logit(logit_theta);
  
  for (revisit_idx in 1:n_revisits) {
    if (z_obs[revisit_idx] > 0){
      
      target +=
        detection_y_lpmf(
          a_obs[unit_start[revisit_idx]:unit_stop[revisit_idx]] |
          y[unit_start[revisit_idx]:unit_stop[revisit_idx]],
          k_subsamples[unit_start[revisit_idx]:unit_stop[revisit_idx]],
          log_psi, log_theta, log1m_theta, logit_p,
          n_samples_per_revisit[revisit_idx]);
      } else {
      target += 
        log_sum_exp(log1m_psi,
                    detection_y_lpmf(
                      a_obs[unit_start[revisit_idx]:unit_stop[revisit_idx]] |
                      y[unit_start[revisit_idx]:unit_stop[revisit_idx]],
                      k_subsamples[unit_start[revisit_idx]:unit_stop[revisit_idx]],
                      log_psi, log_theta, log1m_theta, logit_p,
                      n_samples_per_revisit[revisit_idx]));
        }
  }
}


