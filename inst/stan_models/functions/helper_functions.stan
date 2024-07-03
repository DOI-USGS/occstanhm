real detection_lpmf(
  array[] int a,
  real log_psi,
  real logit_theta,
  int n_samples_per_revisit
  ) {
    real target_temp = 0;

    for (sample_idx in 1:n_samples_per_revisit) {
          target_temp +=
            binomial_logit_lupmf(a[sample_idx] | 1, logit_theta);
        }
    
    return log_psi + target_temp;
}

real detection_y_lpmf(
  array[] int a_obs,
  array[] int y,
  array[] int k_subsamples,
  real log_psi,
  real log_theta,
  real log1m_theta,
  real logit_p,
  int n_samples_per_revisit
  ) {
    real target_temp = 0;

    for (sample_idx in 1:n_samples_per_revisit) {
      if (a_obs[sample_idx] > 0) {
        target_temp += log_theta +
          binomial_logit_lupmf(y[sample_idx] | k_subsamples[sample_idx],
                               logit_p);
      } else {
        target_temp +=
          log_sum_exp(log1m_theta,
            log_theta +
            binomial_logit_lupmf(y[sample_idx] | k_subsamples[sample_idx],
                                logit_p));
      }
    } 
    return log_psi + target_temp;
}

real revisit_loop(
  int n_revisits_per_site,
  array[] int z_obs,
  real log_psi,
  real log1m_psi) {
    real target_temp = 0;
    // loop over revisits
    for (revisit_idx in 1:n_revisits_per_site) {
      if (z_obs[revisit_idx] > 0){
        target_temp += log_psi;
        } else {
        target_temp += log1m_psi;
        }
      }
  return target_temp;
}

real unit_loop(
  array[] int n_revisits_per_site, //needs to be length number of units
  int start,
  int end,
  array[] int z_obs,
  vector log_psi,
  vector log1m_psi,
  array[] int unit_start,
  array[] int unit_stop) {
  real target_temp = 0;

    for (unit_idx in 1:(1 + end - start)) {
      int idx_non_slice = unit_idx - 1 + start;
      // loop over revisits
       target_temp += revisit_loop(
         n_revisits_per_site[unit_idx],
         z_obs[unit_start[idx_non_slice]:unit_stop[idx_non_slice]],
         log_psi[idx_non_slice],
         log1m_psi[idx_non_slice]);
      }
  return target_temp;
}

real revisit_theta_loop(
  int n_revisits_per_unit,
  array[] int z_obs,
  array[] int a,
  array[] int revisit_start_in_unit,
  array[] int revisit_stop_in_unit,
  real log_psi,
  real logit_theta,
  real log1m_psi,
  array[] int n_samples_per_revisit
  ) {
    real target_temp = 0;
    // loop over revisits
    for (revisit_idx in 1:(n_revisits_per_unit)) {
      if (z_obs[revisit_idx] > 0) {
        target_temp +=
          detection_lpmf(
            a[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]] |
            log_psi, logit_theta,
            n_samples_per_revisit[revisit_idx]);
      } else {
        target_temp +=
          log_sum_exp(log1m_psi,
            detection_lpmf(
              a[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]] |
              log_psi, logit_theta,
              n_samples_per_revisit[revisit_idx]));
      }
    }
  return target_temp;
}

real unit_loop_theta(
  array[] int n_revisits_per_unit,
  int start,
  int end,
  array[] int z_obs,
  array[] int a,
  vector log_psi,
  vector log1m_psi,
  vector logit_theta,
  array[] int unit_revisit_start,
  array[] int unit_revisit_stop,
  array[] int unit_dat_start,
  array[] int unit_dat_stop,
  array[] int revisit_start_in_unit,
  array[] int revisit_stop_in_unit,
  array[] int n_samples_per_revisit
) {
  real target_temp = 0.0;

  for (unit_idx in 1:(1 + end - start)) {
    int idx_non_slice = unit_idx - 1 + start;
    target_temp += 
      revisit_theta_loop(n_revisits_per_unit[unit_idx],
                         z_obs[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                         a[unit_dat_start[idx_non_slice]:unit_dat_stop[idx_non_slice]],
                         revisit_start_in_unit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                         revisit_stop_in_unit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                         log_psi[idx_non_slice],
                         logit_theta[idx_non_slice],
                         log1m_psi[idx_non_slice],
                         n_samples_per_revisit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]]
                         );

  }
  return target_temp;
}

real revisit_p_loop(
  int n_revisits_per_unit,
  array[] int z_obs,
  array[] int a_obs,
  array[] int y,
  array[] int k_subsamples,
  array[] int revisit_start_in_unit,
  array[] int revisit_stop_in_unit,
  real log_psi,
  real log_theta,
  real log1m_psi,
  real log1m_theta,
  real logit_p,
  array[] int n_samples_per_revisit
  ) {
    real target_temp = 0;
    // loop over revisits
    for (revisit_idx in 1:(n_revisits_per_unit)) {
      if (z_obs[revisit_idx] > 0) {
        target_temp +=
          detection_y_lpmf(
            a_obs[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]] |
            y[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]],
            k_subsamples[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]],
            log_psi, log_theta, log1m_theta, logit_p,
            n_samples_per_revisit[revisit_idx]);
      } else {
        target_temp +=
          log_sum_exp(log1m_psi,
            detection_y_lpmf(
              a_obs[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]] |
              y[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]],
              k_subsamples[revisit_start_in_unit[revisit_idx]:revisit_stop_in_unit[revisit_idx]],
              log_psi, log_theta, log1m_theta, logit_p,
              n_samples_per_revisit[revisit_idx]));
      }
    }
  return target_temp;
}

real unit_loop_p(
  array[] int n_revisits_per_unit,
  int start,
  int end,
  array[] int z_obs,
  array[] int a_obs,
  array[] int y,
  array[] int k_subsamples,
  vector log_psi,
  vector log1m_psi,
  vector log_theta,
  vector log1m_theta,
  vector logit_p,
  array[] int unit_revisit_start,
  array[] int unit_revisit_stop,
  array[] int unit_dat_start,
  array[] int unit_dat_stop,
  array[] int revisit_start_in_unit,
  array[] int revisit_stop_in_unit,
  array[] int n_samples_per_revisit
) {
  real target_temp = 0.0;

  for (unit_idx in 1:(1 + end - start)) {
    int idx_non_slice = unit_idx - 1 + start;
    target_temp += 
      revisit_p_loop(n_revisits_per_unit[unit_idx],
                     z_obs[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                     a_obs[unit_dat_start[idx_non_slice]:unit_dat_stop[idx_non_slice]],
                     y[unit_dat_start[idx_non_slice]:unit_dat_stop[idx_non_slice]],
                     k_subsamples[unit_dat_start[idx_non_slice]:unit_dat_stop[idx_non_slice]],
                     revisit_start_in_unit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                     revisit_stop_in_unit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]],
                     log_psi[idx_non_slice],
                     log_theta[idx_non_slice],
                     log1m_psi[idx_non_slice],
                     log1m_theta[idx_non_slice],
                     logit_p[idx_non_slice],
                     n_samples_per_revisit[unit_revisit_start[idx_non_slice]:unit_revisit_stop[idx_non_slice]]
                     );

  }
  return target_temp;
}

real detection_y_theta_2_lpmf(
  array[] int a_obs,
  array[] int y,
  array[] int k_subsamples,
  vector log_theta,
  vector log1m_theta,
  vector logit_p,
  int n_samples_per_grouping
  ) {
    real target_temp = 0;

    for (sample_idx in 1:n_samples_per_grouping) {
      if (a_obs[sample_idx] > 0) {
        target_temp +=
          log_theta[sample_idx] +
            binomial_logit_lupmf(y[sample_idx] | k_subsamples[sample_idx],
                                 logit_p[sample_idx]);
      } else {
        target_temp +=
          log_sum_exp(log1m_theta[sample_idx],
            log_theta[sample_idx] +
            binomial_logit_lupmf(y[sample_idx] | k_subsamples[sample_idx],
                                logit_p[sample_idx]));
      }
    }
    return target_temp;
}

real unit_loop_p_2(
  array[] int n_samples_per_grouping,
  int start,
  int end,
  array[] int a_obs,
  array[] int y,
  array[] int k_subsamples,
  vector log_theta,
  vector log1m_theta,
  vector logit_p,
  array[] int grouping_start,
  array[] int grouping_stop
) {
  real target_temp = 0.0;

  for (group_idx in 1:(1 + end - start)) {
    int idx_non_slice = group_idx - 1 + start;
    
    target_temp +=
      detection_y_theta_2_lpmf(
                       a_obs[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]] |
                       y[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]],
                       k_subsamples[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]],
                       log_theta[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]],
                       log1m_theta[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]],
                       logit_p[grouping_start[idx_non_slice]:grouping_stop[idx_non_slice]],
                       n_samples_per_grouping[group_idx]
                       );

  }
  return target_temp;
}
