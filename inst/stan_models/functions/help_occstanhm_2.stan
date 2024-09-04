real detection_bern_lpmf(
  array[] int y,
  int n_samples,
  vector logit_p
  ) {
    real target_temp = 0;
    
    for (sample_idx in 1:n_samples) {
        target_temp += 
          bernoulli_logit_lupmf(y[sample_idx] |
                               logit_p[sample_idx]);
      }
    return target_temp;
}

real revisit_z_loop(
  int n_revisits,
  array[] int z_obs,
  array[] int y,
  vector log_psi,
  vector log1m_psi,
  vector logit_p,
  array[] int n_samples,
  array[] int revisit_dat_start,
  array[] int revisit_dat_stop) {

    real target_temp = 0;
    // loop over revisits
    for (revisit_idx in 1:n_revisits) {
      
      int z_obs_idx = z_obs[revisit_idx];
      int n_samples_idx = n_samples[revisit_idx];
      real log_psi_idx = log_psi[revisit_idx];
      real log1m_psi_idx = log_psi[revisit_idx];
      int revisit_dat_start_idx = revisit_dat_start[revisit_idx];
      int revisit_dat_stop_idx = revisit_dat_stop[revisit_idx];

      array[n_samples_idx] int y_idx = y[revisit_dat_start_idx:revisit_dat_stop_idx];
      vector[n_samples_idx] logit_p_idx = logit_p[revisit_dat_start_idx:revisit_dat_stop_idx];
      
      if (z_obs_idx > 0){
        target_temp += log_psi_idx + 
          detection_bern_lpmf(y_idx | n_samples_idx,
                              logit_p_idx);

        } else {
        target_temp +=
          log_sum_exp(log1m_psi[revisit_idx],
                      log_psi_idx + 
                      detection_bern_lpmf(y_idx | n_samples_idx,
                                          logit_p_idx));
        }
      }
  return target_temp;
}

real unit_loop_psi(
  array[] int n_revisits,
  int start,
  int end,
  array[] int z_obs,
  array[] int y_obs,
  array[] int n_samples,
  vector log_psi,
  vector log1m_psi,
  vector logit_p,
  array[] int unit_spp_dat_start,
  array[] int unit_spp_dat_stop,
  array[] int revisit_dat_start,
  array[] int revisit_dat_stop,
  array[] int revisit_unit_start,
  array[] int revisit_unit_stop
) {
  real target_temp = 0.0;
  
  for (group_idx in 1:(1 + end - start)) {
    int idx_non_slice = group_idx - 1 + start;
    
    // slice out out unit-level variables
    int unit_spp_dat_start_idx = unit_spp_dat_start[idx_non_slice];
    int unit_spp_dat_stop_idx = unit_spp_dat_stop[idx_non_slice];
    
    int revisit_unit_start_idx = revisit_unit_start[idx_non_slice];
    int revisit_unit_stop_idx = revisit_unit_stop[idx_non_slice];
    int n_revisits_idx = n_revisits[group_idx];
    
    // slice out revisit-level variables
    array[n_revisits_idx] int revisit_dat_start_idx = revisit_dat_start[revisit_unit_start_idx:revisit_unit_stop_idx];
    array[n_revisits_idx] int revisit_dat_stop_idx = revisit_dat_stop[revisit_unit_start_idx:revisit_unit_stop_idx];

    array[n_revisits_idx] int z_obs_idx = z_obs[revisit_unit_start_idx:revisit_unit_stop_idx];
    array[n_revisits_idx] int n_samples_idx = n_samples[revisit_unit_start_idx:revisit_unit_stop_idx];
    
    vector[n_revisits_idx] log_psi_idx = log_psi[revisit_unit_start_idx:revisit_unit_stop_idx];
    vector[n_revisits_idx] log1m_psi_idx = log1m_psi[revisit_unit_start_idx:revisit_unit_stop_idx];

    // slice out to revisit-level
    array[sum(n_samples_idx)] int y_obs_idx = y_obs[unit_spp_dat_start_idx:unit_spp_dat_stop_idx];
    vector[sum(n_samples_idx)] logit_p_idx = logit_p[unit_spp_dat_start_idx:unit_spp_dat_stop_idx];

    target_temp +=
      revisit_z_loop(n_revisits_idx,
                     z_obs_idx,
                     y_obs_idx,
                     log_psi_idx,
                     log1m_psi_idx,
                     logit_p_idx,
                     n_samples_idx,
                     revisit_dat_start_idx,
                     revisit_dat_stop_idx);

  }
  return target_temp;
}


/* faster model code */

real revisit_z_loop_fast_lpmf(
  int n_revisits,
  array[] int z_obs,
  array[] int k,
  array[] int y,
  vector log_psi,
  vector log1m_psi,
  vector logit_p,
  array[] int revisit_dat_start,
  array[] int revisit_dat_stop) {

    real target_temp = 0;
    // loop over revisits
    for (revisit_idx in 1:n_revisits) {
      int z_obs_idx = z_obs[revisit_idx];
      real log_psi_idx = log_psi[revisit_idx];
      real log1m_psi_idx = log_psi[revisit_idx];
      int revisit_dat_start_idx = revisit_dat_start[revisit_idx];
      int revisit_dat_stop_idx = revisit_dat_stop[revisit_idx];

      int y_idx = y[revisit_idx];
      int k_idx = k[revisit_idx];
      real logit_p_idx = logit_p[revisit_idx];
      // am here, need to remove sample index

      if (z_obs_idx > 0){
        target_temp += log_psi_idx +
                       binomial_logit_lupmf(y_idx | k_idx, logit_p_idx);

        } else {
        target_temp +=
          log_sum_exp(log1m_psi[revisit_idx],
                      log_psi_idx + 
                      binomial_logit_lupmf(y_idx | k_idx, logit_p_idx));
        }
      }
  return target_temp;
}

real unit_loop_psi_fast(
  array[] int n_revisits,
  int start,
  int end,
  array[] int z_obs,
  array[] int k,
  array[] int y_obs,
  vector log_psi,
  vector log1m_psi,
  vector logit_p,
  array[] int unit_spp_dat_start,
  array[] int unit_spp_dat_stop,
  array[] int revisit_dat_start,
  array[] int revisit_dat_stop,
  array[] int revisit_unit_start,
  array[] int revisit_unit_stop
) {
  real target_temp = 0.0;
  
  for (group_idx in 1:(1 + end - start)) {
    int idx_non_slice = group_idx - 1 + start;
    
    // slice out out unit-level variables
    int unit_spp_dat_start_idx = unit_spp_dat_start[idx_non_slice];
    int unit_spp_dat_stop_idx = unit_spp_dat_stop[idx_non_slice];
    
    int revisit_unit_start_idx = revisit_unit_start[idx_non_slice];
    int revisit_unit_stop_idx = revisit_unit_stop[idx_non_slice];
    int n_revisits_idx = n_revisits[group_idx];
    
    // slice out revisit-level variables
    array[n_revisits_idx] int revisit_dat_start_idx = revisit_dat_start[revisit_unit_start_idx:revisit_unit_stop_idx];
    array[n_revisits_idx] int revisit_dat_stop_idx = revisit_dat_stop[revisit_unit_start_idx:revisit_unit_stop_idx];

    array[n_revisits_idx] int z_obs_idx = z_obs[revisit_unit_start_idx:revisit_unit_stop_idx];

    vector[n_revisits_idx] log_psi_idx = log_psi[revisit_unit_start_idx:revisit_unit_stop_idx];
    vector[n_revisits_idx] log1m_psi_idx = log1m_psi[revisit_unit_start_idx:revisit_unit_stop_idx];

    // slice out to revisit-level
    array[n_revisits_idx] int y_obs_idx = y_obs[unit_spp_dat_start_idx:unit_spp_dat_stop_idx];
    array[n_revisits_idx] int k_idx = k[unit_spp_dat_start_idx:unit_spp_dat_stop_idx];
    vector[n_revisits_idx] logit_p_idx = logit_p[unit_spp_dat_start_idx:unit_spp_dat_stop_idx];

    target_temp +=
      revisit_z_loop_fast_lpmf(n_revisits_idx |
                               z_obs_idx,
                               k_idx,
                               y_obs_idx,
                               log_psi_idx,
                               log1m_psi_idx,
                               logit_p_idx,
                               revisit_dat_start_idx,
                               revisit_dat_stop_idx);
  }
  return target_temp;
}



