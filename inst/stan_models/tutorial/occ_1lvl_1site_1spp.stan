data {
  int<lower=0> n_revisits;
  array[n_revisits] int z;
}
parameters {
  real logit_psi;
}
model {
  real log_psi = log_inv_logit(logit_psi);
  real log1m_psi = log1m_inv_logit(logit_psi);
  
  for (revisit_idx in 1:n_revisits) {
    if (z[revisit_idx] > 0){
      target += log_psi;
      } else {
      target += log1m_psi;
      }
  }
}

