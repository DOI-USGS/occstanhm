/* creates a cholesky_factor using a scale-able matrix */
matrix custom_cholesky(
		       int n_pred,
		       vector r_2,
		       row_vector l_omega) {       
  {
    matrix[n_pred, n_pred] L_Omega = rep_matrix(0, n_pred, n_pred);
    int start = 1;
    int end = 2;

    L_Omega[1,1] = 1.0;
    L_Omega[2,1] = 2.0 * r_2[1] - 1.0;
    L_Omega[2,2] = sqrt(1.0 - square(L_Omega[2,1]));
    for(k in 2:(n_pred - 1)) {
      int kp1 = k + 1;
      row_vector[k] l_row = segment(l_omega, start, k);
      real scale = sqrt(r_2[k] / dot_self(l_row));
      L_Omega[kp1, 1:k] = l_row[1:k] * scale;
      L_Omega[kp1, kp1] = sqrt(1.0 - r_2[k]);
      start = end + 1;
      end = start + k - 1;
    }

    return L_Omega;
  }
}
