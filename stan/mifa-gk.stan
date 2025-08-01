functions {
  vector gk_qf_vec(vector z, real g, real k) {
      return (1 + 0.8 * tanh(0.5 * g * z)) .* z .* (1 + z .* z) .^ k;
  }
}

data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0, upper=p> q;
  int<lower=1> K;
  matrix[n, p] X;
}

transformed data{
  vector[n] one_n = ones_vector(n);
  vector[K] alpha0_p = rep_vector(10, K);
}

parameters {
  simplex[K] prop;
  array[K] vector[p] mu;
  array[K] matrix[p, q] Lambda;
  array[K] vector<lower=0>[p] psi;
  array[K] matrix[n, q] Z;
  array[K] vector[q] g;
  array[K] vector<lower=0>[q] kappa;
}

model {
  array[K] matrix[n, q] Y;
  array[K] matrix[n, p] eta;
  vector[K] log_prop = log(prop);

  prop ~ dirichlet(alpha0_p);

  for (k in 1:K) {
    mu[k] ~ std_normal();
    to_vector(Lambda[k]) ~ std_normal();
    to_vector(Z[k]) ~ std_normal();
    psi[k] ~ inv_gamma(0.5, 0.5);
    g[k] ~ std_normal();
    kappa[k] ~ lognormal(0, 1);
    for(j in 1:q){
      Y[k][, j] = gk_qf_vec(Z[k][, j], g[k][j], kappa[k][j]);
    }
    eta[k] = one_n * mu[k]' + Y[k] * Lambda[k]';
  }

  for (i in 1:n) {
    vector[K] log_lik = log_prop;
    for (k in 1:K) {
      for(l in 1:p){
        log_lik[k] += normal_lpdf(X[i, l] | eta[k][i, l], sqrt(psi[k][l]));
      }
    }
    target += log_sum_exp(log_lik);
  }
}
