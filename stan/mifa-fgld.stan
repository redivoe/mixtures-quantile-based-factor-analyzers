functions{
  real fgldvec_qf(real u, vector theta) {
    if(u < 1e-10)
      return -1e100;
    else if(u > 1 - 1e-10)
      return 1e100;
    else
      return theta[1] + theta[2]*((1 - theta[3])*log(u) - theta[3]*log1m(u) + theta[4]*u);
  }

  vector fgldvec_qf(vector u, vector theta) {
    int n = num_elements(u);
    vector[n] out;
    for (i in 1:n) {
      out[i] = fgldvec_qf(u[i], theta);
    }
    return out;
  }

  real h1(real delta, real kappa){
    return 1 + (pi() ^ 2 / 3 - 4) * delta * (1 - delta) + kappa / 2 * (1 + kappa / 6);
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
  array[K] matrix<lower=0, upper=1>[n, q] U;
  array[K] vector<lower=0, upper=1>[q] delta;
  array[K] vector<lower=0>[q] kappa;
}

model {
  array[K] matrix[n, p] eta;
  vector[K] log_lik;
  array[K] vector[q] alpha;
  array[K] vector[q] beta;
  array[K] matrix[n, q] Y;

  prop ~ dirichlet(alpha0_p);

  for(k in 1:K){
    mu[k] ~ std_normal();
    to_vector(Lambda[k]) ~ std_normal();
    to_vector(U[k]) ~ uniform(0, 1);
    psi[k] ~ inv_gamma(0.5, 0.5);
    delta[k] ~ uniform(0, 1);
    kappa[k] ~ lognormal(0, 1);
    for(j in 1:q){
      beta[k][j] = sqrt(1 / h1(delta[k][j], kappa[k][j]));
      alpha[k][j] = - beta[k][j] * (2 * delta[k][j] - 1 + kappa[k][j] / 2);
      Y[k][, j] = fgldvec_qf(U[k][, j], to_vector({alpha[k][j], beta[k][j], delta[k][j], kappa[k][j]}));
    }
    eta[k] = one_n * mu[k]' + Y[k] * Lambda[k]';
  }

  vector[K] log_prop = log(prop);

  for(i in 1:n){
    log_lik = log_prop;
    for(k in 1:K){
      for(l in 1:p){
        log_lik[k] += normal_lpdf(X[i, l] | eta[k][i, l], sqrt(psi[k][l]));
      }
    }
    target += log_sum_exp(log_lik);
  }
}
