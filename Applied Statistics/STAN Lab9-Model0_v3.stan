// Model 0 - $\theta_i$ is same in each region $= \theta$
data {
  int<lower=1> N;         // number of observations
  int<lower=0> y[N];      // observed no. of deaths
  vector[N] offset;       // expected no. of deaths
  vector[N] xc;           // Centered percent of males engaged in agriculture/forestry and fisheries
  int<lower=0> region[N]; // Region inside GDR
}

parameters {
  real theta;
}

model {

  theta ~ normal (0, 1);
 
  y ~ poisson(theta*offset);
}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {
    y_rep[n] = poisson_rng(theta*offset[n]);       // Model 2

    log_lik[n] = poisson_lpmf( y[n] | theta*offset[n]);
  }

}
