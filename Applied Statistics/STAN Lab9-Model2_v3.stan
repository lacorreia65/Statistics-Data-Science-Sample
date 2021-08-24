// Model 2 - $\theta_i$ is different in each region and modeled hierarchically
data {
  int<lower=1> N;         // number of observations
  int<lower=0> y[N];      // observed no. of deaths
  vector[N] offset;       // expected no. of deaths
  vector[N] xc;           // Centered percent of males engaged in agriculture/forestry and fisheries
  int<lower=0> region[N]; // Region inside GDR
}

parameters {
  vector [N] alpha;     
  real mu;              
  real beta;            
  real<lower=0> sigma_mu;
}

model {
  vector[N] log_lambda;

  for (i in 1:N) {
    log_lambda[i] = alpha[i] + beta*xc[i] + offset[i];
  }
  
  alpha ~ normal (mu, sigma_mu);

  mu ~normal (0, 1);
  beta ~ normal(0, 1);
  sigma_mu ~ normal (0, 1);

  y ~ poisson_log(log_lambda);
}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {
    y_rep[n] = poisson_log_rng(alpha[n] + beta*xc[n] + offset[n]);       // Model 2

    log_lik[n] = poisson_log_lpmf( y[n] | alpha[n] + beta*xc[n] + offset[n]);
  }

}
