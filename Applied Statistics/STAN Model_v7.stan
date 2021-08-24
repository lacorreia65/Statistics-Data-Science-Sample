data {
  int<lower=1> N;        // number of observations
  int<lower=1> J;        // number of counties
  int <lower=1,upper=J> ctynb[N];   // County number
  vector[N] y;    // log(activity)
  vector[N] x;    // floor
  vector[J] u;    // log(uppm)
}

parameters {
  real alpha[J];        // alphas
  real beta;            // beta
  real gamma0;          // 1st component of group-level mean
  real gamma1;          // 2nd component of group-level mean
  real<lower=0> sigma;  // error sd for Gaussian likelihood
  real<lower=0> sigma_alpha;  // error sd for alphas

}

model {
  vector[N] y_hat;
  vector[J] alpha_hat;
  for (i in 1:N) {
    y_hat[i] = alpha [ctynb[i]] + x[i] * beta;
  }
  for(j in 1:J) {
    alpha_hat[j] = gamma0 + gamma1*u[j];
  }
  alpha ~ normal (alpha_hat, sigma_alpha);
  beta ~ normal(0, 1) ;
  sigma ~ normal (0, 1) ;
  sigma_alpha ~ normal (0, 1) ;
  y ~ normal (y_hat, sigma) ;
}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[J] alpha_n;     


  for (n in 1:N) {
    for (j in 1:J) {
      alpha_n[j] = gamma0 + gamma1 * u[j] + beta * x[n];
    }   
    y_rep[n] = normal_rng(mean(alpha_n), sigma);
  }

}
