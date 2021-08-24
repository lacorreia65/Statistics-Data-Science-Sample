data {
  int<lower=1> N;
  int<lower=1> J; // number of counties
  int<lower=1,upper=J> county[N];
  vector[J] u;
  vector[N] x;
  vector[N] y;
}
parameters {
  vector[J] alpha;
  real beta;
  real mu_alpha;
  real gamma0;
  real gamma1;
  real<lower=0> sigma;
  real<lower=0> sigma_alpha;
}
model {
  vector[N] y_hat;
  vector[J] alpha_hat;
  for (i in 1:N)
    y_hat[i] = alpha[county[i]] + x[i] * beta;
  
  for(j in 1:J)
    alpha_hat[j] = gamma0 + gamma1*u[j];

  alpha ~ normal(alpha_hat, sigma_alpha);
  beta ~ normal(0, 1);
  gamma0 ~ normal(0,1);
  gamma1 ~ normal(0, 1);
  sigma ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);

  y ~ normal(y_hat, sigma);
}