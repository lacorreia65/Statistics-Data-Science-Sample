data {
  int<lower=1> N;                       // number of observations
  vector[N] beds;                       // number of beds per roomtype
  vector[N] accom;                      // accomodates per roomtype
  vector[N] baths;                      // number of bathrooms
  vector[N] y;                          // log of price
}

parameters {
  real beta[4];              // beta coefficients
  real<lower=0> sigma;       // sigma of log_price
}

model {

  beta ~ normal(0, 1);
 
  y ~ normal(beta[1] + beta[2]*beds + beta[3]*accom + beta[4]*baths, sigma);

}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {  
    y_rep[n] = normal_rng(beta[1] + beta[2]*beds[n] + beta[3]*accom[n] + beta[4]*baths[n], sigma);            // Model 1-v3

    log_lik[n] = normal_lpdf( y[n] | beta[1] + beta[2]*beds[n] + beta[3]*accom[n] + beta[4]*baths[n], sigma);
  }

}
