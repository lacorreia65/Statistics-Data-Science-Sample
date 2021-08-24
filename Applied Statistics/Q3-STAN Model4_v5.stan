data {
  int<lower=1> N;                       // number of observations
  int<lower=1> L;                       // Number of Locations
  int<lower=1> R;                       // Number of room types
  vector[N] beds;                       // number of beds per roomtype
  vector[N] accom;                      // accomodates per roomtype
  vector[N] baths;                      // number of bathrooms
  int<lower=1> locat[N];                // location of the unit
  int<lower=1> room_tp[N];              // room type of the unit
  vector[N] y;                          // log of price
}

parameters {
  real eta_loc[L];           // eta location
  real eta_rt[R];            // eta room type
  real beta[4];              // beta coefficients
  real<lower=0> sigma;       // sigma of log_price
  real<lower=0> sigma_eloc;  // sigma eta localization of the unit
  real<lower=0> sigma_rt;    // sigma eta room type of the unit
}

model {
  vector[N] y_hat;

  beta ~ normal(0, 1);
  sigma ~ normal (0,1);
  sigma_eloc ~ normal (0, 1);
  sigma_rt ~ normal (0, 1);
  eta_loc ~ normal (0, sigma_eloc);
  eta_rt ~ normal (0, sigma_rt);
 
  for (i in 1:N)
      y_hat[i] = beta[1] + eta_loc[locat[i]] + eta_rt[room_tp[i]] +
                           beta[2]*beds[i] + beta[3]*accom[i] +  beta[4]*baths[i];

  y ~ normal(y_hat, sigma);

}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {  
    y_rep[n] = normal_rng(beta[1] + eta_loc[locat[n]] + eta_rt[room_tp[n]] +
                           beta[2]*beds[n] + beta[3]*accom[n] +  beta[4]*baths[n], sigma);            // Model 4-v4

    log_lik[n] = normal_lpdf( y[n] | beta[1] + eta_loc[locat[n]] + eta_rt[room_tp[n]] +
                           beta[2]*beds[n] + beta[3]*accom[n] +  beta[4]*baths[n], sigma);
  }

}
