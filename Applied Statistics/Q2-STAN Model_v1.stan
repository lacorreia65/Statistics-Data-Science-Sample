data {
  int<lower=1> N;                       // number of observations
  int<lower=1> C;                       // Number of countries
  int<lower=1> R;                       // Number of regions
  vector[N] log_gdp;                    // log of GDP of country
  vector[N] log_gfr;                    // log of GFR of country
  vector[N] sab;                        // SAB of that country
  int<lower=1> country[N];              // country of observation
  int<lower=1> region[N];               // region of observation
  vector[N] y;                          // log of PMNA
}

parameters {
  real eta_country[C];           // eta country
  real eta_region[R];            // eta region
  real beta[4];                  // beta coefficients
  real<lower=0> sigma;           // sigma of log_pmna
  real<lower=0> sigma_country;   // sigma eta country
  real<lower=0> sigma_region;    // sigma eta region
}

model {
  vector[N] y_hat;

  beta ~ normal(0, 1);
  sigma ~ normal (0,1);
  sigma_country ~ normal (0, 1);
  sigma_region ~ normal (0, 1);
  eta_country ~ normal (0, sigma_country);
  eta_region~ normal (0, sigma_region);
 
  for (i in 1:N)
      y_hat[i] = beta[1] + eta_country[country[i]] + eta_region[region[i]] +
                           beta[2]*log_gdp[i] + beta[3]*log_gfr[i] +  beta[4]*sab[i];

  y ~ normal(y_hat, sigma);

}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {  
    y_rep[n] = normal_rng(beta[1] + eta_country[country[n]] + eta_region[region[n]] +
                            beta[2]*log_gdp[n] + beta[3]*log_gfr[n] +  beta[4]*sab[n], sigma);            // Model 4-v4

    log_lik[n] = normal_lpdf( y[n] | beta[1] + eta_country[country[n]] + eta_region[region[n]] +
                                       beta[2]*log_gdp[n] + beta[3]*log_gfr[n] +  beta[4]*sab[n], sigma);
  }

}
