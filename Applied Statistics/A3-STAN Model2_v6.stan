data {
  int<lower=1> N;             // number of observations
  int<lower=0,upper=1> y[N];  // switch well (0 = no; 1 = yes)
  vector[N] d;                // distance
  vector[N] a;                // arsenic
}

transformed data {
  vector[N] difdist;          // Deviance from mean Distance
  vector[N] diflogars;           // Deviance from mean log(Arsenic)
  vector[N] inter;            // Vector for interaction
  real dbar;                  // Mean of distances
  real logabar;               // Mean of log(Arsenic) concentration

  dbar = mean(d);                    // Calculate mean of distance
  logabar = mean(log(a));            // Calculate mean of log(arsenic) concentration
  difdist = d - dbar;                // Calculate deviation on Distance
  diflogars = log(a) - logabar;      // Calculate deviantion on log(Arsenic) concentration
  inter = difdist .* diflogars ;     // Define interaction between deviances of dist and log(arsenic)
}

parameters {
  real beta[4];               // Coefficients betas
}

model {
  beta ~ normal(0, 1);
  y ~ bernoulli_logit(beta[1] +
                      beta[2] * difdist +
                      beta[3] * diflogars +
                      beta[4] * inter);  // Model 2
}

generated quantities {
  vector[N] y_rep;     // replications from posterior predictive dist
  vector[N] log_lik;   // pointwise log-likelihood for LOO

  for (n in 1:N) {
    y_rep[n] = bernoulli_logit_rng(beta[1] + 
                    	  	   beta[2] * difdist[n] +
                    	  	   beta[3] * diflogars[n] +
                    	  	   beta[4] * inter[n]);       // Model 2

    log_lik[n] = bernoulli_logit_lpmf( y[n] | beta[1] + beta[2] * difdist[n] +
                                      beta[3] * diflogars[n] + beta[4] * inter[n]);
  }

}
