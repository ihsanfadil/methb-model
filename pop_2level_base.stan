
functions {
  vector methb_model(
    
    int ndays,
    vector dose, 
    real logit_total_rate, 
    real alpha

    ){
    vector[ndays] methb_total; 
    real steady_state = 0.01;
    real total_rate = inv_logit(logit_total_rate);
    real prod_rate = steady_state * total_rate;
    real red_rate = (1 - steady_state) * total_rate;

    methb_total[1] = steady_state;

    for(t in 2:ndays) {
      real dose_effect = exp(alpha * dose[t-1]);
      methb_total[t] = methb_total[t-1]
                       + ( prod_rate * dose_effect * ( 1 - methb_total[t-1] ) )
                       - ( red_rate * methb_total[t-1] );
    }
    return methb_total;
  }
}

data {
  int<lower = 1> n_patients;
  int ndays[n_patients];
  int<lower = 1> nobs[n_patients];
  vector<lower = 0, upper = 1>[sum(nobs)] methb;
  vector<lower = 0>[sum(ndays)] dose;
  int<lower = 1, upper = sum(ndays)> ind_obs[sum(nobs)];
  int<lower = 1, upper = n_patients> patient_id[sum(nobs)];

  // Covariates
  int<lower = 1> sex[n_patients]; // Male, female
  int<lower = 1> age_cat[n_patients]; // Age categories
  int<lower = 1> site[n_patients]; // Study sites
  int<lower = 1, upper = 3> dose_group[n_patients]; // Dose categories: Low, intermediate, high

  // Sizes
  int<lower = 1> n_sites;
  int<lower = 1> n_age_cats;
}

parameters {
  real logit_total_rate;
  real logit_mu_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma;
  vector[n_patients] log_alpha_p;
}

transformed parameters {
  vector<lower = 0>[sum(ndays)] pred_methb;
  {
    int pos = 1;
    for (p in 1:n_patients) {
      vector[ndays[p]] temp_methb;
      real alpha;
      alpha = inv_logit(logit_mu_alpha) * exp(log_alpha_p[p]); 
      temp_methb = methb_model(ndays[p],
                               dose[pos:pos+ndays[p]-1],
                               logit_total_rate,
                               alpha);
      pred_methb[pos:pos+ndays[p]-1] = temp_methb;
      pos += ndays[p];
    }
  }
}

model {
  // Priors
  logit_mu_alpha ~ normal(3, 1);
  logit_total_rate ~ normal(-4, 1);
  sigma_alpha ~ normal(1, .5) T[0,];
  sigma ~ normal(.01, .005) T[0,];
  log_alpha_p ~ normal(0, sigma_alpha);

  // Likelihood
  methb ~ normal(pred_methb[ind_obs], sigma);
}

generated quantities {
  vector[sum(nobs)] log_lik;
  vector[sum(nobs)] residuals;
  vector[n_patients] alpha_individual;
  
  for (i in 1:sum(nobs)) {
    real pred = pred_methb[ind_obs[i]];
    log_lik[i] = normal_lpdf(methb[i] | pred_methb[ind_obs[i]], sigma);
    residuals[i] = methb[i] - pred;
  }
  for (p in 1:n_patients) {
    alpha_individual[p] = inv_logit(logit_mu_alpha) * exp(log_alpha_p[p]);
  }
}
