// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	// number of observational units
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=2> n_random;   // number of random effects
  int n_fixed;
  int is_fixed[1,n_fixed];
  int is_random[n_random]; // which parameters to model person-specific
  int<lower=1> N_obs_id[N]; // number of observations for each unit
  vector[N_obs] y; 	// D*N_obs array of observations

  // handling of missing values
  int n_miss;             // total number of missings across D
  int pos_miss[n_miss];   // array of missings' positions

  // model adaptations based on user inputs:
  // - fixing parameters to constant values
  // - 1. innovation variances
  int<lower=0,upper=1> innos_rand;

  // - time-invariant variables:
  // covariates as predictors of random effects
  int<lower=1> n_cov;           // number of covariates - minimum of 1 for intercepts
  int n_cov_bs;
  int n_cov_mat[n_cov_bs, 2];
  matrix[N, n_cov] W;  // predictors of individual parameters
  // outcome prediction
  int n_out;                 // number of outcome variables
  int n_out_bs[n_out,1];     // number of predictors per outcome
  int n_out_bs_max;          // number of predictors per outcome
  int n_out_bs_sum;          // number of predictors per outcome
  int n_out_b_pos[n_out,n_out_bs_max]; // index positions
  int n_z;              // number of additional time-invariant as outcome predictors
  matrix[N, n_z] Z;     // observations of Z
  vector[N] out[n_out];        // outcome

  // priors
  matrix[n_random,2] prior_gamma;
  matrix[n_random,2] prior_sd_R;
  real prior_LKJ;
  matrix[1-innos_rand,2] prior_sigma;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;
}

parameters {
  vector[n_random] b_free[N];            // person-specific parameter
  vector<lower=0>[n_random] sd_R;        // random effect SD
  vector<lower=0>[1-innos_rand] sigma;   // SDs of fixed innovation variances
  cholesky_factor_corr[n_random] L;      // cholesky factor of random effects correlation matrix
  vector[n_miss] y_impute;               // vector to store imputed values
  row_vector[n_random] gammas;           // fixed effect (intercepts)
  vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  vector[n_out] alpha_out;
  vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  vector[n_out_bs_sum] b_out_pred;           // regression coefs of out prediction
}

transformed parameters {
  matrix[N, n_random] bmu;               // gammas of person-specific parameters
  matrix[N,n_pars] b;
  vector<lower = 0>[N] sd_noise;
  matrix[n_cov, n_random] b_re_pred_mat = rep_matrix(0, n_cov, n_random);

  // REs regressed on covariates
  b_re_pred_mat[1,] = gammas;
  if(n_cov>1){
     for(i in 1:n_cov_bs){
     b_re_pred_mat[n_cov_mat[i,1],n_cov_mat[i,2]] = b_re_pred[i];
    }
  }
  // calculate population means (intercepts) of person-specific parameters
  bmu = W * b_re_pred_mat;

  // create array of (person-specific) parameters to use in model
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }

  // transformation of log-innovation variances if modeled as person-specific
  if(innos_rand == 0){
      sd_noise = rep_vector(sigma[1],N);
    } else {
      sd_noise = sqrt(exp(b[,(n_random)]));
    }
}


model {
  int pos = 1;       // initialize position indicator
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sd_R, L);
  vector[N_obs] y_merge;
  y_merge = y;              // add observations

  // add imputed values for missings on each indicator
  if(n_miss>0){
    y_merge[pos_miss] = y_impute;
  }


  // (Hyper-)Priors
  gammas ~ normal(prior_gamma[,1],prior_gamma[,2]);
  sd_R ~ cauchy(prior_sd_R[,1], prior_sd_R[,2]);
  L ~ lkj_corr_cholesky(prior_LKJ);

  if(innos_rand == 0){
    sigma ~ cauchy(prior_sigma[1,1], prior_sigma[1,2]);
  }
  if(n_cov > 1){
    b_re_pred ~ normal(prior_b_re_pred[,1], prior_b_re_pred[,2]);
  }
  if(n_out > 0){
    alpha_out ~ normal(prior_alpha_out[,1], prior_alpha_out[,2]);
    b_out_pred ~ normal(prior_b_out[,1], prior_b_out[,2]);
    sigma_out ~ cauchy(prior_sigma_out[,1], prior_sigma_out[,2]);
  }


  for (pp in 1:N) {
    // store number of observations per person
    obs_id = N_obs_id[pp];

    // individual parameters from multivariate normal distribution
    b_free[pp, is_random] ~ multi_normal_cholesky(bmu[pp, 1 : n_random], SIGMA);

    // local variable declaration: array of predicted values
    {
    vector[obs_id-1] mus;

    // create latent mean centered versions of observations
    vector[obs_id] y_cen;
    y_cen = y_merge[pos:(pos+obs_id-1)] - b[pp,1];

    // use build predictor matrix to calculate latent time-series means
    mus =  b[pp,1] + b[pp,2] * y_cen[1:(obs_id-1)];

    // sampling statement
    y_merge[(pos+1):(pos+(obs_id-1))] ~ normal(mus, sd_noise[pp]);
    }

    // update index variables
    pos = pos + obs_id;
  }

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    int k = 1;
    matrix[N,n_random+n_z] b_z = append_col(b,Z);
    for(i in 1:n_out){
      int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
      out[i,] ~ normal(alpha_out[i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred,k,n_bs), sigma_out[i]);
      k = k + n_bs; // update index
    }
  }
}

generated quantities{
  matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
//  matrix[n_random,n_random] bcov; // random coefficients covariance matrix
  // create random coefficient matrices
  bcorr = multiply_lower_tri_self_transpose(L);
//  bcov = quad_form_diag(bcorr, sd_R[1:n_random]);
}
