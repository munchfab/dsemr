// autoregressive DSEM with manifest variables
data {
  int<lower=1> N; 	// number of observational units
  int<lower=1> D; 	// number of time-varying constructs
  int<lower=1, upper=3> maxLag; // maximum lag
  int<lower=1> N_obs; 	// observations in total: N * TP
  int<lower=1> n_pars;
  int<lower=D> n_random;   // number of random effects
  int n_fixed;
  int is_fixed[1,n_fixed];
  int is_random[n_random];  // which parameters to model person-specific
  int<lower=1> N_obs_id[N]; // number of observations for each unit
  vector[N_obs] y[D]; 	    // D*N_obs array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  int n_miss_D[D];                 // missings per D
  int pos_miss_D[D,max(n_miss_D)]; // array of missings' positions

  // censored models
  real censL_val;
  int n_censL;                      // total number of obs at LB across D
  int n_censL_D[D];                 // obs at LB per D
  int pos_censL_D[D,max(n_censL_D)]; // array of obs at LBs' positions
  real censR_val;
  int n_censR;                      // total number of obs at LB across D
  int n_censR_D[D];                 // obs at LB per D
  int pos_censR_D[D,max(n_censR_D)]; // array of obs at LBs' positions

  // model adaptations based on user inputs:
  // - fixing parameters to constant values:
  // - innovation variances
  int<lower=0,upper=1> innos_rand[D];
  int n_innos_fix;
  int innos_fix_pos[D];
  int innos_pos[D];

  // - dynamic model specification per D
  int<lower=1> N_pred[D];     // Number of predictors per dimension
  int<lower=0> D_pred[D,max(N_pred)];    // matrix to index predictors to use per dimension
  int<lower=0> Lag_pred[D,max(N_pred)];  // matrix to index lag of used predictors
  int Dpos1[D];  // index positions of dynamic effect parameters
  int Dpos2[D];
  int n_int;
  int D_pred2[D,max(N_pred)];    // matrix to index predictors to use per dimension
  int Lag_pred2[D,max(N_pred)];  // matrix to index lag of used predictors

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
  matrix[n_fixed,2] prior_b_fix;
  matrix[n_innos_fix,2] prior_sigma;
  matrix[n_cov_bs,2] prior_b_re_pred;
  matrix[n_out,2] prior_alpha_out;
  matrix[n_out_bs_sum,2] prior_b_out;
  matrix[n_out,2] prior_sigma_out;

  // - covariances of innovations:
  int n_inno_covs; // number of potential innovation covs to include
  int n_obs_cov;   // total number of residuals
  int inno_cov_pos[1,n_inno_covs];
  int<lower=-1,upper=1> inno_cov_load[2];
}


parameters {
  vector[n_random] b_free[N];            // person-specific parameter
  vector<lower=0>[n_random] sd_R;        // random effect SD
  vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  cholesky_factor_corr[n_random] L;           // cholesky factor of random effects correlation matrix
  vector[n_miss] y_impute;               // vector to store imputed values
  row_vector[n_random] gammas;           // fixed effect (intercepts)
  vector[n_cov_bs] b_re_pred;            // regression coefs of RE prediction
  vector[n_fixed] b_fix;
  vector[n_out] alpha_out;               // outcome precition intercepts
  vector<lower=0>[n_out] sigma_out;      // residual SD(s) of outcome(s)
  vector[n_out_bs_sum] b_out_pred;       // regression coefs of out prediction
  vector[n_obs_cov] eta_cov[n_inno_covs];
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<lower=censR_val>[n_censR] y_impute_censR;
}

transformed parameters {
  matrix[N, n_random] bmu;     // gammas of person-specific parameters
  matrix[N,n_pars] b;
  vector[N] sd_noise[D];
  vector[N] sd_inncov[n_inno_covs];
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
  if(n_fixed>0){
    for(i in 1:n_fixed){
      b[,is_fixed[1,i]] = rep_vector(b_fix[i],N);
    }
  }

  // transformation of log-innovation variances if modeled as person-specific
  for(i in 1:D){
    if(innos_rand[i] == 0){
      sd_noise[i,] = rep_vector(sigma[innos_fix_pos[i]],N);
    } else {
      sd_noise[i,] = sqrt(exp(b[,innos_pos[i]]));
    }
  }
  // transform log innovation covarainces
  if(n_inno_covs > 0){
      for(i in 1:n_inno_covs){
        sd_inncov[i,1:N] = sqrt(exp(to_vector(b[,inno_cov_pos[1,i]])));
      }
    }
 }

model {
  int pos = 1;       // initialize position indicator
  int p_miss = 1;    // running counter variable to index positions on y_impute
  int p_censL = 1;
  int p_censR = 1;
  int pos_cov = 1;   // covariance position
  int obs_id = 1;    // declare local variable to store variable number of obs per person
  vector[N_obs] y_merge[D];
  matrix[n_random, n_random] SIGMA = diag_pre_multiply(sd_R, L);

  y_merge = y;      // add observations
  if(n_miss>0){
    for(i in 1:D){
    // add imputed values for missings on each indicator
    y_merge[i,pos_miss_D[i,1:n_miss_D[i]]] = segment(y_impute, p_miss, n_miss_D[i]);
    p_miss = p_miss + n_miss_D[i];    // update counter for next indicator i+1
  }
  }
  // replace values at censor thresholds
  if(n_censL>0){
    for(i in 1:D){
    // add imputed values for observations at floor (threshold for censoring)
    y_merge[i,pos_censL_D[i,1:n_censL_D[i]]] = segment(y_impute_censL, p_censL, n_censL_D[i]);
    p_censL = p_censL + n_censL_D[i];    // update counter for next indicator i+1
    }
  }
  if(n_censR>0){
    for(i in 1:D){
    // add imputed values for observations at ceiling (threshold for censoring)
    y_merge[i,pos_censR_D[i,1:n_censR_D[i]]] = segment(y_impute_censR, p_censR, n_censR_D[i]);
    p_censR = p_censR + n_censR_D[i];    // update counter for next indicator i+1
    }
  }

  // (Hyper-)Priors
  gammas ~ normal(prior_gamma[,1],prior_gamma[,2]);
  sd_R ~ cauchy(prior_sd_R[,1], prior_sd_R[,2]);
  L ~ lkj_corr_cholesky(prior_LKJ);

  if(n_innos_fix>0){
    sigma ~ cauchy(prior_sigma[,1], prior_sigma[,2]);
  }

  if(n_cov > 1){
    b_re_pred ~ normal(prior_b_re_pred[,1], prior_b_re_pred[,2]);
  }
  if(n_out > 0){
    alpha_out ~ normal(prior_alpha_out[,1], prior_alpha_out[,2]);
    b_out_pred ~ normal(prior_b_out[,1], prior_b_out[,2]);
    sigma_out ~ cauchy(prior_sigma_out[,1], prior_sigma_out[,2]);
  }

  if(n_fixed > 0){
    b_fix ~ normal(prior_b_fix[,1],prior_b_fix[,2]);
  }

  for (pp in 1:N) {
    // store number of observations per person
    obs_id = (N_obs_id[pp]);

    // individual parameters from (multivariate) normal distribution
    if(n_random == 1){
      b_free[pp,1] ~ normal(bmu[pp,1], sd_R[1]);
    } else {
      b_free[pp, 1:n_random] ~ multi_normal_cholesky(bmu[pp, 1 : n_random], SIGMA);
    }

    vector[obs_id-maxLag] eta_cov_id[n_inno_covs];
    if(n_inno_covs>0){
      for(i in 1:n_inno_covs){
          eta_cov_id[i,] = segment(eta_cov[i,], pos_cov, (obs_id-maxLag));
          eta_cov_id[i,] ~ normal(0, sd_inncov[i,pp]);
      }
    }

    // local variable declaration: array of predicted values
    {
    vector[obs_id-maxLag] mus[D];

    // create latent mean centered versions of observations
    vector[obs_id] y_cen[D];

    for(d in 1:D){ // start loop over dimensions
      y_cen[d,] = y_merge[d,pos:(pos+obs_id-1)] - b[pp,d];
    }

    for(d in 1:D){ // start loop over dimensions
      // build prediction matrix for specific dimensions
      int n_cols; // matrix dimensions
      n_cols = (n_inno_covs>0 && d<3) ? N_pred[d]+n_inno_covs : N_pred[d];

      {
      matrix[(obs_id-maxLag),n_cols] b_mat;
      vector[n_cols] b_use;
      for(nd in 1:N_pred[d]){ // start loop over number of predictors in each dimension
         int lag_use = Lag_pred[d,nd];
         if(D_pred2[d,nd] == -99){
          b_mat[,nd] = y_cen[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)];
         } else {
          int lag_use2 = Lag_pred2[d,nd];
          b_mat[,nd] = y_cen[D_pred[d, nd],(1+maxLag-lag_use):(obs_id-lag_use)] .*
                       y_cen[D_pred2[d, nd],(1+maxLag-lag_use2):(obs_id-lag_use2)];
         }
      }
      b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);

      if(n_inno_covs>0 && d < 3){  // add latent factor scores
         b_mat[,(N_pred[d]+1)] = eta_cov_id[1,]; // add innovation covariance factor scores
         b_use[N_pred[d]+1] = inno_cov_load[d];
         }

      mus[d,] = b[pp,d] + b_mat * b_use;
      }

      // sampling statement
      y_merge[d,(pos+maxLag):(pos+(obs_id-1))] ~ normal(mus[d,], sd_noise[d,pp]);
      }
    }
    // update index variables
    pos = pos + obs_id;
    pos_cov = pos_cov + obs_id - maxLag;
  } // end loop over subjects

  // outcome prediction: get expectations of outcome values
  if(n_out > 0){
    int k = 1;
    matrix[N,n_random+n_z] b_z = append_col(b[,is_random],Z);
    for(i in 1:n_out){
      int n_bs = n_out_bs[i,1];      // number of predictors for each outcome
      out[i,] ~ normal(alpha_out[i] + b_z[,n_out_b_pos[i,1:n_bs]] * segment(b_out_pred,k,n_bs), sigma_out[i]);
      k = k + n_bs; // update index
    }
  }
}


generated quantities{
  matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
  bcorr = multiply_lower_tri_self_transpose(L);
}
