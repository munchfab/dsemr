#' Prepare data for Stan
#'
#' @param model `data.frame`. Output of model-Functions.
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param ts Character. The variable(s) in `data` that
#' contain the time-series construct. If multiple variables are provided in a
#' character vector, a vector autoregressive model is fit.
#' @param covariates Named character vector. An optional named vector of
#' characters to refer to predictors of random effects as specified in the `model`.
#' Note that specifying `covariates` is only necessary if the respective
#' variable name(s) in `data` differ from the variables names specified in `model`.
#' @param outcomes Named character vector. Similar to `covariates`, an optional named vector of
#' characters to refer to outcome predicted by random effects as specified in the `model`.
#' Note that specifying `outcomes` is only necessary if the respective
#' variable name(s) in `data` differ from the outcome variable name(s) specified in `model`
#' @param outcome_pred_btw tba.
#' @param center_covs tba.
#' @param complete_obs tba.
#'
#' @return A `list` object that can be passed to \code{\link[rstan]{sampling}}.
#' @export
#' @noRd
#'
#' @details
#' For internal use only.
#'
#' @examples
#' \donttest{
#' # build simple vector-autoregressive model
#' var_model <- mlts_model(q = 2)
#'
#' # pass data and model to VARprepare()
#' stan_list <- VARprepare(
#'   model = var_model,
#'   data = ts_data,
#'   ts = c("Y1", "Y2")
#' )
#' }
#'
#'
VARprepare <- function(model, data, ts, covariates = NULL, outcomes = NULL,
                       outcome_pred_btw = NULL, center_covs = T, complete_obs = F
){


  # data entered here should be preprocessed prior (e.g., handling of missing values...)

  # data specific information: -------------------------------------------------
  N = length(unique(data$num_id))         # number of subjects
  D = length(ts)                          # number of time-varying constructs
  N_obs = nrow(data)                      # total number of observations (obs)
  N_obs_id = data.frame(table(data$num_id))$Freq # number of obs per subject

  ## handling of missing values
  n_miss = sum(is.na(data[,ts]))           # overall number of NAs
  n_miss_D = sapply(ts, FUN = function(x){ # store number of NAs per D
    sum(is.na(data[,x]))
  })
  n_miss_D = as.array(n_miss_D)
  # NA positions as matrix
  pos_miss_D = matrix(0, nrow = D, ncol = max(n_miss_D), byrow = T)
  for(i in 1:D){
    if(n_miss_D[i] > 0){
      pos_miss_D[i,1:n_miss_D[i]] = which(is.na(data[,ts[i]]))
    }
  }

  # model specific information: ------------------------------------------------
  infos = mlts_model_eval(model)
  maxLag = infos$maxLag


  # alternative handling of missing values via indexing -------------------------
  if(complete_obs == T){

  #### Experimental at this point....
  ids = unique(data$num_id)
  # add row id to data
  data$row_id = 1:nrow(data)
  N_obs_id_new = c()
  N_obs_id_use = c()
  # get number of obs to use as DV
  for(p in 1:N){
    dat_id = data[data$num_id == p,c("num_id","row_id",ts)]
    # create lagged variables
    for(d in 1:length(ts)){
      for(i in 1:maxLag){
        dat_id[,paste(ts[d],"__lag",i)] = c(rep(NA,times = i),dat_id[,ts[d]][1:(nrow(dat_id)-i)])
      }
    }
    # check for complete rows
    if(length(ts) == 1){
      row_y_complete = which(!is.na(dat_id[,3:(2+length(ts))]))
    } else {
      row_y_complete = which(rowSums(is.na(dat_id[,3:(2+length(ts))])) == 0)
    }
    row_use_id = which(rowSums(is.na(dat_id[,3:ncol(dat_id)])) == 0)

    # get lagged obs
    rows_lag = c(row_use_id)
    for(z in 1:maxLag){
      rows_lag = c(rows_lag, (row_use_id-z))
    }
    # filter complete rows with complete lagged observations
    row_y_complete = row_y_complete[row_y_complete %in% rows_lag]

    N_obs_id_use[p] = length(row_use_id)
    N_obs_id_new[p] = length(row_y_complete)
  }

  # get position indices on Y in the data
  N_obs_use_max = max(N_obs_id_use)
  N_obs_id_max = max(N_obs_id_new)
  N_obs_y_pos = array(data = 0, dim = c(N, N_obs_id_max))
  N_obs_yt_pos = array(data = 0, dim = c(N, N_obs_use_max))
  N_obs_ylag_pos = array(data = 0, dim = c(N, maxLag, N_obs_use_max))

  for(p in 1:N){
    dat_id = data[data$num_id == p,c("num_id","row_id", ts)]
    # create lagged variables
    for(d in 1:length(ts)){
      for(i in 1:maxLag){
        dat_id[,paste(ts[d],"__lag",i)] = c(rep(NA,times = i),dat_id[,ts[d]][1:(nrow(dat_id)-i)])
      }
    }
    # check for complete rows
    if(length(ts) == 1){
      row_y_complete = which(!is.na(dat_id[,3:(2+length(ts))]))
    } else {
      row_y_complete = which(rowSums(is.na(dat_id[,3:(2+length(ts))])) == 0)
    }
    row_use_id = which(rowSums(is.na(dat_id[,3:ncol(dat_id)])) == 0)

    # get lagged obs
    rows_lag = c(row_use_id)
    for(z in 1:maxLag){
      rows_lag = c(rows_lag, (row_use_id-z))
    }
    # filter complete rows with complete lagged observations
    row_y_complete = row_y_complete[row_y_complete %in% rows_lag]

    dat_id_cen = dat_id[row_y_complete,]
    dat_id_cen$row_id_cen = 1:nrow(dat_id_cen)

    # get positions of Yt
    N_obs_y_pos[p, 1:N_obs_id_new[p]] = dat_id$row_id[row_y_complete]
    N_obs_yt_pos[p, 1:N_obs_id_use[p]] = dat_id$row_id[row_use_id]
    for(z in 1:maxLag){
      rows = dat_id_cen$row_id_cen[dat_id_cen$row_id %in% (N_obs_yt_pos[p, 1:N_obs_id_use[p]] - z)]
      N_obs_ylag_pos[p,z,1:N_obs_id_use[p]] = rows
    }
   }
  n_obs_cov_new = sum(N_obs_id_use)

  } else {
    N_obs_id_new = 0
    N_obs_id_max = 0
    N_obs_id_use = 0
    N_obs_use_max = 0
    N_obs_y_pos = 0
    N_obs_yt_pos = 0
    N_obs_ylag_pos = 0
    n_obs_cov_new = 0
  }

  # END EXPERIMENTAL SECTION -----------------------------------

  ## replace missing values with -Inf
  data[,ts][is.na(data[,ts])] = -Inf
  y = t(data[,ts])   # store observations

  # ----

  n_pars = infos$n_pars       # no of dynamic parameters (including means, CRs, innovation variance)
  n_random = infos$n_random   # no of individual (random) effects
  n_fixed = infos$n_fixed     # no of fixed dynamic parameters (AR and CR effects)
  is_random = as.array(infos$is_random) # which parameters (in order of Fixed effect parameters in model) to model as random effect
  is_fixed = infos$is_fixed   # a matrix of n_fixed x 1, indicating which parameter to model as constant
  re_pars = infos$re_pars     # subset of model of random effect parameters

  # innovations
  innos_rand = as.array(infos$innos_rand) # Are innovation variance(s) random?
  n_innos_fix = infos$n_innos_fix
  innos_fix_pos = as.array(infos$innos_fix_pos)
  innos_pos = as.array(infos$innos_pos)

  # innovation covariances
  n_inno_covs = infos$n_inno_covs
  n_inno_cors = infos$n_inno_cors
  n_inno_cov_fix = infos$n_inno_cov_fix
  inno_cov_pos = infos$inno_cov_pos
  inno_cov_load = infos$inno_cov_load

  # NEEDS TO BE UPDATED ----------------------------
  n_obs_cov = ifelse(n_inno_covs == 0, 0, N_obs - (N*infos$maxLag))
  # -----------------------------------------


  # dynamic model specification by dimension (D)
  maxLag = infos$maxLag
  N_pred = infos$N_pred
  D_pred = infos$D_pred
  Lag_pred = infos$Lag_pred
  Dpos1 = as.array(infos$Dpos1)
  Dpos2 = as.array(infos$Dpos2)
  # covariate(s) as predictor(s) of random effects
  RE.PREDS = infos$RE.PREDS
  n_cov = infos$n_cov
  n_cov_bs = infos$n_cov_bs
  n_cov_mat = infos$n_cov_mat
  W = matrix(NA, nrow = N, ncol = n_cov) # build matrix of predictors
  W[,1] = 1
  if(n_cov > 1){
    for(i in 2:n_cov){
      cov.use = which(covariates == unique(RE.PREDS$re_preds[RE.PREDS$pred_no == i-1]))
      for(p in 1:N){
        W[p,i] = unique(data[data$num_id == p, names(covariates)[cov.use]])
      }
      if(center_covs==T){
        W[,i] = W[,i] - mean(W[,i])
      }
    }
  }
  # outcome prediction - REs as predictors
  n_out = infos$n_out
  out_vars = infos$out_var
  out_vars_data = unlist(lapply(out_vars, function(x){
    names(outcomes)[which(outcomes == x)]
  }))
  n_out_bs = infos$n_out_bs
  n_out_b_pos = infos$n_out_b_pos
  n_out_bs_max = infos$n_out_bs_max
  n_out_bs_sum = infos$n_out_bs_sum

  # outcome data
  out = matrix(nrow = n_out, ncol = N, data = NA)
  if(n_out>0){
    for(i in 1:n_out){
      for(p in 1:N){
        out[i,p] = unique(data[data$num_id==p,out_vars_data[i]])
      }
    }
  }

  n_z_vars = infos$n_z_vars
  n_z_vars_data = unlist(lapply(n_z_vars, function(x){
    names(outcome_pred_btw)[which(outcome_pred_btw == x)]
  }))
  n_z = infos$n_z
  Z = matrix(NA, nrow = N, ncol = n_z)
  if(n_z > 0){
    for(i in 1:n_z){
      for(p in 1:N){
        Z[p,i] = unique(data[data$num_id==p, n_z_vars_data[i]])
      }
    }
  }

  # priors
  prior_gamma = infos$prior_gamma
  prior_sd_R = infos$prior_sd_R
  prior_LKJ = infos$prior_LKJ
  prior_sigma = infos$prior_sigma
  prior_b_re_pred = infos$prior_b_re_pred
  prior_b_out = infos$prior_b_out
  prior_alpha_out = infos$prior_alpha_out
  prior_sigma_out = infos$prior_sigma_out
  prior_b_fix = infos$prior_b_fix

  # get SDs of latent variables
  standardized = 0   # will be overwritten in mlts_fit if requested

  # combine all information
  standata = rstan::nlist(
    # data specifications
    N, D,maxLag, N_obs, N_obs_id, y, n_miss, n_miss_D, pos_miss_D,
    # addons for faster missing data appoach via indexing
    N_obs_id_new, N_obs_id_max, N_obs_id_use, N_obs_use_max,
    N_obs_y_pos, N_obs_yt_pos,N_obs_ylag_pos, n_obs_cov_new,

    # model specifications
    n_pars, n_random, n_fixed, is_random, is_fixed,
    N_pred, D_pred, Lag_pred, Dpos1, Dpos2,
    ## innovations
    innos_rand, innos_pos, n_innos_fix,
    ## innovation covariances
    innos_fix_pos, n_obs_cov, n_obs_cov, inno_cov_pos,
    n_inno_covs,n_inno_cors, n_inno_cov_fix, inno_cov_load,

    n_cov, n_cov_bs, n_cov_mat, W,
    n_out, n_out_bs, n_out_bs_max, n_out_bs_sum, n_out_b_pos, out, n_z, Z,
    # priors
    prior_LKJ, prior_gamma, prior_sd_R, prior_sigma, prior_b_re_pred,
    prior_b_fix,
    prior_b_out, prior_alpha_out, prior_sigma_out,

    # add indicator names to print in summary
    standardized, ts
  )

  # for latent models:
  if(infos$isLatent == TRUE){
    standata$D = infos$q
    standata$D_np = as.array(infos$p)
    standata$n_p = nrow(infos$indicators)
    standata$is_SI = as.array(unlist(lapply(X = 1:nrow(infos$indicators), FUN = function(x){
      ifelse(sum(infos$indicators$q == infos$indicators$q[x]) == 1, 1, 0)
      })))
    standata$D_pos_is_SI = as.array(unlist(lapply(X = 1:standata$D, FUN = function(x){
      ifelse(infos$p[x] == 1, which(infos$indicators$q == x), 0)
    })))

    standata$D_perP = as.array(as.integer(infos$indicators$q))
    standata$n_etaW_free = sum(infos$p > 1)

    # alternative missing data indexing
    n_miss = sum(data[,ts] == -Inf)           # overall number of NAs
    n_miss_p = as.array(sapply(ts, FUN = function(x){sum(data[,x] == -Inf)}))
    pos_miss_p = matrix(0, nrow = standata$n_p, ncol = max(n_miss_p), byrow = T)
    for(i in 1:standata$n_p){
      if(n_miss_p[i] > 0){
        pos_miss_p[i,1:n_miss_p[i]] = which(data[,ts[i]] == -Inf)
      }
    }
    standata$n_miss <- n_miss
    standata$n_miss_p <- n_miss_p
    standata$pos_miss_p <- pos_miss_p

    standata$n_loadBfree <- infos$n_loadBfree
    standata$n_loadWfree <- infos$n_loadWfree
    standata$n_alphafree <- infos$n_alphafree
    standata$n_sigmaBfree <- infos$n_sigmaBfree
    standata$n_sigmaWfree <- infos$n_sigmaWfree
    standata$pos_loadBfree <- as.array(infos$pos_loadBfree)
    standata$pos_loadWfree <- as.array(infos$pos_loadWfree)
    standata$pos_alphafree <- as.array(infos$pos_alphafree)
    standata$pos_sigmaBfree <- as.array(infos$pos_sigmaBfree)
    standata$pos_sigmaWfree <- as.array(infos$pos_sigmaWfree)

    standata$n_YB_free <- infos$n_YB_free
    standata$YB_free_pos <- as.array(infos$YB_free_pos)
    standata$mu_etaB_pos <- as.array(infos$mu_etaB_pos)
    standata$mu_is_etaB <- as.array(infos$mu_is_etaB)

    standata$prior_alpha = infos$prior_alpha
    standata$prior_loadB = infos$prior_loadB
    standata$prior_loadW = infos$prior_loadW
    standata$prior_sigmaB = infos$prior_sigmaB
    standata$prior_sigmaW = infos$prior_sigmaW
  }


  return(standata)

}
