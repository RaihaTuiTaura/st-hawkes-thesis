
##############################################################

# 2D index to 1D index
convert_dist <- function (x, y, dist_size) {
  n = dist_size

  k_all = sapply(1:length(y), function(l){
    y_l = as.numeric(y[l])
    i = max(x,y_l)
    j = min(x,y_l)
    valid = (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
    k = (2 * n - j) * (j - 1) / 2 + (i - j)
    k[!valid] = NA_real_
    k
  })
  return(k_all)
}

mcsapply<-function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.preschedule = TRUE,
                    mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                    mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL )
{
  answer <- mclapply(X = X, FUN = FUN, ...,mc.preschedule = mc.preschedule,
                     mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores,
                     mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive, affinity.list = affinity.list)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

decay_rbf = function(dists,sigma){
  exp(-dists^2/(2*sigma^2))
}

decay_geometric = function(t,h,p){
  h$beta[p] * (1-h$beta[p])^(t-1)
}


#############################################################

# Calculate conditional intensity function
lambda_cond_multidim_dt_st = function(theta,decay_fn,decay_space_fn,data){
  M = length(theta$mu)
  f=function(m){
    y_m = lambda_cond_multidim_dt_st_m(theta,m,decay_fn,decay_space_fn,data)
    return(y_m)}
  L=lapply(1:M,f)
  return(L)
}

# For each combination of variates
lambda_cond_multidim_dt_st_m = function(theta,m,decay_fn,decay_space_fn,data){

  h = theta$h
  mu = theta$mu
  sigma = theta$sigma

  M =length(mu)
  decay_dists_vals  = lapply(1:M, function(l) {p=(m-1)*M+l;
    t(sapply(1:nrow(data$space_dists[[l]]), function(k) {
      decay_space_fn(data$space_dists[[l]][k,],sigma[p]) * data$Times_events_index[[l]]$n
    }))
  })

  g=function(l){
    p = (m-1)*M+l
    decay_space_vals =  unlist(lapply(data$Times_obs, function(t){
      t_ind = which(data$Times_events_index[[l]]$times %in% (t-data$smax[p]):(t-1))
      t_ind_times = data$Times_events_index[[l]]$times[t_ind]
      sapply(1:nrow(decay_dists_vals[[l]]), function(k) {
        sum(decay_dists_vals[[l]][k,t_ind] * decay_fn((t-t_ind_times),h,p),na.rm=T)
      })
    }))
    L_m = h$alpha[p]*decay_space_vals
    return(L_m)
  }

  L_m=lapply(1:M,g)
  L_m=matrix(unlist(L_m), nrow=length(unlist(L_m))/M, ncol=M)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m=mu[m] + R_m
  return(y_m)

}

# Calculate conditional intensity function
lambda_cond_multidim_dt = function(theta,decay_fn,data){
  M = length(theta$mu)
  f=function(m){
    y_m = lambda_cond_multidim_dt_m(theta,m,decay_fn,data)
    return(y_m)}
  L=lapply(1:M,f)
  return(L)
}

# For each combination of variates
lambda_cond_multidim_dt_m = function(theta,m,decay_fn,data){

  h = theta$h
  mu = theta$mu

  M =length(mu)

  g=function(l){
    p = (m-1)*M+l
    decay_vals =  sapply(data$Times_obs, function(t){
      t_ind = which(data$Times_events[[l]]$times %in% (t-data$smax[p]):(t-1))
      t_ind_times = data$Times_events[[l]]$times[t_ind]
      t_ind_n = data$Times_events[[l]]$n[t_ind]
      sum(decay_fn((t-t_ind_times),h,p) * t_ind_n, na.rm=T)
    })
    L_m = h$alpha[p]*decay_vals
    return(L_m)
  }

  L_m=lapply(1:M,g)
  L_m=matrix(unlist(L_m), nrow=length(unlist(L_m))/M, ncol=M)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m=mu[m] + R_m
  return(y_m)

}


#############################################################

# Log likelihood
log_likelihood_multidim_dt_st=function(data,decay_fn,decay_space_fn,theta_vec,theta=NULL,op='sum'){

  #for optim (geometric only)
  theta = list(mu=theta_vec[1:M], h=list(alpha=theta_vec[(M+1):(M+M^2)], beta=theta_vec[(M+M^2+1):(M+2*M^2)]), sigma=theta_vec[((M+2*M^2)+1):((M+3*M^2))])
  
  M = length(theta$mu)
  eval_lambda = lambda_cond_multidim_dt_st(theta,decay_fn,decay_space_fn,data)
  L= vapply(1:M,function(m){sum(dpois(data$Times_all[[m]],lambda=eval_lambda[[m]],log=TRUE))},1)

  if (min(L)==-Inf | sum(is.na(L))>0){L=-10000}
  if (op=='sum'){L = sum(L)};
  return(L)
}


# Log likelihood
log_likelihood_multidim_dt=function(data,decay_fn,theta_vec,theta=NULL,op='sum'){

  #for optim (geometric only)
  theta = list(mu=theta_vec[1:M], h=list(alpha=theta_vec[(M+1):(M+M^2)], beta=theta_vec[(M+M^2+1):(M+2*M^2)]))

  M = length(theta$mu)
  eval_lambda = lambda_cond_multidim_dt(theta,decay_fn,data)
  L= vapply(1:M,function(m){sum(dpois(data$Times_all[[m]],lambda=eval_lambda[[m]],log=TRUE))},1)

  if (min(L)==-Inf | sum(is.na(L))>0){L=-10000}
  if (op=='sum'){L = sum(L)};
  return(L)
}
