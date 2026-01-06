
##############################################################

# 2D index to 1D index
convert_dist <- function (x, y, dist_size) {
  n = dist_size
  k_all = sapply(1:nrow(y), function(l){
    y_l = as.numeric(y[l,])
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
 
decay_powerlaw = function(x, pars){
  p=pars[1]
  c=pars[2]
  return((p-1)/c*(1+(x/c))^-p)
}


decay_rbf = function(dists,pars){
  exp(-dists^2/(2*pars^2))
}

decay_geometric = function(t,pars){
  pars * (1-pars)^(t-1)
}


calc_bic = function(loglik, n_pars){
  -2*loglik + log(n_pars)*n_pars
}

#############################################################


################### no spatial or temporal SE + spatially varying mu


lambda_cond_multidim_dt_icar_only = function(theta,decay_fn,data,n_locations){
  
  h = theta$h
  mu_components = theta$mu
  mu = exp(mu_components$mu_intercept + mu_components$phi)
  
  M =1
  
  y_m=matrix(mu,nrow=(n_locations*length(data$Times_obs)),ncol=M)
  return(y_m)
  
}

# Log likelihood
log_likelihood_dt_icar_only=function(data,decay_fn,theta_vec,n_locations,op='sum'){
  
  
  #for optim
  theta = list(mu=list(mu_intercept=theta_vec[1],
                       phi=theta_vec[2:(2+n_locations-1)]),
               prop_0=theta_vec[(2+n_locations)])
  
  
  M = 1
  eval_lambda = lambda_cond_multidim_dt_icar_only(theta,decay_fn,data,n_locations)
  L= vapply(1:M,function(m){sum(dzipois(data$Times_all,lambda=eval_lambda,pstr0=theta$prop_0,log=TRUE))},1)
  
  if (min(L)==-Inf | sum(is.na(L))>0){L=-10000}
  if (op=='sum'){L = sum(L)};
  return(L)
}


################### no spatial SE + spatially varying mu

lambda_cond_multidim_dt_nospatialSE_icar = function(theta,decay_fn,data,n_locations){
  
  h = theta$h
  mu_components = theta$mu
  mu = exp(mu_components$mu_intercept + mu_components$phi)
  
  M =1
  m=1
  y_ti  = lapply(1:M, function(l) {
    data$Times_events_index[[l]]$n
  })
  
  centroids = data$space_utilised$unique_space
  
  g=function(l){
    p = (m-1)*M+l
    decay_vals =  unlist(lapply(data$Times_obs, function(t){
      sapply(1:n_locations, function(k){
        t_ind = which(data$Times_events_index[[l]]$times %in% (t-smax[p]):(t-1) & 
                        data$Times_events_index[[l]]$centroids == centroids[k])
        t_ind_times = data$Times_events_index[[l]]$times[t_ind]
        sum(y_ti[[l]][t_ind] * decay_fn((t-t_ind_times),h$pars_time), na.rm=T)
       })
    }))
    L_m = h$alpha[p]*decay_vals
    return(L_m)
  }
  
  L_m=lapply(1:M,g)
  L_m=matrix(unlist(L_m), nrow=length(unlist(L_m))/M, ncol=M)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m= mu + R_m
  return(y_m)
  
}

# Log likelihood
log_likelihood_dt_nospatialSE_icar=function(data,decay_fn,theta_vec,n_locations,op='sum'){
  
  if (desc=="nospatialSE_icar_powerlaw"){
    theta = list(mu=list(mu_intercept=theta_vec[1],
                         phi=theta_vec[2:(2+n_locations-1)]),
                 h=list(alpha=theta_vec[(2+n_locations)],
                        pars_time=theta_vec[((2+n_locations)+1:2)]),
                 prop_0=theta_vec[((2+n_locations)+3)])
    
  } else if (desc=="nospatialSE_icar"){
    theta = list(mu=list(mu_intercept=theta_vec[1],
                         phi=theta_vec[2:(2+n_locations-1)]),
                 h=list(alpha=theta_vec[(2+n_locations)],
                        pars_time=theta_vec[((2+n_locations)+1)]),
                 prop_0=theta_vec[((2+n_locations)+2)])
    
  } 
  
  M = 1
  eval_lambda = lambda_cond_multidim_dt_nospatialSE_icar(theta,decay_fn,data,n_locations)
  L= vapply(1:M,function(m){sum(dzipois(data$Times_all,lambda=eval_lambda,pstr0=theta$prop_0,log=TRUE))},1)
  
  if (min(L)==-Inf | sum(is.na(L))>0){L=-10000}
  if (op=='sum'){L = sum(L)};
  return(L)
}

################### spatiotemporal SE + spatially varying mu


lambda_cond_multidim_dt_st_icar = function(theta,decay_fn,decay_space_fn,data,n_locations){
  
  h = theta$h
  mu_components = theta$mu
  mu = exp(mu_components$mu_intercept + mu_components$phi)
  
  M =1
  m=1

  
  decay_dists_vals  = lapply(1:M, function(l) {p=(m-1)*M+l;
  t(sapply(1:nrow(data$space_dists[[l]]), function(k) {
    decay_space_fn(data$space_dists[[l]][k,],h$pars_space) * data$Times_events_index[[l]]$n
  }))
  })

  g=function(l){
    p = (m-1)*M+l
    decay_space_vals =  unlist(lapply(data$Times_obs, function(t){
      t_ind = which(data$Times_events_index[[l]]$times %in% (t-data$smax[p]):(t-1))
      t_ind_times = data$Times_events_index[[l]]$times[t_ind]
      sapply(1:nrow(decay_dists_vals[[l]]), function(k) {
        sum(decay_dists_vals[[l]][k,t_ind] * decay_fn((t-t_ind_times),h$pars_time),na.rm=T)
      })
    }))
    L_m = h$alpha[p]*decay_space_vals
    return(L_m)
  }
  
  L_m=lapply(1:M,g)
  L_m=matrix(unlist(L_m), nrow=length(unlist(L_m))/M, ncol=M)
  if (is.vector(L_m)==TRUE){L_m=matrix(L_m,nrow=1,ncol=M)}
  R_m =rowSums(L_m)
  y_m= mu + R_m
  return(y_m)
  
}

# Log likelihood
log_likelihood_dt_st_icar=function(data,decay_fn,decay_space_fn,theta_vec,n_locations,op='sum'){
  
  if (desc=="ST-SE_icar_powerlaw"){
    theta = list(mu=list(mu_intercept=theta_vec[1],
                         phi=theta_vec[2:(2+n_locations-1)]),
                 h=list(alpha=theta_vec[(2+n_locations)],
                        pars_time=theta_vec[((2+n_locations)+1:2)],
                        pars_space=theta_vec[((2+n_locations)+3:4)]),
                 prop_0=theta_vec[((2+n_locations)+5)])
  
    
  } else if (desc=="ST-SE_icar"){
    theta = list(mu=list(mu_intercept=theta_vec[1],
                         phi=theta_vec[2:(2+n_locations-1)]),
                 h=list(alpha=theta_vec[(2+n_locations)],
                        pars_time=theta_vec[((2+n_locations)+1)],
                        pars_space=theta_vec[((2+n_locations)+2)]),
                 prop_0=theta_vec[((2+n_locations)+3)])
    
  } 
  
  M = 1
  eval_lambda = lambda_cond_multidim_dt_st_icar(theta,decay_fn,decay_space_fn,data,n_locations)
  L= vapply(1:M,function(m){sum(dzipois(data$Times_all,lambda=eval_lambda,pstr0=theta$prop_0,log=TRUE))},1)
  
  if (min(L)==-Inf | sum(is.na(L))>0){L=-10000}
  if (op=='sum'){L = sum(L)};
  return(L)
}
