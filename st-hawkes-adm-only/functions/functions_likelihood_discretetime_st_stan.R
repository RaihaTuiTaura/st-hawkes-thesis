
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
