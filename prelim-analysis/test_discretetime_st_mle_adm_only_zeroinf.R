
library(tidyverse)
library(readxl)
library(GGally)
library(magrittr)
library(reshape2)
library(parallel)
library(ape)
library(sf)
library(INLA)
library(VGAM)

source('prelim-analysis/functions/functions_likelihood_discretetime_st_mle_zeroinf.R')
source("prelim-analysis/functions/nb_data_funs.R")


RNGkind("L'Ecuyer-CMRG")
options(mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
smax = smax_p = as.numeric(args[1])
print(paste0("smax=",smax))
print(paste0("is.numeric=",is.numeric(smax)))

#############################################################
################# Read in data ###################
#############################################################

# Load data
agg_data = readRDS("weekly_counts.rds")


#############################################################
################# Parameters ###################
#############################################################


time_aggregation = "weeks"
M=1


countries = c("Bangladesh", "Sri Lanka", "Nepal", "Pakistan")
iso3codes = c("BGD","LKA","NPL","PAK")



#############################################################
#################### MLEs for ST model #####################
#############################################################

run_mle = function(country){
  
  set.seed(NULL)
  set.seed(3956+which(countries==country))
  
  ISO3C = iso3codes[which(countries==country)]
  
  # Filter data
  sth_asia_st = agg_data$data %>%
    filter(Country == country)
  
  # Keep relevant columns
  sth_asia_st %<>%
    dplyr::select(WEEK_COUNT,EVENT_TYPE,ADM_code,lat,long,EVENT_COUNTS)
  colnames(sth_asia_st) = c("times","EVENT_TYPE","centroids","lat_centroid","long_centroid","n")
  
  data_country_space = sth_asia_st %>%
    filter(n>0)
  
  unique_types_individual = sort(unique(data_country_space$EVENT_TYPE))
  unique_types = list()
  for (i in 1:length(unique_types_individual)){
    unique_types = c(unique_types,unique_types_individual[i])
  }
  
  print(unique_types)
  
  space_all =  data.frame(unique(sth_asia_st %>% dplyr::select(centroids,lat_centroid,long_centroid)))
  colnames(space_all) = c("unique_space","lat_centroid","long_centroid")
  
  unique_space =  unique(sth_asia_st$centroids)
  
  Tmax = max(sth_asia_st$times)
  Tinf = max(smax)
  
  # create object to store data
  data = list(
    Tinf = Tinf,
    Tmax = Tmax,
    smax = smax,
    M = M)
  
  # make neighbourhood graph and get neighbours 
  #c("BGD","LKA","NPL","PAK")
  if(ISO3C %in% c("BGD","LKA")){
    boundaries = st_read(paste0("geopackages/gadm41_",ISO3C,".gpkg"),layer="ADM_ADM_1") %>%
      filter(GID_1 %in% sth_asia_st$centroids) %>%
      dplyr::select(GID_1,geom)
    ###### THIS WAS FILTERING ON EVENT TYPE. CHECK OTHER CODE BASED ON THIS 
  } else {
    boundaries = st_read(paste0("geopackages/gadm41_",ISO3C,".gpkg"),layer="ADM_ADM_2") %>%
      filter(GID_2 %in% sth_asia_st$centroids) %>%
      dplyr::select(GID_2,geom)
  }
  colnames(boundaries) = c("ADM_code","geom")
  boundaries = boundaries[match(space_all$unique_space,boundaries$ADM_code),]
  if (identical(boundaries$ADM_code,space_all$unique_space) == F){
    print("Issue with order of neighbourhoods")
  }
  
  nb=poly2nb(boundaries)
  ## bym2 model requires symmetric neighborhood structure
  validateNb(nb)
  ## check that nb object is fully connected.
  !(isDisconnected(nb))
  # calculate scaling factor
  data$scaling_factor = scale_nb_components(nb)[1];
  # entire graph represented as node pairs
  nbs=nb2graph(nb);
  data$N_edges = nbs$N_edges
  data$node1 = nbs$node1
  data$node2 = nbs$node2
  n_locations = nrow(boundaries)
  
  data$space_utilised = space_all %>%
    mutate(space_index = 1:nrow(space_all)) %>%
    dplyr::select(unique_space,space_index)
  
  dists = dist(space_all[,2:3])
  dist_size=attr(dists,"Size")
  
  
  res_all_types = lapply(1:length(unique_types), function(type_ind){
    
    type = unique_types[[type_ind]]
    print(type)
    M = length(type)
    smax = rep(smax_p,M^2)
    
    data_country_space_type = data_country_space %>%
      filter(EVENT_TYPE %in% type)
    
    Times_events = lapply(1:M,function(m) {
      data_country_space_type %>%
        filter(EVENT_TYPE == type[m]) %>%
        dplyr::select(times,centroids,n) %>%
        arrange(times,centroids)
    })

    
    Times_all = lapply(1:M,function(m) data.frame(times=rep(1:Tmax,rep(length(unique_space),Tmax)),centroids=rep(unique_space,Tmax)))
    Times_all = lapply(1:M,function(m) left_join(Times_all[[m]], Times_events[[m]]))
    Times_all = lapply(1:M,function(m) {Times_all[[m]][is.na(Times_all[[m]][,3]),3] = 0; return(Times_all[[m]][Times_all[[m]]$times>Tinf,3])})
    
    data$Times_all = Times_all[[1]]
    data$Times_obs=times=(Tinf+1):Tmax
    
    if (!desc=="nospatialSE_constBL"){
      data$Times_events_index = lapply(1:M,function(m) {
        Times_events[[m]] %>% 
          left_join(data$space_utilised, by=c("centroids"="unique_space"))
      })
      
      data$space_dists = lapply(1:M, function(m){
        t(sapply(1:nrow(data$space_utilised), function(k){
          y=data$Times_events_index[[m]][,4,drop=F]
          ind=convert_dist(k,y,dist_size)
          dists_ind =dists[ind]
          dists_ind[is.na(dists_ind)] = 0
          return(dists_ind)
        }))
      })
    } else {
      data$Times_events_index = lapply(1:M,function(m) {
        Times_events[[m]] %>% group_by(times) %>% summarise(n = sum(n)) %>% ungroup() 
      })
    }
  
    if (desc=="ST-SE_icar_powerlaw"){
      
      init_theta =  c(0.5,rep(0,n_locations),0.5,1.5,0.5,1.5,0.5,0.5)
      lhood_fun = log_likelihood_dt_st_icar
      lower = c(-Inf,rep(-Inf,n_locations),0,1,0,1,0,0)
      upper = c(Inf,rep(Inf,n_locations),Inf,Inf,Inf,Inf,Inf,1)
      
    }  else if (desc=="ST-SE_icar"){
      
      init_theta =  c(0.5,rep(0,n_locations),0.5,0.5,0.5,0.5)
      lhood_fun = log_likelihood_dt_st_icar
      lower = c(-Inf,rep(-Inf,n_locations),0,0,0,0)
      upper = c(Inf,rep(Inf,n_locations),Inf,1,Inf,1)
      
    } else if (desc=="nospatialSE_icar_powerlaw"){
      init_theta =  c(0.5,rep(0,n_locations),0.5,1.5,0.5,0.5)
      lhood_fun = log_likelihood_dt_nospatialSE_icar
      lower = c(-Inf,rep(-Inf,n_locations),0,1,0,0)
      upper = c(Inf,rep(Inf,n_locations),Inf,Inf,Inf,1)
      
    } else if (desc=="nospatialSE_icar"){
      
      init_theta = c(0.5,rep(0,n_locations),0.5,0.5,0.5)
      lhood_fun = log_likelihood_dt_nospatialSE_icar
      lower = c(-Inf,rep(-Inf,n_locations),0,0,0)
      upper = c(Inf,rep(Inf,n_locations),Inf,1,1)
      
    }  else if (desc=="icar_only"){
      init_theta =  c(0.5,rep(0,n_locations),0.5)
      lhood_fun = log_likelihood_dt_icar_only
      lower = c(-Inf,rep(-Inf,n_locations),0)
      upper = c(Inf,rep(Inf,n_locations),1)
      
    } 
    
    if (desc %in% c("ST-SE_icar","ST-SE_icar_powerlaw","ST-SE_constBL")){
      res_st = optim(init_theta, lhood_fun, data = data, decay_fn=decay_fn, decay_space_fn=decay_space_fn,n_locations=n_locations,
                     control=list(fnscale=-1,pgtol=1e-5), lower=lower,upper=upper, method="L-BFGS-B")
      
    } else if (desc %in% c("nospatialSE_icar","nospatialSE_icar_powerlaw","nospatialSE_constBL","icar_only")) {
      res_st = optim(init_theta, lhood_fun, data = data, decay_fn=decay_fn, n_locations=n_locations,
                     control=list(fnscale=-1,pgtol=1e-5), lower=lower,upper=upper, method="L-BFGS-B")
      
    }
    
  mle_st = c(type,res_st$par)
  print(res_st$par)
  return(list(ests=mle_st,data=data,optim_res=res_st))
  })
  
  ind=which(countries==country)
  if (desc == "icar_only"){
    saveRDS(res_all_types, paste0("prelim-analysis/outputs-zeroinf/",desc,"_mle_",time_aggregation,"_adm_only_",country,"_zeroinf.rds"))
    
  } else {
    saveRDS(res_all_types, paste0("prelim-analysis/outputs-zeroinf/",desc,"_mle_smax",smax,"_",time_aggregation,"_adm_only_",country,"_zeroinf.rds"))
    
  }
  return(c(paste0("saved output for ",desc)))
}


decay_fn = decay_geometric
decay_space_fn = decay_rbf
start.time = proc.time()
desc="nospatialSE_icar"
mclapply(countries, function(country) run_mle(country))
print(desc)
print(proc.time()-start.time)

decay_fn = decay_geometric
decay_space_fn = decay_rbf
start.time = proc.time()
desc="ST-SE_icar"
mclapply(countries, function(country) run_mle(country))
print(desc)
print(proc.time()-start.time)

decay_fn = decay_powerlaw
decay_space_fn = NULL
start.time = proc.time()
desc="nospatialSE_icar_powerlaw"
mclapply(countries, function(country) run_mle(country))
print(desc)
print(proc.time()-start.time)

decay_fn = decay_powerlaw
decay_space_fn = decay_powerlaw
start.time = proc.time()
desc="ST-SE_icar_powerlaw"
mclapply(countries, function(country) run_mle(country))
print(desc)
print(proc.time()-start.time)

# # baseline only model
# # note: just ran this separately to above
# decay_fn = NULL
# decay_space_fn = NULL
# smax=smax_p=0
# start.time = proc.time()
# desc="icar_only"
# mclapply(countries, function(country) run_mle(country))
# print(desc)
# print(proc.time()-start.time)
