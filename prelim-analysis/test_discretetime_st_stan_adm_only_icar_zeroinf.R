# working directory should be the parent folder of the github repo 'st-hawkes-thesis'
# this script requires geopackages files obtained from https://gadm.org/download_country.html
# also requires output from MLE estimates

ind <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("================================",ind)
ind <- as.numeric(ind)

library(parallel)
library(dplyr)
library(readr)
library(magrittr)
library(readxl)
library(sf)
library(cmdstanr)
library(INLA)
library(spdep)

source('prelim-analysis/functions/functions_likelihood_discretetime_st_stan_zeroinf.R')
RNGkind("L'Ecuyer-CMRG")
options(mc.cores = parallel::detectCores())

#############################################################
################# Read in data ###################
#############################################################

# Load data
agg_data = readRDS("weekly_counts.rds")


#############################################################
################# Parameters ###################
#############################################################

decay_fn = decay_geometric
decay_space_fn = decay_rbf

#read in experiment grid:
grid <- read_csv("prelim-analysis/experiment_grid_ext.csv")
grid <- grid[ind,]

desc = grid$desc
seed_id = grid$seed_id
country = grid$country
ISO3C = grid$ISO3C
iter_sampling=grid$iter_sampling
burn_in=grid$burn_in
smax_p = grid$smax_p
n_chains= grid$n_chains
type_ind = grid$type_ind
time_aggregation = grid$time_aggregation
program_file_name = grid$program_file_name

# Filter data
sth_asia_st = agg_data$data %>%
  filter(Country == country) 

# Keep relevant columns
sth_asia_st %<>%
  select(WEEK_COUNT,EVENT_TYPE,ADM_code,lat,long,EVENT_COUNTS)
colnames(sth_asia_st) = c("times","EVENT_TYPE","centroids","lat_centroid","long_centroid","n")

data_country_space = sth_asia_st %>%
  filter(n>0)


unique_types_individual = sort(unique(data_country_space$EVENT_TYPE))
unique_types = list()
for (i in 1:length(unique_types_individual)){
  unique_types = c(unique_types,unique_types_individual[i])
}
type=unique_types[[type_ind]]

print(type)
M = length(type)
smax = rep(smax_p,M^2)

data_country_space_type = data_country_space %>%
  filter(EVENT_TYPE %in% type)

unique_space =  unique(sth_asia_st$centroids)

Times_events = lapply(1:M,function(m) {
  data_country_space_type %>%
    filter(EVENT_TYPE == type[m]) %>%
    select(times,centroids,n)
})

Tmax = max(agg_data$data$WEEK_COUNT)
Tinf = max(smax)

Times_all = lapply(1:M,function(m) data.frame(times=rep(1:Tmax,rep(length(unique_space),Tmax)),centroids=rep(unique_space,Tmax)))
Times_all = lapply(1:M,function(m) left_join(Times_all[[m]], Times_events[[m]]))
Times_all = lapply(1:M,function(m) {Times_all[[m]][is.na(Times_all[[m]][,3]),3] = 0; return(Times_all[[m]][Times_all[[m]]$times>Tinf,3])})


space_all =  data.frame(unique(sth_asia_st %>% select(centroids,lat_centroid,long_centroid)))
colnames(space_all) = c("unique_space","lat_centroid","long_centroid")
dists = dist(space_all[,2:3])

Times_obs=times=(Tinf+1):Tmax

data = list(M=M)
data$Msq = M^2
data$T_obs = length(Times_obs)
data$n_locations = length(unique_space)
if (M==1){
  data$n_events = array(sapply(1:M,function(m) {nrow(Times_events[[m]])}),dim=1)
} else {
  data$n_events = sapply(1:M,function(m) {nrow(Times_events[[m]])})
}
data$n_events_all = sum(data$n_events)
data$n_events_all_M = data$n_events_all*M
data$n_rows_total = data$n_locations * data$T_obs
data$Times_obs = Times_obs
data$Times_all = array(unlist(Times_all), c(M,data$n_rows_total))
space_utilised = space_all %>%
  mutate(space_index = 1:nrow(space_all)) %>%
  select(unique_space,space_index)
Times_events_index = lapply(1:M,function(m) {
  left_join(Times_events[[m]], space_utilised, by=c("centroids"="unique_space"))
})
data$Times_events_index_ragged = do.call(rbind,Times_events_index)
t_ind = lapply(1:M^2, function(p){ m =(p-1)%/%M+1; l = (p-1)%%M+1
t(sapply(Times_obs, function(t) {
  t_ind = rep(0,data$n_events[l])
  t_ind[which(Times_events_index[[l]]$times %in% (t-smax[p]):(t-1))] = 1
  return(t_ind)
}))
})
data$t_ind_ragged = do.call(cbind,t_ind)

space_ind = lapply(1:M^2, function(p){ m =(p-1)%/%M+1; l = (p-1)%%M+1
t(sapply(1:length(unique_space), function(k) {
  space_ind = rep(0,data$n_events[l])
  space_ind[which(Times_events_index[[l]]$space_index == k)] = 1
  return(space_ind)
}))
})
data$space_ind_ragged = do.call(cbind,space_ind)

t_ind_times = lapply(1:M,function(m) {
  Times_events_index[[m]]$times
})
data$t_ind_times_ragged = do.call(c,t_ind_times)
dist_size=attr(dists,"Size")
space_dists = lapply(1:M, function(m){
  t(sapply(1:nrow(space_utilised), function(k){
    y=Times_events_index[[m]][,"space_index",drop=F]
    ind=convert_dist(k,y,dist_size)
    dists_ind =dists[ind]
    dists_ind[is.na(dists_ind)] = 0
    return(dists_ind)
  }))
})
data$space_dists_ragged = do.call(cbind,space_dists)

# make neighbourhood graph and get neighbours 
#c("BGD","LKA","NPL","PAK")
if(ISO3C %in% c("LKA","BGD")){
  boundaries = st_read(paste0("geopackages/gadm41_",ISO3C,".gpkg"),layer="ADM_ADM_1") %>%
    filter(GID_1 %in% sth_asia_st$centroids) %>%
    select(GID_1,geom)
  
} else {
  boundaries = st_read(paste0("geopackages/gadm41_",ISO3C,".gpkg"),layer="ADM_ADM_2") %>%
    filter(GID_2 %in% sth_asia_st$centroids) %>%
    select(GID_2,geom)
}
colnames(boundaries) = c("ADM_code","geom")
boundaries = boundaries[match(space_all$unique_space,boundaries$ADM_code),]
if (identical(boundaries$ADM_code,space_all$unique_space) == F){
  print("Issue with order of neighbourhoods")
}

source("prelim-analysis/functions/nb_data_funs.R")
nb=poly2nb(boundaries)
## bym2 model requires symmetric neighborhood structure
validateNb(nb);
## check that nb object is fully connected.
!(isDisconnected(nb))
# entire graph represented as node pairs
nbs=nb2graph(nb);
data$N_edges = nbs$N_edges;
data$node1 = nbs$node1;
data$node2 = nbs$node2;

#-----------------------------
#--------- HMC
#------------------------------

print(country)
print(type)

# other file names
res_filename = paste0("./prelim-analysis/outputs-stan-zeroinf/","res_",desc,"_",time_aggregation,"_adm_only_",country,"_",gsub('/', '',type),"_smax",smax,"_zeroinf_transformed_",ind,".rds")
data_filename = paste0("./prelim-analysis/outputs-stan-zeroinf/","data_",desc,"_",time_aggregation,"_adm_only_",country,"_",gsub('/', '',type),"_smax",smax,"_zeroinf_transformed_",ind,".rds")
output_csv_dir ="./prelim-analysis/outputs-stan-zeroinf"
output_csv_basename = paste0("chains_",desc,"_",time_aggregation,"_adm_only_",country,"_",gsub('/', '',type),"_smax",smax,"_zeroinf_transformed_",ind,"_")
if (grepl("ST_SE_icar",desc)){
  mle_results = paste0("./prelim-analysis/outputs-zeroinf/ST-SE_icar_mle_smax",smax,"_",time_aggregation,"_adm_only_",country,"_zeroinf.rds")
} else if (grepl("nospatialSE_icar",desc)){
  mle_results = paste0("./prelim-analysis/outputs-zeroinf/nospatialSE_icar_mle_smax",smax,"_",time_aggregation,"_adm_only_",country,"_zeroinf.rds")
}

file <- file.path(program_file_name)
mod <- cmdstan_model(file)

# load mles (needed for desc="*_mleprior" and desc="ST_SE_icar")
mles = readRDS(mle_results)
mles_type = mles[[which(unique_types==type)]][[3]]$par
if (grepl("ST_SE_icar",desc)){
  names(mles_type) = c("mu_intercept",paste0("phi",1:data$n_locations),"alpha","beta","sigma","prop_0")
} else if (grepl("nospatialSE_icar",desc)){
  names(mles_type) = c("mu_intercept",paste0("phi",1:data$n_locations),"alpha","beta","prop_0")
}


# add small amount onto 0 values of prop_0
if(mles_type["prop_0"]<1e-6){
  mles_type["prop_0"] = 1e-6
}
# add small amount onto 0 values of alpha
if(mles_type["alpha"]<1e-6){
  mles_type["alpha"] = 1e-6
}
# add small amount onto 0 values of beta
if(mles_type["beta"]<1e-6){
  mles_type["beta"] = 1e-6
}
if (grepl("ST_SE_icar",desc)){
  
  # add small amount onto 0 values of sigma
  if(mles_type["sigma"]<1e-6){
    mles_type["sigma"] = 1e-6
  }
} 

# subtract small amount from large beta
if(mles_type["beta"]>0.99999){
  mles_type["beta"] = 0.99999
}
# subtract small amount from large prop_0
if(mles_type["prop_0"]>0.99999){
  mles_type["prop_0"] = 0.99999
}

if (grepl("mleprior",desc)){
  
  data$mu_sd = 0.1
  data$mu_mean = as.numeric(mles_type['mu_intercept'])
  data$alpha_sd = 0.1
  data$alpha_mean = log(as.numeric(mles_type['alpha'])) - data$alpha_sd^2/2
  
  if (grepl("ST_SE_icar",desc)){
    data$sigma_sd = 0.1
    data$sigma_mean = log(as.numeric(mles_type['sigma'])) - data$sigma_sd^2/2
  }
  
  data$beta_a = 100
  if (as.numeric(mles_type['beta']) == 0){
    data$beta_b = data$beta_a * (1-0.01)/0.01
  } else {
    data$beta_b = data$beta_a * (1-as.numeric(mles_type['beta']))/as.numeric(mles_type['beta'])
    
  }
  data$prop_0_a = 100
  if (as.numeric(mles_type['prop_0']) == 0){
    data$prop_0_b = data$prop_0_a * (1-0.01)/0.01
  } else {
    data$prop_0_b = data$prop_0_a * (1-as.numeric(mles_type['prop_0']))/as.numeric(mles_type['prop_0'])
    
  }
}

# transform (do after above since we are putting prior still on the constrained parameters)
# take log of alpha and sigma
mles_type["alpha"] = log(mles_type["alpha"])

# take logit of beta and prop_0
mles_type["beta"] = qlogis(mles_type["beta"])
mles_type["prop_0"] = qlogis(mles_type["prop_0"])
if (grepl("ST_SE_icar",desc)){
  mles_type["sigma"] = log(mles_type["sigma"])
  names(mles_type) = c("mu_intercept",paste0("phi",1:data$n_locations),"alpha_log","beta_logit","sigma_log","prop_0_logit")
  
}else if (grepl("nospatialSE_icar",desc)){
  names(mles_type) = c("mu_intercept",paste0("phi",1:data$n_locations),"alpha_log","beta_logit","prop_0_logit")
  
}



fit <- mod$sample(data = data,
                  iter_warmup = burn_in,
                  iter_sampling = iter_sampling,
                  chains = n_chains,
                  seed = seed_id,
                  show_messages = TRUE)


#save
fit$save_object(file=res_filename)
fit$save_output_files(dir = output_csv_dir, basename = output_csv_basename, timestamp = FALSE, random = FALSE)
saveRDS(data, data_filename)

sessionInfo()

message("================================",ind)
