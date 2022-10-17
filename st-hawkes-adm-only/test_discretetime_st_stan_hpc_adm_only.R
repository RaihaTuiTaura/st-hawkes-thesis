
.libPaths(Sys.getenv("R_LIB_USER"))
ind <- commandArgs(TRUE)
message("================================",ind)
ind <- as.numeric(ind)


library(parallel)
library(tidyverse)
library(magrittr)
library(readxl)
library(rstan)

sessionInfo()

source('functions/functions_likelihood_discretetime_st_stan.R')


#read in experiment grid:
grid <- read_csv("experiment_grid_ext.csv")
grid <- grid[ind,]

RNGkind("L'Ecuyer-CMRG")
options(mc.cores = parallel::detectCores())

#############################################################
################# Read in data ###################
#############################################################

# Load data (that includes population)
agg_data = readRDS("./conflict_ADM.Rdata")


#############################################################
################# Parameters ###################
#############################################################

country = grid$country
decay_fn = decay_geometric
decay_space_fn = decay_rbf
smax_p = grid$smax_p
time_aggregation = grid$time_aggregation
smax_p = grid$smax_p
end_year = 2015


# Filter data
sth_asia_st = agg_data$events %>%
  filter(YEAR < end_year & Country == country)

# Convert year and month to numeric value
start_year = min(sth_asia_st$YEAR)
start_month = min(sth_asia_st %>% filter(YEAR==start_year) %>% select(MONTH))
sth_asia_st %<>%
  mutate(times = case_when(
    (YEAR-start_year) == 0 ~ MONTH-start_month+1,
    (YEAR-start_year) > 0  ~ ((12-start_month+1) + (YEAR-start_year-1)*12 + MONTH)
  )) %>%
  select(times,EVENT_TYPE,ADM2_code,LATITUDE,LONGITUDE,EVENT_COUNTS)
colnames(sth_asia_st) = c("times","EVENT_TYPE","centroids","lat_centroid","long_centroid","n")

data_country_space = sth_asia_st %>%
  filter(n>0)


unique_types_individual = unique(data_country_space$EVENT_TYPE)
unique_types = list()
for (i in 1:length(unique_types_individual)){
  unique_types = c(unique_types,unique_types_individual[i])
}
unique_types[[(length(unique_types_individual)+1)]] = c("Protests","Riots")
type=unique_types[[grid$type_ind]]

print(type)
M = length(type)
smax = rep(smax_p,M^2)

data_country_space_type = data_country_space %>%
  filter(EVENT_TYPE %in% type)

unique_space =  unique(data_country_space_type$centroids)

if (grid$validation=="T"){
  remove_space = sample(unique_space,3,replace=F)
  data_country_space_type %<>% filter(!(centroids %in% remove_space))
  unique_space =  unique(data_country_space_type$centroids)
}

Times_events = lapply(1:M,function(m) {
  data_country_space_type %>%
    filter(EVENT_TYPE == type[m]) %>%
    select(times,centroids,n)
})

Tmax = max(data_country_space_type$times)
Tinf = max(smax)

Times_all = lapply(1:M,function(m) data.frame(times=rep(1:Tmax,rep(length(unique_space),Tmax)),centroids=rep(unique_space,Tmax)))
Times_all = lapply(1:M,function(m) left_join(Times_all[[m]], Times_events[[m]]))
Times_all = lapply(1:M,function(m) {Times_all[[m]][is.na(Times_all[[m]][,3]),3] = 0; return(Times_all[[m]][Times_all[[m]]$times>Tinf,3])})


space_all =  data.frame(unique(data_country_space_type %>% select(centroids,lat_centroid,long_centroid)))
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
data$Times_all = Times_all
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
t_ind_times = lapply(1:M,function(m) {
  Times_events_index[[m]]$times
})
data$t_ind_times_ragged = do.call(c,t_ind_times)
dist_size=attr(dists,"Size")
space_dists = lapply(1:M, function(m){
  t(sapply(1:nrow(space_utilised), function(k){
    y=as.vector(Times_events_index[[m]][,4])
    ind=convert_dist(k,y,dist_size)
    dists_ind =dists[ind]
    dists_ind[is.na(dists_ind)] = 0
    return(dists_ind)
  }))
})
data$space_dists_ragged = do.call(cbind,space_dists)


#-----------------------------
#--------- HMC
#------------------------------

fit <- stan(file = "program_adm_only.stan",
            data = data,
            iter = grid$N_MCMC,
            warmup = grid$burn_in,
            thin = 1,
            chains = grid$n_chains,
            verbose = F,
            seed = grid$seed_id,
            refresh = max(5000/10, 1),
            control = list(adapt_delta = 0.99)
)


message("================================",ind)
saveRDS(fit, paste0("acled_st_hawkes_adm_only_",time_aggregation,"_",ind,".rds"))
