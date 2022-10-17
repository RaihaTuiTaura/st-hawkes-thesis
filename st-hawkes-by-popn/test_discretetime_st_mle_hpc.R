
.libPaths(Sys.getenv("R_LIB_USER"))
ind <- commandArgs(TRUE)
message("================================",ind)
ind <- as.numeric(ind)

library(tidyverse)
library(readxl)
library(GGally)
library(magrittr)
library(reshape2)
library(parallel)
library(ape)

source('functions/functions_likelihood_discretetime_st_mle.R')


#read in experiment grid:
grid <- read_csv("experiment_grid_ext_mle.csv")
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
end_year = 2015
M=1



#############################################################
#################### MLEs for ST model #####################
#############################################################

# Filter data
sth_asia_st = agg_data$events %>%
  mutate(POPULATION_psk = POPULATION/AREA_km2) %>%
  filter(YEAR < end_year & Country == country)

# Add binary classifier (low/high) for population per square km
sth_asia_st %<>%
  mutate(pop_index_1 = POPULATION_psk <= quantile(sth_asia_st$POPULATION_psk,probs=0.75),
         pop_index_2 = POPULATION_psk > quantile(sth_asia_st$POPULATION_psk,probs=0.75))

head(sth_asia_st)

# Convert year and month to numeric value
start_year = min(sth_asia_st$YEAR)
start_month = min(sth_asia_st %>% filter(YEAR==start_year) %>% select(MONTH))
sth_asia_st %<>%
  mutate(times = case_when(
    (YEAR-start_year) == 0 ~ MONTH-start_month+1,
    (YEAR-start_year) > 0  ~ ((12-start_month+1) + (YEAR-start_year-1)*12 + MONTH)
  )) %>%
  select(times,EVENT_TYPE,ADM2_code,LATITUDE,LONGITUDE,pop_index_1,pop_index_2,EVENT_COUNTS)
colnames(sth_asia_st) = c("times","EVENT_TYPE","centroids","lat_centroid","long_centroid","pop_index_1","pop_index_2","n")

pop_index_data = sth_asia_st %>% select(times,EVENT_TYPE,centroids,pop_index_1,pop_index_2)

data_country_space = sth_asia_st %>%
  filter(n>0)

unique_types_individual = unique(data_country_space$EVENT_TYPE)
unique_types = list()
for (i in 1:length(unique_types_individual)){
  unique_types = c(unique_types,unique_types_individual[i])
}

print(unique_types)

res_all_types = lapply(1:length(unique_types), function(type_ind){
  type = unique_types[[type_ind]]
  print(type)
  M = length(type)
  smax = rep(smax_p,M^2)

  init_theta = c(rep(1,M),rep(1,M),rep(0.5,M^2),rep(0.5,M^2),rep(0.5,M^2),rep(1,M^2))

  data_country_space_type = data_country_space %>%
    filter(EVENT_TYPE %in% type)

  unique_space =  unique(data_country_space_type$centroids)

  Times_events = lapply(1:M,function(m) {
    data_country_space_type %>%
      filter(EVENT_TYPE == type[m]) %>%
      select(times,centroids,n)
  })

  Tmax = max(data_country_space_type$times)
  Tinf = max(smax)

  Times_all = lapply(1:M,function(m) data.frame(times=rep(1:Tmax,rep(length(unique_space),Tmax)),centroids=rep(unique_space,Tmax)))
  Times_all = lapply(1:M,function(m) left_join(Times_all[[m]], Times_events[[m]]))

  pop_index_1 = lapply(1:M,function(m) {pop=left_join(Times_all[[m]][Times_all[[m]]$times>Tinf,],(pop_index_data%>%filter(EVENT_TYPE %in% type)%>%select(times,centroids,pop_index_1)),
                                                  by=c("times"="times","centroids"="centroids")); return(as.numeric(pop$pop_index_1))})
  pop_index_2 = lapply(1:M,function(m) {pop=left_join(Times_all[[m]][Times_all[[m]]$times>Tinf,],(pop_index_data%>%filter(EVENT_TYPE %in% type)%>%select(times,centroids,pop_index_2)),
                                                  by=c("times"="times","centroids"="centroids")); return(as.numeric(pop$pop_index_2))})

  Times_all = lapply(1:M,function(m) {Times_all[[m]][is.na(Times_all[[m]][,3]),3] = 0; return(Times_all[[m]][Times_all[[m]]$times>Tinf,3])})


  space_all =  data.frame(unique(data_country_space_type %>% select(centroids,lat_centroid,long_centroid)))
  colnames(space_all) = c("unique_space","lat_centroid","long_centroid")
  dists = dist(space_all[,2:3])


  Times_obs=lapply(1:M,function(m){times=(Tinf+1):Tmax})

  data = list(Times_obs = Times_obs,
              Times_all = Times_all,
              pop_index_1 = pop_index_1,
              pop_index_2 = pop_index_2,
              Tinf = Tinf,
              Tmax = Tmax,
              smax = smax,
              dists = dists,
              M = M)

  data$space_utilised = space_all %>%
    mutate(space_index = 1:nrow(space_all)) %>%
    select(unique_space,space_index)
  data$Times_events_index = lapply(1:M,function(m) {
    left_join(Times_events[[m]], data$space_utilised, by=c("centroids"="unique_space"))
  })
  data$t_ind = lapply(1:M^2, function(p){ m =(p-1)%/%M+1; l = (p-1)%%M+1
    t(sapply(Times_obs[[m]], function(t) {
      t_ind = rep(0,nrow(data$Times_events_index[[l]]))
      t_ind[which(data$Times_events_index[[l]]$times %in% (t-smax[p]):(t-1))] = 1
      return(t_ind)
    }))
  })
  data$t_ind_times = lapply(1:M,function(m) {
    data$Times_events_index[[m]]$times
  })
  dist_size=attr(dists,"Size")
  data$space_dists = lapply(1:M, function(m){
    t(sapply(1:nrow(data$space_utilised), function(k){
      y=as.vector(data$Times_events_index[[m]][,4])
      ind=convert_dist(k,y,dist_size)
      dists_ind =dists[ind]
      dists_ind[is.na(dists_ind)] = 0
      return(dists_ind)
    }))
  })

  res_st = optim(init_theta, log_likelihood_multidim_dt_st, data = data, decay_fn=decay_fn, decay_space_fn=decay_space_fn,
                 control=list(fnscale=-1,pgtol=1e-5), lower=c(rep(0,M),rep(0,M),rep(0,M^2),rep(0,M^2),rep(0,M^2),rep(0,M^2)),
                 upper=c(rep(Inf,M),rep(Inf,M),rep(Inf,M^2),rep(Inf,M^2),rep(1,M^2),rep(Inf,M^2)), method="L-BFGS-B")
  mle_st = c(type,res_st$par)
  print(res_st$par)
  return(list(mle_st,data,res_st))
})


message("================================",ind)
saveRDS(res_all_types, paste0("acled_st_hawkes_mle_",time_aggregation,"_pop_",ind,".rds"))
