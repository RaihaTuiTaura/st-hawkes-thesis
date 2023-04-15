# Description

This is a repository to fit spatiotemporal discrete-time Hawkes processes to ACLED conflict data (https://acleddata.com). These anlayses make up Chapter 7 of my PhD thesis titled 'Bayesian Approaches for Modelling Discrete-Time Self-Exciting Processes and Their Applications'. Preprint coming soon.

There are two sets of analyses: with and without population included in the model. Both analyses are temporally discretised at the monthly level, and spatially discretised at level 2 spatial administration boundaries (ADM).

The following details describe the folder 'st-hawkes-adm-only' which contains the analysis exluding population information. However, a similar structure is followed in the folder 'st-hawkes-by-popn' which includes a coarse indicator of population in the model.

# Data 

The file 'conflict_ADM.Rdata' contains the processed data used for the analysis.

# Code

The code is set up to perform each model run as a separate job on a high performance computing cluster. Each row in the file 'experiment_grid_ext.csv' contains the parameter settings for each unique job. The shell script 'run-file-adm-only.sh' creates a job using the template file 'generic-job-adm-only.sh' and submits the job to the queue.

The file 'program_adm_only.stan' contains Stan code for model with an exponential temporal triggering kernel and RBF spatial triggering kernel. The R script 'test_discretetime_st_stan_adm_only.R' prepares the data and runs the Stan model. This is the main result from this anlaysis.

There are several other examples produced, based on the MLEs. The R script 'test_discretetime_st_mle_adm_only.R' prepares the data and estimates MLEs for the model. 'test_discretetime_st_temporal_monthly.R' prepares the data and estimates the MLEs of a temporal Hawkes process with no spatial component. Lastly, 'test_discretetime_st_temp_comparison.R' produces a comparison model to compare temporal and spatiotemporal Hawkes processes (essentially fits a temporal model at the ADM level for each location).
