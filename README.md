# Description

This is a repository to fit spatiotemporal discrete-time Hawkes processes to ACLED conflict data (https://acleddata.com). These anlayses make up Chapter 7 of my PhD thesis titled 'Bayesian Approaches for Modelling Discrete-Time Self-Exciting Processes and Their Applications' and the preprint is available on [arXiv](https://arxiv.org/abs/2408.14940).

# Data 

The file 'weekly_counts.rds' contains the processed data used for the analysis.

# Code

## MLEs

We use estimates of the MLEs to inform the priors in our Bayesian model. 

The script 'test_discretetime_st_mle_adm_only_zeroinf.R' does this for all conflict types and countries, in addition to considering several different model specifications: namely, no self-excitation, only temporal self-excitation, and temporal and spatial self-excitation.

The script 'run-mles.sh' runs 'test_discretetime_st_mle_adm_only_zeroinf.R' for a number of different maximum excitation times (12, 26 and 52 weeks).

## Bayesian model (using Stan)

The code is set up to perform each model run as a separate job on a high performance computing cluster. Each row in the file 'experiment_grid_ext.csv' contains the parameter settings for each unique job. The shell script 'run-stan.sh' creates an array job (should reflect each scenario in  'experiment_grid_ext.csv') and submits these jobs to the queue.

The file 'program_adm_only_icar_cmdstanr_zeroinf_mleprior.stan' contains Stan code for model with an exponential temporal triggering kernel and RBF spatial triggering kernel (temporal and spatial self-excitation), whereas 'program_adm_only_icar_cmdstanr_zeroinf_mleprior.stan' contains Stan code for a model with an exponential temporal triggering kernel and no spatial triggering kernel (only temporal self-excitation).

The R script 'test_discretetime_st_stan_adm_only_zeroinf.R' prepares the data and runs the Stan model. This is the main result from this analysis.
