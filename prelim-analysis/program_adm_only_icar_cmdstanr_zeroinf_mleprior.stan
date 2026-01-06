// This code has been modified to include an ICAR prior on mu based on code from  
// Morris et al. (2019): Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan
// note the my function hasn't been adjusted to account for multiple dimensions

functions {


  // Geometric decay function
  vector geometric_decay(array[] int t, real beta_p, int len_t){

    array[len_t] real arg;
    vector[len_t] decay;

    for (i in 1:len_t){
      arg[i] = (1 - beta_p);
      decay[i] = beta_p * arg[i]^(t[i]-1);
    }
    return decay;
  }

  // Radial basis function kernel
  matrix rbf_kernel(matrix space_dists, real sigma, int n_events, int n_locations){

    matrix[n_locations,n_events] K;

    for (i in 1:n_events){
      for (j in 1:n_locations){
        K[j,i] = exp(-square(space_dists[j,i]) / (2*sigma^2));
      }
    }
    return K;
  }

  // Calculate conditional intensity function for each m
  vector lambda_cond_multidim_dt_m(matrix space_dists_ragged, matrix Times_events_index_ragged, array[] int Times_obs,
                                  matrix t_ind_ragged, array[] int t_ind_times_ragged, vector mu_vec, array[] real alpha, real beta,  array[] real sigma,
                                  int M, int n_locations, int T_obs, int n_rows_total, int m, array[] int n_events, int n_events_all){
    int p;
    int row_ind;
    matrix[n_rows_total,M] L_m;
    vector[n_rows_total] R_m;
    vector[n_rows_total] y_m;
    row_vector[M] ones;
    matrix[n_locations,n_events_all] decay_dists_vals;
    vector[n_locations] decay_space_vals;
    vector[n_events_all] decay_times_vals;
    array[n_events_all] int times;
    int pos_start = 1;
    int pos_end = 0;
    matrix[n_locations,n_events_all] rbf;
    int mu_ind;


    for (l in 1:M){p = (m-1) * M + l;
      pos_end = pos_end + n_events[l];
      rbf[,pos_start:pos_end] = rbf_kernel(space_dists_ragged[,pos_start:pos_end], sigma[p], n_events[l], n_locations);

      for (i in 1:n_locations){
        decay_dists_vals[i,pos_start:pos_end] = to_row_vector(col(Times_events_index_ragged[pos_start:pos_end],3)) .* rbf[i,pos_start:pos_end];
      }


      for (t in 1:T_obs){
        for (j in pos_start:pos_end){
          times[j] = Times_obs[t] - t_ind_times_ragged[j];
        }
        decay_times_vals[pos_start:pos_end] = geometric_decay(times[pos_start:pos_end],beta,n_events[l]) .* to_vector(t_ind_ragged[t,pos_start:pos_end]);
        for (i in 1:n_locations){
          decay_space_vals[i] = sum(to_vector(decay_dists_vals[i,pos_start:pos_end]) .* decay_times_vals[pos_start:pos_end] );
          row_ind = (t-1)*n_locations + i;
          L_m[row_ind,l] = alpha[p]*decay_space_vals[i];
        }
      }
      pos_start = pos_start + n_events[l];
    }

    for (i in 1:M){
      ones[i] = 1;
    }
    R_m = (ones * L_m')';
    for (row_ind_y in 1:n_rows_total){
      //print("Adding mu to intensity");
      mu_ind = row_ind_y % n_locations;
      if (mu_ind == 0)
        mu_ind = n_locations;
      //print(mu_ind);
      y_m[row_ind_y] = mu_vec[mu_ind] + R_m[row_ind_y];
      //print("Added mu to intensity");
    }

    return y_m;
  }
  
  real icar_normal_lpdf(vector phi, int n_locations, array[] int node1, array[] int node2) {
    // print("Printing nodes and phis");
    // print(node1);
    // print(node2);
    // print(phi[node1]);
    // print(phi[node2]);
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
 
  // // Count number of zeros in data 
  // int num_zeros(array[] int y) {
  //   int sum = 0;
  //   print(y)
  //   for (n in 1:size(y)) {
  //     sum += (y[n] == 0);
  //   }
  //   return sum;
  // }
  
  real zero_inf_poisson_lpmf(int y, real prop_0, real lambda){
    real loglik;

    if (y == 0) {
      loglik = log_sum_exp(log(prop_0),
                            log1m(prop_0)
                              + poisson_lpmf(y | lambda));
    } else {
      loglik = log1m(prop_0)
                  + poisson_lpmf(y | lambda);
    }
    
    return loglik;
  }
 
}

// The input data.
data {
  int<lower = 0> M;
  int<lower = 0> Msq;
  int<lower = 0> T_obs;
  int<lower = 0> n_locations;
  array[M] int<lower=0> n_events;
  int<lower = 0> n_events_all;
  int<lower = 0> n_events_all_M;
  int<lower = 0> n_rows_total;
  array[T_obs] int<lower=0> Times_obs;
  array[M, n_rows_total] int<lower=0> Times_all;
  matrix[n_events_all, 4] Times_events_index_ragged;
  matrix<lower = 0>[T_obs,n_events_all_M] t_ind_ragged;
  array[n_events_all] int<lower=0> t_ind_times_ragged;
  matrix<lower = 0>[n_locations,n_events_all] space_dists_ragged;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=n_locations> node1;  // node1[i], node2[i] neighbors
  array[N_edges] int<lower=1, upper=n_locations> node2;  // node1[i] < node2[i]
  real<lower = 0> mu_sd;
  real mu_mean;
  real<lower = 0> alpha_sd;
  real alpha_mean;
  real<lower = 0> sigma_sd;
  real sigma_mean;
  real<lower = 0> beta_a;
  real<lower = 0> beta_b;
  real<lower = 0> prop_0_a;
  real<lower = 0> prop_0_b;
}

// transformed data {
//   int<lower=0> N_zero = num_zeros(Times_all[1]);
//   array[n_rows_total - N_zero] int<lower=1> y_nonzero;
//   int N_nonzero = 0;
//   for (n in 1:n_rows_total) {
//     if (Times_all[n:,1] == 0) continue;
//     N_nonzero += 1;
//     y_nonzero[N_nonzero] = Times_all[n:,1];
//   print(N_zero)
//   print(y_nonzero)
//   print(N_nonzero)
//   }
// }

// The parameters accepted by the model.
parameters {
  real mu_intercept;                // overall baseline risk
  sum_to_zero_vector[n_locations] phi;         // spatial effects
  array[Msq] real alpha_log;
  real beta_logit;
  array[Msq] real sigma_log;
  real prop_0_logit;
}

transformed parameters {
  array[Msq] real<lower=0> alpha = exp(alpha_log);
  real<lower=0, upper=1> beta = inv_logit(beta_logit);
  array[Msq] real<lower=0> sigma = exp(sigma_log);
  real<lower=0, upper=1> prop_0 = inv_logit(prop_0_logit);
  jacobian += alpha_log ;
  jacobian += (-beta_logit - 2*log(1+exp(-beta_logit))) ; 
  jacobian += sigma_log ;
  jacobian += (-prop_0_logit - 2*log(1+exp(-prop_0_logit))) ; //for alpha, beta, sigma and prop_0
}

// The model to be estimated.
model {

  array[M] vector[n_rows_total] L;
  int row_ind;
  int pos;
  vector[n_locations] mu_vec;
  // prior for mus
  //print("Calculating mu_vec");
  mu_vec = exp(mu_intercept + phi);
  //print("Calculated mu_vec");
  //print(mu_vec);
  mu_intercept ~ normal(mu_mean, mu_sd);
  phi ~ icar_normal(n_locations, node1, node2);


  // priors for alphas, betas and sigmas
  for (p in 1:Msq){
    alpha[p] ~ lognormal(alpha_mean,alpha_sd);
    sigma[p] ~ lognormal(sigma_mean,sigma_sd);
  }
  beta ~ beta(beta_a,beta_b);
  prop_0 ~ normal(prop_0_a,prop_0_b);

  pos = 1;
  for (m in 1:M){

    matrix[T_obs,n_events_all] t_ind_m;
    t_ind_m = block(t_ind_ragged, 1, pos, T_obs, n_events_all);


    // Calculate rate for the Poisson distribution
    L[m] = lambda_cond_multidim_dt_m(space_dists_ragged, Times_events_index_ragged, Times_obs, t_ind_m, t_ind_times_ragged,
                                  mu_vec, alpha, beta, sigma, M, n_locations, T_obs, n_rows_total, m, n_events, n_events_all);
    pos = pos + n_events_all;
  }

  // y^M_t ~ Poisson(L_m(t))
  for (m in 1:M){
    for (t in 1:T_obs){
      for (i in 1:n_locations){
        row_ind = (t-1)*n_locations + i;
        Times_all[m,row_ind] ~ zero_inf_poisson(prop_0,L[m,row_ind]);
      }
    }
  }
}
