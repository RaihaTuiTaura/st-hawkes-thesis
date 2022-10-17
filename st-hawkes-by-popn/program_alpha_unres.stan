//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {

  // Geometric decay function
  vector geometric_decay(int[] t, real beta_p, int len_t){

    real arg[len_t];
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
        K[j,i] =  exp(-square(space_dists[j,i]) / (2*sigma^2));
      }
    }
    return K;
  }

  // Calculate conditional intensity function for each m
  vector lambda_cond_multidim_dt_m(matrix space_dists_ragged, matrix Times_events_index_ragged, int[] Times_obs,
                                  matrix t_ind_ragged, int[] t_ind_times_ragged, vector[] pop_index_1, vector[] pop_index_2,
                                  real[] mu_1, real[] mu_2, real[] alpha_1, real[] alpha_2, real[] beta,  real[] sigma,
                                  int M, int n_locations, int T_obs, int n_rows_total, int m, int[] n_events, int n_events_all){
    int p;
    int row_ind;
    matrix[n_rows_total,M] L_m;
    vector[n_rows_total] R_m;
    vector[n_rows_total] y_m;
    row_vector[M] ones;
    matrix[n_locations,n_events_all] decay_dists_vals;
    vector[n_locations] decay_space_vals;
    vector[n_events_all] decay_times_vals;
    int times[n_events_all];
    int pos_start = 1;
    int pos_end = 0;
    matrix[n_locations,n_events_all] rbf;


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
        decay_times_vals[pos_start:pos_end] = geometric_decay(times[pos_start:pos_end],beta[p],n_events[l]) .* to_vector(t_ind_ragged[t,pos_start:pos_end]);
        for (i in 1:n_locations){
          decay_space_vals[i] = sum(to_vector(decay_dists_vals[i,pos_start:pos_end]) .* decay_times_vals[pos_start:pos_end] );
          row_ind = (t-1)*n_locations + i;
          //L_m[row_ind,l] = (alpha_1[p]*(pop_index_1[m,row_ind]+pop_index_2[m,row_ind]))*decay_space_vals[i];
          L_m[row_ind,l] = (alpha_1[p]*pop_index_1[m,row_ind] + alpha_2[p]*pop_index_2[m,row_ind])*decay_space_vals[i];
        }
      }
      pos_start = pos_start + n_events[l];
    }

    for (i in 1:M){
      ones[i] = 1;
    }
    R_m = (ones * L_m')';
    for (row_ind_y in 1:n_rows_total){
      //y_m[row_ind_y] = (mu_1[m]*(pop_index_1[m,row_ind_y]+pop_index_2[m,row_ind_y])) + R_m[row_ind_y];
      y_m[row_ind_y] = (mu_1[m]*pop_index_1[m,row_ind_y] + mu_2[m]*pop_index_2[m,row_ind_y]) + R_m[row_ind_y];
    }

    return y_m;
  }
}

// The input data.
data {
  int<lower = 0> M;
  int<lower = 0> Msq;
  int<lower = 0> T_obs;
  int<lower = 0> n_locations;
  int<lower = 0> n_events[M];
  int<lower = 0> n_events_all;
  int<lower = 0> n_events_all_M;
  int<lower = 0> n_rows_total;
  int<lower = 0> Times_obs[T_obs];
  int<lower = 0> Times_all[M,n_rows_total];
  matrix<lower = 0>[T_obs,n_events_all_M] t_ind_ragged;
  int<lower = 0> t_ind_times_ragged[n_events_all];
  matrix[n_events_all, 4] Times_events_index_ragged;
  matrix<lower = 0>[n_locations,n_events_all] space_dists_ragged;
  vector[n_rows_total] pop_index_1[M];
  vector[n_rows_total] pop_index_2[M];
}


// The parameters accepted by the model.
parameters {
  real<lower = 0> mu_1[M];
  real<lower = 0> mu_2[M];
  real<lower = 0> alpha_1[Msq];
  real<lower = 0> alpha_2[Msq];
  real<lower = 0, upper=1> beta[Msq];
  real<lower = 0> sigma[Msq];
}

// The model to be estimated.
model {

  vector[T_obs] L[M];
  int row_ind;
  int pos;


  // prior for mus
  for (m in 1:M){
    mu_1[m] ~ gamma(2,2);
    mu_2[m] ~ gamma(2,2);
  }

  // priors for alphas, betas and sigmas
  for (p in 1:Msq){
    alpha_1[p] ~ gamma(2,2);
    alpha_2[p] ~ gamma(2,2);
    beta[p] ~ uniform(0,1);
    sigma[p] ~ inv_gamma(5, 5);
  }

  pos = 1;
  for (m in 1:M){

    matrix[T_obs,n_events_all] t_ind_m;
    t_ind_m = block(t_ind_ragged, 1, pos, T_obs, n_events_all);


    // Calculate rate for the Poisson distribution
    L[m] = lambda_cond_multidim_dt_m(space_dists_ragged, Times_events_index_ragged, Times_obs, t_ind_m, t_ind_times_ragged,pop_index_1,pop_index_2,
                                  mu_1, mu_2, alpha_1, alpha_2, beta, sigma, M, n_locations, T_obs, n_rows_total, m, n_events, n_events_all);
    pos = pos + n_events_all;
  }

  // y^M_t ~ Poisson(L_m(t))
  for (m in 1:M){
    for (t in 1:T_obs){
      for (i in 1:n_locations){
        row_ind = (t-1)*n_locations + i;
        Times_all[m,row_ind] ~ poisson(L[m,row_ind]);
      }
    }
  }
}
