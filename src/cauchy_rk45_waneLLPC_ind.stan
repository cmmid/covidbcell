functions {
  vector myODE_waneLLPB(real t, vector y, real param1, real param2, real param3, real param4, real prop_asc,
    real deltaAg, real deltaM, real deltaP, real deltaL, real deltaA, real deltaG) {
    vector[6] dydt;  // Array to store the derivatives
    // Define the ODEs
    dydt[1] = exp(- deltaAg * t) * param1 - exp(- deltaAg * t) * y[1] * param2 - y[1] * deltaM ;  // memory b-cell conc.
    dydt[2] = exp(- deltaAg * t) * y[1] * param2 * prop_asc - y[2] * deltaP ;  // plasmablast conc.
    dydt[3] = y[2] * param3 - y[3] * deltaA ;  // antibodies from plasmablasts
    dydt[4] = exp(- deltaAg * t) * y[1] * param2 * (1 - prop_asc) - y[4] * deltaG; // germinal center
    dydt[5] = y[4] * deltaG - y[5] * deltaL;  //  long-lived plasma cells
    dydt[6] = y[5] * param4 - y[6] * deltaA;  //  antibodies from long-lived plasma cells
    return dydt;
  }
}
data {
  int<lower=1> N;  // Number of observations
  array[N] real N_id; // IDs for teh number of observations

  int<lower=1> N_vac;
  array[N] int x_vac;

  int<lower=1> N_time;
  array[N] int x_time;

  int<lower=1> N_age;
  array[N] int x_age;

  vector[N_vac] prop_vac;
  vector[N_time] prop_time;
  vector[N_age] prop_age;

  int T;
  array[N, T] int y_times;  // Array of observation times
  array[N] int t_max ;
  array[N] int t_length ;

  array[N] vector[6] y0_init;  // Initial values of y (state variables)
  array[N, T, 3] real y_obs;  // Array of observed values of y (state variables)

  array[N_vac] vector[6] y0_init_vac;
  array[N_time] vector[6] y0_init_time;
  array[N_age] vector[6] y0_init_age;

  //int<lower = 0, upper = 1> flag_llps_ab;
  //int<lower = 0, upper = 1> flag llps_decay;

}
transformed data {
  array[250] real t_full = linspaced_array(250, 1, 250);
  real relative_tolerance = 1e-6;
  real absolute_tolerance = 1e-5;
  int max_num_steps = 100000;
  int min_steps_big = 750;
  int min_steps_small = 250;
  int N_ll = 0;
  for (j in 1:N) {
    for (i in 1:t_length[j]) {
      // on observed dynamics
      if (y_obs[j, i, 1] > -1) {
        N_ll = N_ll + 1;
      }
      if (y_obs[j, i, 2] > -1) {
        N_ll = N_ll + 1;
      }
      if (y_obs[j, i, 3] > -1) {
        N_ll = N_ll + 1;
      }
    }
  }
  real deltaM = 1000;
}
parameters {
    real<lower=-10, upper = 10> param1;  // Parameters to be estimated
    real<lower=-10, upper = 10> param2;  // Parameters to be estimated
    real<lower=-10, upper = 10> param3;  // Parameters to be estimated
    real<lower=-10, upper = 10> param4;  // Parameters to be estimated

    vector[N_vac] v_1;  // Parameters to be estimated
    vector[N_vac] v_3;  // Parameters to be estimated
    vector[N_vac] v_4;  // Parameters to be estimated

    vector[N_time] t_1;  // Parameters to be estimated
    vector[N_time] t_3;  // Parameters to be estimated
    vector[N_time] t_4;  // Parameters to be estimated

    vector[N_age] a_1;  // Parameters to be estimated
    vector[N_age] a_3;  // Parameters to be estimated
    vector[N_age] a_4;  // Parameters to be estimated

    real<lower = 0>  sigmav1; 
    real<lower = 0>  sigmat1; 
    real<lower = 0>  sigmaa1; 
    real<lower = 0>  sigmav3; 
    real<lower = 0>  sigmat3; 
    real<lower = 0>  sigmaa3; 
    real<lower = 0>  sigmav4; 
    real<lower = 0>  sigmat4; 
    real<lower = 0>  sigmaa4; 

    vector[N] ind_1;
    vector[N] ind_3;
    vector[N] ind_4;
    real<lower = 0>  sigmai1; 
    real<lower = 0>  sigmai3; 
    real<lower = 0>  sigmai4; 


    real<lower = 1.0 / 30, upper = 1> deltaAg;
    real<lower = 1.0, upper =6> deltaP;
    real<lower = 365, upper = 365 * 10> deltaL;
    real<lower = 1, upper = 50> deltaA;
    real<lower = 10, upper = 28> deltaG;

    real<lower = 0, upper = 1> prop_asc;

    real<lower = 0> sigma_M;
    real<lower = 0> sigma_P;
    real<lower = 0> sigma_A;
}
transformed parameters {
  array[N, T] vector[6] y_hat;  // Array to store the predicted values of y


  // Solve the ODE system using the initial conditions and parameter values
  
    array[11] real theta;
    theta[5] = prop_asc;

    theta[6] = 1.0/(deltaAg * 30);
    theta[7] = 1.0/(deltaM);
    theta[8] = 1.0/(deltaP);
    theta[9] = 1.0/(deltaL);
    theta[10] = 1.0/(deltaA);
    theta[11] = 1.0/(deltaG);

    //array[N] real ts;
    for (i in 1:N) { // individuals 

      theta[1] = inv_logit(param1 + v_1[x_vac[i]] * sigmav1 + t_1[x_time[i]] * sigmat1 + a_1[x_age[i]] * sigmaa1 + ind_1[i] * sigmai1) * 2;
      theta[2] = inv_logit(param2); 
      theta[3] = inv_logit(param3 + v_3[x_vac[i]] * sigmav3 + t_3[x_time[i]] * sigmat3 + a_3[x_age[i]] * sigmaa3 + ind_3[i] * sigmai3) * 6;
      theta[4] = inv_logit(param4 + v_4[x_vac[i]] * sigmav4 + t_4[x_time[i]] * sigmat4 + a_4[x_age[i]] * sigmaa4 + ind_4[i] * sigmai4) * 0.5;

      
      if (
        (theta[1] > 0 && theta[1] < 20) &&
        (theta[2] > 0 && theta[2] < 1) && 
        (theta[3] > -1 && theta[3] < 20) &&
        (theta[4] > -1 && theta[4] < 20) ) {

        y_hat[i, 1:t_length[i]] = ode_rk45(myODE_waneLLPB, y0_init[i], 0, y_times[i, 1:t_length[i]],
          theta[1], theta[2], theta[3], theta[3] * theta[4], theta[5], theta[6], theta[7], theta[8], theta[9], theta[10], theta[11]);
        
      } else{
        break;
      }
    }
}
model {
  // Priors on parameters
  param1 ~ normal(0, 1.814);
  param2 ~ normal(0, 1.814);
  param3 ~ normal(0, 1.814);
  param4 ~ normal(0, 1.814);

  prop_asc ~ beta(60, 40); // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2211779

  deltaAg ~ uniform(1.0 / 30, 1);  // Specify priors for deltaAg (10 days antigen in body)
  deltaP ~ normal(4, 1);  // Specify priors for deltaP (~2 days) https://www.frontiersin.org/articles/10.3389/fimmu.2019.00721/full
  deltaA ~ normal(30, 5);  // Specify priors for deltaA (half life 100 days)
  deltaG ~ normal(18, 2);  // Prior on time in germinal center (~ 2 weeks)
  deltaL ~ normal(730, 200);  // Specify priors for deltaP (~2 days)

  v_1 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmav1 ~ exponential(3);
  v_3 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmav3 ~ exponential(3);
  v_4 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmav4 ~ exponential(3);

  t_1 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmat1 ~ exponential(3);
  t_3 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmat3 ~ exponential(3);
  t_4 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmat4 ~ exponential(3);

  a_1 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmaa1 ~ exponential(3);
  a_3 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmaa3 ~ exponential(3);
  a_4 ~ normal(0, 1);  // Specify priors for B-cell profileration
  sigmaa4 ~ exponential(3);

  ind_1 ~ normal(0, 1);
  ind_3 ~ normal(0, 1);
  ind_4 ~ normal(0, 1);
  sigmai1 ~ exponential(3);
  sigmai3 ~ exponential(3);
  sigmai4 ~ exponential(3);

  sigma_M ~ exponential(1);
  sigma_P ~ exponential(1);
  sigma_A ~ exponential(1);


  if (theta[1] < 0 || theta[2] < -1 || theta[3] < -1  || theta[4] < -1 ||
    theta[1] > 20 || theta[2] > 20 || theta[3] > 20 || theta[4] > 20 ) {
    target += -2147483647;
  } else {
  // Likelihood
    for (j in 1:N) {
      for (i in 1:t_length[j]) {
        // on observed dynamics
        if (y_obs[j, i, 1] > -1) {y_obs[j, i, 1] ~ normal(y_hat[j, i, 1], sigma_M);};  // Specify the likelihood distribution for y_obs given y_hat
        if (y_obs[j, i, 2] > -1) {y_obs[j, i, 2] ~ normal(y_hat[j, i, 2], sigma_P);};  // Specify the likelihood distribution for y_obs given y_hat
        if (y_obs[j, i, 3] > -1) {y_obs[j, i, 3] ~ normal(y_hat[j, i, 3] + y_hat[j, i, 6], sigma_A);};  // Specify the likelihood distribution for y_obs given y_hat
      }
    }
  }
}
generated quantities {

  array[N_vac, 250] vector[6] y_full_vac;  // Array to store the predicted values of y
  array[N_time, 250] vector[6] y_full_time_cha;  // Array to store the predicted values of y
  array[N_time, 250] vector[6] y_full_time_bnt;  // Array to store the predicted values of y
  array[N_age, 250] vector[6] y_full_age_cha;  // Array to store the predicted values of y
  array[N_age, 250] vector[6] y_full_age_bnt;  // Array to store the predicted values of y

  array[4] real theta_main;
  array[N_vac, 4] real theta_v;
  array[N_time, 2, 4] real theta_t;
  array[N_age, 2, 4] real theta_a;

  theta_main = rep_array(0, 4);
  theta_v = rep_array(0, N_vac, 4);
  theta_t = rep_array(0, N_time, 2, 4);
  theta_a = rep_array(0, N_age, 2, 4);

  real deltaAg_i;
  real deltaM_i;
  real deltaP_i;
  real deltaL_i;
  real deltaA_i;
  real deltaG_i;

  deltaAg_i = 1.0/(deltaAg * 30);
  deltaM_i = 1.0/(deltaM);
  deltaP_i = 1.0/(deltaP);
  deltaL_i = 1.0/(deltaL);
  deltaA_i = 1.0/(deltaA);
  deltaG_i = 1.0/(deltaG);

  real temp_1;
  real temp_2;
  real temp_3;
  real temp_4;

  // Marginal posterior distirbutions (to compare with data)

  for (i in 1:N_vac) {
    for (j in 1:N_time) {
      for (k in 1:N_age) {
        temp_1 = inv_logit(param1 + v_1[i] * sigmav1 + t_1[j] * sigmat1 + a_1[k] * sigmaa1) * 2;
        temp_2 = inv_logit(param2);
        temp_3 = inv_logit(param3 + v_3[i] * sigmav3 + t_3[j] * sigmat3 + a_3[k] * sigmaa3) * 6;
        temp_4 = inv_logit(param4 + v_4[i] * sigmav4 + t_4[j] * sigmat4 + a_4[k] * sigmaa4) * 0.5;

        theta_main[1] += temp_1 * prop_vac[i] * prop_time[j] * prop_age[k];
        theta_main[2] += temp_2 * prop_vac[i] * prop_time[j] * prop_age[k];
        theta_main[3] += temp_3 * prop_vac[i] * prop_time[j] * prop_age[k];
        theta_main[4] += temp_4 * prop_vac[i] * prop_time[j] * prop_age[k];

        theta_v[i, 1] += temp_1 * prop_time[j] * prop_age[k];
        theta_v[i, 2] += temp_2 * prop_time[j] * prop_age[k];
        theta_v[i, 3] += temp_3 * prop_time[j] * prop_age[k];
        theta_v[i, 4] += temp_4 * prop_time[j] * prop_age[k];

        theta_t[j, i, 1] += temp_1 * prop_age[k];
        theta_t[j, i, 2] += temp_2 * prop_age[k];
        theta_t[j, i, 3] += temp_3 * prop_age[k];
        theta_t[j, i, 4] += temp_4 * prop_age[k];

        theta_a[k, i, 1] += temp_1 * prop_time[j];
        theta_a[k, i, 2] += temp_2 * prop_time[j];
        theta_a[k, i, 3] += temp_3 * prop_time[j];
        theta_a[k, i, 4] += temp_4 * prop_time[j];

      }
    }
  }

  for (i in 1:N_vac) {
      profile("ODE_gp") {
        y_full_vac[i] = ode_rk45(myODE_waneLLPB, y0_init_vac[i], 0, t_full,
          theta_v[i, 1], theta_v[i, 2], theta_v[i, 3], theta_v[i, 3]*theta_v[i, 4], prop_asc, deltaAg_i, deltaM_i, deltaP_i, deltaL_i, deltaA_i, deltaG_i);
      }

  }

    for (j in 1:N_time) {

      profile("ODE_gp") {
      y_full_time_cha[j] = ode_rk45(myODE_waneLLPB, y0_init_time[j], 0, t_full,
        theta_t[j, 1, 1], theta_t[j, 1, 2], theta_t[j, 1, 3], theta_t[j, 1, 3]*theta_t[j, 1, 4], prop_asc, deltaAg_i, deltaM_i, deltaP_i,  deltaL_i, deltaA_i, deltaG_i);
      y_full_time_bnt[j] = ode_rk45(myODE_waneLLPB, y0_init_time[j], 0, t_full,
        theta_t[j, 2, 1], theta_t[j, 2, 2], theta_t[j, 2, 3], theta_t[j, 2, 3]*theta_t[j, 2, 4], prop_asc, deltaAg_i, deltaM_i, deltaP_i,  deltaL_i, deltaA_i, deltaG_i);
      }
  }


    for (k in 1:N_age) {
      profile("ODE_gp") {
      y_full_age_cha[k] = ode_rk45(myODE_waneLLPB, y0_init_age[k], 0, t_full,
        theta_a[k, 1, 1], theta_a[k,1, 2], theta_a[k, 1, 3], theta_a[k, 1, 3]*theta_a[k, 1, 4], prop_asc, deltaAg_i, deltaM_i, deltaP_i,  deltaL_i, deltaA_i, deltaG_i);
      y_full_age_bnt[k] = ode_rk45(myODE_waneLLPB, y0_init_age[k], 0, t_full,
        theta_a[k, 2, 1], theta_a[k, 2, 2], theta_a[k, 2, 3], theta_a[k, 2, 3]*theta_a[k, 2, 4], prop_asc, deltaAg_i, deltaM_i, deltaP_i, deltaL_i,  deltaA_i, deltaG_i);
      }
  }

  vector[N_ll] log_lik;
  log_lik = rep_vector(0, N_ll);
  int k = 1;
  for (j in 1:N) {
    for (i in 1:t_length[j]) {
      // on observed dynamics
      if (y_obs[j, i, 1] > -1) {
        log_lik[k] = normal_lpdf(y_obs[j, i, 1] | y_hat[j, i, 1], sigma_M);
        k = k + 1;
      }
      if (y_obs[j, i, 2] > -1) {
        log_lik[k] = normal_lpdf(y_obs[j, i, 2] | y_hat[j, i, 2], sigma_P);
        k = k + 1;
      }
      if (y_obs[j, i, 3] > -1) {
        log_lik[k] = normal_lpdf(y_obs[j, i, 3] | y_hat[j, i, 3] + y_hat[j, i, 6], sigma_A);
        k = k + 1;
      }
    }
  }
}

