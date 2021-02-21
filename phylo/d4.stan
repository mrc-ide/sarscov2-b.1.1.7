data {
  int<lower=1> N_t;  // total number of timepoints
  int<lower=1> N_i;  // total number of areas
  int<lower=1> N_w;  // total number of weeks
  // for the binomial model
  int<lower=1> N_B117;  // number of observations
  int<lower=0> Y_B117[N_t,N_i];  // response variable
  int<lower=0> Y_total[N_t,N_i];  // offset variable
  real first_offset;
  // covariate indices
  int<lower=1> t_to_w[N_t];  // time to week
  // logic variables
  int<lower=0, upper=1> prior_only;  // should the likelihood be ignored?
}
transformed data {
  int N_t_m1 = N_t - 1;
  int N_t_m2 = N_t - 2;
}
parameters {
  // covariates
  real firstobs[N_i];
  real <lower=0> rho_sd;          // random effect growth diff sd
  real overall_rho_0[N_t_m1];      // change in overall transmission advantage (zero-centered)
  real rho_0[N_t_m1,N_i];   //  regional transmission advantage (difference from overall, zero-centered)
  real<lower=0> overall_rho_t_sd;      // overall growth difference variance at following weeks
  
}
transformed parameters {
  real overall_rho[N_t_m1];      // overall transmission advantage
  real rho[N_t_m1,N_i];     // random effect transmission advantage
  real<lower=0, upper=1> freq[N_w,N_i];
  real phi[N_w,N_i];
  overall_rho[1] = overall_rho_0[1];
  for(t in 2:N_t_m1)
    overall_rho[t] = overall_rho_0[t] + overall_rho[t-1] ;
  
  for(i in 1:N_i){
    rho[1,i] = rho_0[1,i] + overall_rho[1];
    for(t in 2:N_t_m1)
      rho[t,i] = rho_0[t,i] + overall_rho[t] ;
  }
  // frequency in log-odds for each week and region 
  for(i in 1:N_i){
    phi[1,i] = firstobs[i] + first_offset;
    for(w in 2:(N_w)){
      phi[w,i] = phi[w-1,i] + rho[w-1,i];
    }
  }
  freq = inv_logit(phi);
}
model {
  // likelihood including all constants
  if (!prior_only) {
    for(t in 1:N_t){
      for(i in 1:N_i){
        if(Y_total[t,i] > 0){
          target += binomial_lpmf(Y_B117[t,i] | Y_total[t,i],freq[t_to_w[t],i]);
        }
      }
    }
  }
  for(i in 1:N_i){
    target += normal_lpdf(firstobs[i] | 0, 10.0); //2.0
  }
  // 
  // prior on overall growth difference
  // autocorrelated in time
  // a priori centered at no advantage
  
  // initial rho
  target += normal_lpdf( overall_rho_0[1] | 0.0, 1.0);
  
  // change in rho variance param
  target += exponential_lpdf(overall_rho_t_sd | 10. ); // expected  small drift in growth rate difference per week 
  
  // rw prior for overall rho (_0 is centered)
  target += normal_lpdf(  to_array_1d(overall_rho_0[2:N_t_m1]) |  0.0,  overall_rho_t_sd );
  
  // cor regional and overall growth, variance parameter
  target += exponential_lpdf( rho_sd | 1. ); 
  // cor with overall rho (_0 is offset from overall rho)
  target += normal_lpdf( to_array_1d(rho_0) |  0.0,  rho_sd );
  
  // rw prior for regional rho 
  target += normal_lpdf(
    to_array_1d(rho[2:N_t_m1,]) |
    to_array_1d(rho[1:N_t_m2,]),
    overall_rho_t_sd );
}
generated quantities {
  int B_count[N_t,N_i];
  for(t in 1:N_t)
    for(i in 1:N_i)
      B_count[t,i] = binomial_rng(Y_total[t,i],freq[t_to_w[t],i]);
  
}
