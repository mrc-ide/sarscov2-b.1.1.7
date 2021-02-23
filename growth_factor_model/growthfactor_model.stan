data {
  // sizes
  int<lower=1> N_t;  // total number of timepoints in data
  int<lower=1> N_i;  // total number of areas in data
  int<lower=1> N_age;  // total number of ages in data
  int<lower=1> N_t_par;  // total number of timepoints for model - only for parameter b
  int<lower=1> N_reg;  // total number of regions for model - only for parameter b
  int<lower=1> N_reg_g;  // total number of regions for genome data - normally same as N_reg
  int<lower=1> N_age_par;  // total number of ages for model - only for parameter b
  // data for the poisson/NB models
  int<lower=0> Y_Spositive[N_t, N_i,N_age];  // response variable
  int<lower=0> Y_Snegative[N_t, N_i,N_age];  // response variable
  int<lower=0> Y_Sna[N_t, N_i,N_age];  // response variable
  int<lower=0> genomes_voc[N_t, N_reg_g];  // response variable
  int<lower=0> genomes_nonvoc[N_t, N_reg_g];  // response variable 
  // model covariate indices
  int<lower=1> reg_ind[N_i];  // region index
  int<lower=1> reg_ind_g[N_i];  // region index for genomic data
  int<lower=1> t_ind[N_t];  // week index
  int<lower=1> a_ind[N_age+1];
  // logic variables
  int prior_only;  // should the likelihood be ignored?
  int NB; // if NB==0, sd ~ mean Neg bin, otherwise var~mean
  int Tonly; // whether or not to use T and not s (not in use)
  int infTs; // whether or not to infer T and s 
  int Tminus; //whether or not plus and minus have different generation times
  int Gonly; //only fit genomic data (not S- freq in SGTF). Still models overall pillar 2 counts. -1=fit only SGTF, 0=fit both, 1=fit only genomes
  // Generation time parameters
  real mean_generation_time_T;
  real log_Tg_sd;
  real sq_CV_generation_time_s;
  real log_s_sd;
}

transformed data {
  int<lower=1> N;// = 2*N_i*N_t*N_age;
  real log_first_neg[N_i,N_age];
  real log_first_pos[N_i,N_age];
  int genomes_total[N_t, N_reg_g];
  int Y_total[N_t, N_i,N_age];
  int Y_Sgene[N_t, N_i,N_age];
  // sample size
  if(Gonly==1)
    N = N_i*N_t*N_age + N_t*N_reg_g;
  else if(Gonly==0) 
    N = 2*N_i*N_t*N_age + N_t*N_reg_g;
  else 
    N = 2*N_i*N_t*N_age;
  // get data totals
  for(t in 1:N_t){
    for(a in 1:N_age)
      for(i in 1:N_i){
        Y_total[t,i,a] = Y_Sna[t,i,a] + Y_Spositive[t,i,a] + Y_Snegative[t,i,a];
        Y_Sgene[t,i,a] = Y_Spositive[t,i,a] + Y_Snegative[t,i,a];
      }
    for(i in 1:N_reg_g)
      genomes_total[t,i] = genomes_voc[t,i] +genomes_nonvoc[t,i] ;
  }
  // write initial conditions
  for(a in 1:N_age){
    for(i in 1:N_i){
      if(Y_Sgene[1,i,a]>0){
        log_first_neg[i,a] = log(1 + Y_Snegative[1,i,a] * Y_total[1,i,a] * 1.0/Y_Sgene[1,i,a]);
        log_first_pos[i,a] = log(1 + Y_Spositive[1,i,a] * Y_total[1,i,a] * 1.0/Y_Sgene[1,i,a]);
      }else{
        log_first_neg[i,a] = log(1 + Y_Snegative[1,i,a] );
        log_first_pos[i,a] = log(1 + Y_Spositive[1,i,a] );
      }
    }
  }
}
parameters {
  real rti[N_t-1,N_i,N_age]; // rti - base rate
  real<lower=0> log_growth_Spos_t_sd; // std for base rate
  real b_t[N_t_par,N_reg,N_age_par]; // VOC advantage
  real<lower=0> b_std; // std for first change in b
  real<lower=0> b_t_std; // std for subsequent changes
  real log_T_ratio; // non-VOC generation time / VOC generation time
  // for count models
  real<lower=0> reciprocal_phi; // NB precision parameter
  real<lower=0,upper=1> gamma; // BB precision parameter
  // first obserbations
  real firstobspos[N_i,N_age];
  real firstobsneg[N_i,N_age];
  // generation time
  real<lower=-sq_CV_generation_time_s> s_var;
  real<lower=-mean_generation_time_T> Tg_var;
}

transformed parameters {
  // expected values for count models
  real<upper=20> log_expected_Snegative[N_t, N_i, N_age];
  real<upper=20> log_expected_Spositive[N_t, N_i, N_age];
  real expected_Snegative[N_t, N_i, N_age];
  real expected_Spositive[N_t, N_i, N_age];
  matrix[N_t, N_reg_g] expected_voc;
  matrix[N_t, N_reg_g] expected_nonvoc;
  matrix[N_t, N_reg_g] expected_genome_freq;
  real<upper=20> log_expected_all[N_t, N_i, N_age];
  real expected_all[N_t, N_i, N_age];
  real t1;
  real t2;
  // transformed NB parameter
  real phi = 1. / reciprocal_phi;
  real NB_phi[N_t, N_i, N_age];
  // transformed BB parameters
  real reciprocal_gamma=1.0/gamma-1.0; // BB precision parameter
  real BB_alpha[N_t, N_i, N_age];
  real BB_beta[N_t, N_i, N_age];
  // growth factor variables
  real rminus;
  real exponent_B = 1.0;
  real log_multiplier_A = 0.0;
  // generation time variables
  real Tg = mean_generation_time_T;
  real s = sq_CV_generation_time_s;
  real invTSover7;
  real T_ratio = 1.0;
  // log likelihood to accumulate
  real log_likelihood = 0.0;
  
  // fill in values
  if(Tminus==1) T_ratio = exp(log_T_ratio);
  if(infTs==1){
    Tg += Tg_var; // correction for lognormal mean
    s += s_var;
  }
  invTSover7 = 7.0/(Tg*s);
  // initialise to zero
  for(t in 1:N_t) {
    for(i in 1:N_reg_g) {
      expected_voc[t,i] = 0.0;
      expected_nonvoc[t,i] = 0.0;
    }
  }
  // build expected counts
  for(i in 1:N_i){
    for(a in 1:N_age){
      log_expected_Spositive[1,i,a] = firstobspos[i,a];
      log_expected_Snegative[1,i,a] = firstobsneg[i,a];
      for(t in 1:(N_t-1)){
        if(Tonly==1){
          rminus = (7.0/Tg)*b_t[t_ind[t],reg_ind[i],a_ind[a]] + rti[t,i,a];
        }
        if(Tonly==0){
          exponent_B = exp(s * b_t[t_ind[t],reg_ind[i],a_ind[a]]) / T_ratio;
          log_multiplier_A = exponent_B*invTSover7 - invTSover7 / T_ratio;
          rminus = exponent_B*rti[t,i,a] + log_multiplier_A;
        }
        log_expected_Spositive[t+1,i,a] = log_expected_Spositive[t,i,a] + rti[t,i,a] ;
        log_expected_Snegative[t+1,i,a] = log_expected_Snegative[t,i,a] + rminus ;
      }
    }
  }
  // extract parameters
  expected_Snegative = exp(log_expected_Snegative);
  expected_Spositive = exp(log_expected_Spositive);
  for(t in 1:N_t){
    for(i in 1:N_i){
      for(a in 1:N_age){
        t1 = expected_Spositive[t,i,a];
        t2 = expected_Snegative[t,i,a];
        expected_nonvoc[t,reg_ind_g[i]] += t1;
        expected_voc[t,reg_ind_g[i]] += t2;
        expected_all[t,i,a] = t1+t2;
        NB_phi[t,i,a] = phi*(t1+t2);
        BB_alpha[t,i,a] = reciprocal_gamma * t2/(t1+t2); // expected_freq[t+1,i,a];
        BB_beta[t,i,a] = reciprocal_gamma * t1/(t1+t2); // ( 1.0 - expected_freq[t+1,i,a]);
      }
    }
  }
  log_expected_all = log(expected_all);
  if(NB==0) for(t in 1:N_t) NB_phi[1,,] = rep_array(phi,N_i,N_age); // if NB==0, revert to sd ~ mean Neg bin, otherwise var~mean
  
  // add up log likelihood
  log_likelihood += neg_binomial_2_log_lpmf(to_array_1d(Y_total) | 
           to_array_1d(log_expected_all),to_array_1d(NB_phi));
           
  // if not genome-only model
  if(Gonly != 1)
    log_likelihood += beta_binomial_lpmf(to_array_1d(Y_Snegative) | to_array_1d(Y_Sgene),
      to_array_1d(BB_alpha), to_array_1d(BB_beta));
  
  // if not SGTF-only model
  expected_genome_freq = expected_voc ./(1e-10+expected_voc + expected_nonvoc);
  if(Gonly != -1){
    for(t in 1:N_t)
      log_likelihood += beta_binomial_lpmf(to_array_1d(genomes_voc[t,]) |to_array_1d(genomes_total[t,]) ,
        to_array_1d(reciprocal_gamma*expected_genome_freq[t,]) ,to_array_1d(reciprocal_gamma*(1.0-expected_genome_freq[t,])));
  }
  
}


model {
  // likelihood including all constants
  if(prior_only==0){
    target += log_likelihood;
  }
  // initial condition
  for(a in 1:N_age){
    target += normal_lpdf(firstobspos[,a] | log_first_pos[,a], 0.4) ;
    target += normal_lpdf(firstobsneg[,a] | log_first_neg[,a], 0.4) ;
  }
  // unique growth factor per pair of SGTF observations
  target += exponential_lpdf( log_growth_Spos_t_sd | 10. );
  for(a in 1:N_age){
    target += normal_lpdf(rti[1,,a] | 0.0, .1);
    for(t in 2:(N_t-1)) {
      target += normal_lpdf(rti[t,,a] | rti[t-1,,a],  log_growth_Spos_t_sd);
    }
  }
  // reciprocal precision parameter for NB
  target += exponential_lpdf(reciprocal_phi | 0.1);
  // prior for precision parameter for BB
  target += beta_lpdf(gamma | 1.0, 10.0);
  // std for VOC advantage over time
  target += exponential_lpdf(b_std | 10. );
  target += exponential_lpdf(b_t_std | 10. );
  for(a in 1:N_age_par){
    target += normal_lpdf(b_t[1,,a] | 0,  b_std);
    if(N_t_par>1){
      for(t in 2:(N_t_par)) {
        target += normal_lpdf(b_t[t,,a] | b_t[t-1,,a], b_t_std);
      }
    }
  }
  // priors for generation time values
  target += normal_lpdf(log_T_ratio | 0.0, 0.1);
  s_var ~ normal( 0.0, log_s_sd) T[-sq_CV_generation_time_s,];
  Tg_var ~ normal(0.0, log_Tg_sd) T[-mean_generation_time_T,];
}

generated quantities {
  // generate predictions for return
  real S_negative_pr[N_t, N_i,N_age];
  real S_positive_pr[N_t, N_i,N_age];
  real S_all_pr[N_t, N_i,N_age];
  real genomes_voc_pr[N_t, N_reg_g];  
  real genomes_nonvoc_pr[N_t, N_reg_g]; 
  real WBIC = log_likelihood;
  vector[N_t*(N_i+N_reg_g)] log_lik;
  real log_lik_2d[N_t, N_i+N_reg_g];  
  simplex[N_age] pi_negative[N_t,N_i];
  simplex[N_age] pi_positive[N_t,N_i];
  
  for(i in 1:N_i){
    for(t in 1:N_t){
      log_lik_2d[t,i] = 0;
      for(a in 1:N_age){
        S_negative_pr[t,i,a] = beta_binomial_rng(Y_Sgene[t,i,a],BB_alpha[t,i,a],BB_beta[t,i,a]);
        S_positive_pr[t,i,a] = Y_Sgene[t,i,a]-S_negative_pr[t,i,a];
        S_all_pr[t,i,a] = neg_binomial_2_log_rng(log_expected_all[t,i,a],NB_phi[t,i,a]);
        log_lik_2d[t,i] += neg_binomial_2_log_lpmf(Y_total[t,i,a] |log_expected_all[t,i,a],NB_phi[t,i,a]);
        if(Gonly != 1) {
          log_lik_2d[t,i] += beta_binomial_lpmf(Y_Snegative[t,i,a]|Y_Sgene[t,i,a],
              BB_alpha[t,i,a],BB_beta[t,i,a]);
        }
      }
    }
    for(t in 1:N_t){
      pi_negative[t,i] = softmax( to_vector(log_expected_Snegative[t,i,]) );
      pi_positive[t,i] = softmax( to_vector(log_expected_Spositive[t,i,]) );
    }
  }
  for(t in 1:N_t) {
    for(i in 1:N_reg_g) {
        log_lik_2d[t,i+N_i] = 0;
        if(Gonly != -1)
            log_lik_2d[t,i+N_i] += beta_binomial_lpmf(genomes_voc[t,i]|genomes_total[t,i],
                reciprocal_gamma*expected_genome_freq[t,i],reciprocal_gamma*(1.0-expected_genome_freq[t,i]));
        genomes_voc_pr[t,i]= beta_binomial_rng(genomes_voc[t,i]+genomes_nonvoc[t,i],
            reciprocal_gamma*expected_genome_freq[t,i],reciprocal_gamma*(1.0-expected_genome_freq[t,i]));
        genomes_nonvoc_pr[t,i]=genomes_voc[t,i]+genomes_nonvoc[t,i]-genomes_voc_pr[t,i];
    }
  }

  // cast log likelihood into 1 dim format
  log_lik = to_vector( to_array_1d( log_lik_2d ) );
}
