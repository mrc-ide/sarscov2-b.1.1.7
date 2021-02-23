// generated with brms 2.14.4
functions {
  real sigmoid(real x){
  return 1/(1+exp(-x));
  }  

}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_negative;  // number of observations
  int Y_negative[N_negative];  // response variable
  int trials_negative[N_negative];  // number of trials
  int<lower=1> Ksp_negative;  // number of special effects terms
  int<lower=1> N_Rt;  // number of observations
  vector[N_Rt] Y_Rt;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_Rt] noise_Rt;
  // information about non-missings
  int<lower=0> Nme_Rt;
  int<lower=1> Jme_Rt[Nme_Rt];
  int<lower=1> K_Rt;  // number of population-level effects
  matrix[N_Rt, K_Rt] X_Rt;  // population-level design matrix
  int<lower=1> Ksp_Rt;  // number of special effects terms
  // data for splines
  int Ks_Rt;  // number of linear effects
  matrix[N_Rt, Ks_Rt] Xs_Rt;  // design matrix for the linear effects
  // data for spline s(scale(epiweek),k=4)
  int nb_Rt_1;  // number of bases
  int knots_Rt_1[nb_Rt_1];  // number of knots
  // basis function matrices
  matrix[N_Rt, knots_Rt_1[1]] Zs_Rt_1_1;
  int<lower=1> N_fr;  // number of observations
  //vector[N_fr] Y_fr;  // response variable
  int<lower=0> Nmi_fr;  // number of missings
  int<lower=1> Jmi_fr[Nmi_fr];  // positions of missings
  int<lower=1> K_fr;  // number of population-level effects
  matrix[N_fr, K_fr] X_fr;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
  int additive;  // 1=additive model, 0=multiplicative
}
transformed data {
  int Kc_Rt = K_Rt - 1;
  matrix[N_Rt, Kc_Rt] Xc_Rt;  // centered version of X_Rt without an intercept
  vector[Kc_Rt] means_X_Rt;  // column means of X_Rt before centering
  int Kc_fr = K_fr - 1;
  matrix[N_fr, Kc_fr] Xc_fr;  // centered version of X_fr without an intercept
  vector[Kc_fr] means_X_fr;  // column means of X_fr before centering
  for (i in 2:K_Rt) {
    means_X_Rt[i - 1] = mean(X_Rt[, i]);
    Xc_Rt[, i - 1] = X_Rt[, i] - means_X_Rt[i - 1];
  }
  for (i in 2:K_fr) {
    means_X_fr[i - 1] = mean(X_fr[, i]);
    Xc_fr[, i - 1] = X_fr[, i] - means_X_fr[i - 1];
  }
}
parameters {
  vector[N_Rt] Yl_Rt;  // latent variable
  vector[Kc_Rt] b_Rt;  // population-level effects
  real Intercept_Rt;  // temporary intercept for centered predictors
  vector[Ksp_Rt] bsp_Rt;  // special effects coefficients
  vector[Ks_Rt] bs_Rt;  // spline coefficients
  // parameters for spline s(scale(epiweek),k=4)
  // standarized spline coefficients
  vector[knots_Rt_1[1]] zs_Rt_1_1;
  real<lower=0> sds_Rt_1_1;  // standard deviations of spline coefficients
  real<lower=0> sigma_Rt;  // residual SD
  vector[Nmi_fr] Ymi_fr;  // estimated missings
  vector[Kc_fr] b_fr;  // population-level effects
  real Intercept_fr;  // temporary intercept for centered predictors
  real<lower=0> sigma_fr;  // residual SD
}
transformed parameters {
  // actual spline coefficients
  vector[knots_Rt_1[1]] s_Rt_1_1;
  // compute actual spline coefficients
  s_Rt_1_1 = sds_Rt_1_1 * zs_Rt_1_1;
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // vector combining observed and missing responses
    //vector[N_fr] Yl_fr = Y_fr;
    // initialize linear predictor term
    vector[N_negative] mu_negative = rep_vector(0.0, N_negative);
    // initialize linear predictor term
    vector[N_Rt] mu_Rt = Intercept_Rt + Xc_Rt * b_Rt + Xs_Rt * bs_Rt + Zs_Rt_1_1 * s_Rt_1_1;
    // initialize linear predictor term
    vector[N_fr] mu_fr = Intercept_fr + Xc_fr * b_fr;
    //Yl_fr[Jmi_fr] = Ymi_fr;
    for (n in 1:N_negative) {
      // add more terms to the linear predictor
      mu_negative[n] +=  Ymi_fr[n];
    }
    // add more terms to the linear predictor
    if(additive){
      for (n in 1:N_Rt) {
        mu_Rt[n] += bsp_Rt[1] * sigmoid(Ymi_fr[n]);
      }
    }else{
      for (n in 1:N_Rt) {
        mu_Rt[n] *= (1 - sigmoid(Ymi_fr[n]) + bsp_Rt[1] * sigmoid(Ymi_fr[n]));
      }
    }
    target += binomial_logit_lpmf(Y_negative | trials_negative, mu_negative);
    target += normal_lpdf(Yl_Rt | mu_Rt, sigma_Rt);
    target += normal_lpdf(Ymi_fr | mu_fr, sigma_fr);
  }
  // priors including all constants
  target += normal_lpdf(Y_Rt[Jme_Rt] | Yl_Rt[Jme_Rt], noise_Rt[Jme_Rt]);
  for(n in 1:Kc_Rt)
    target += normal_lpdf(b_Rt[n] | 0, .1);
  target += normal_lpdf(Intercept_Rt | 0, .5);
  target += std_normal_lpdf(sds_Rt_1_1 );
  target += std_normal_lpdf(zs_Rt_1_1);
  target += std_normal_lpdf(bs_Rt);
  target += std_normal_lpdf(bsp_Rt);
  target += exponential_lpdf(sigma_Rt | 20.);
  for(n in 1:Kc_fr)
    target += normal_lpdf(b_fr[n] | 0, .5);
  target += std_normal_lpdf(Intercept_fr );
  target += exponential_lpdf(sigma_fr | 10.);
}
generated quantities {
  // actual population-level intercept
  real b_Rt_Intercept = Intercept_Rt - dot_product(means_X_Rt, b_Rt);
  // actual population-level intercept
  real b_fr_Intercept = Intercept_fr - dot_product(means_X_fr, b_fr);
    vector[N_Rt] mu_Rt = Intercept_Rt + Xc_Rt * b_Rt + Xs_Rt * bs_Rt + Zs_Rt_1_1 * s_Rt_1_1;
    // add more terms to the linear predictor
    if(additive){
      for (n in 1:N_Rt) {
        mu_Rt[n] += bsp_Rt[1] * sigmoid(Ymi_fr[n]);
      }
    }else{
      for (n in 1:N_Rt) {
        mu_Rt[n] *= (1 - sigmoid(Ymi_fr[n]) + bsp_Rt[1] * sigmoid(Ymi_fr[n]));
      }
    }
}