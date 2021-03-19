// Model for FSD data, marginal model for Z

functions {
  
  
  /* weighted median, ported from .wquant function from loo package
  *   Vehtari A, Gabry J, Magnusson M, Yao Y, Bürkner P, Paananen T, Gelman A (2020). 
  “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” 
  R package version 2.3.1, <URL: https://mc-stan.org/loo>.
  * Args:
  *   xo: vector of values
  *   wo: weights
  * Returns:
  *   a weighted median
  */
  real weighted_median(vector xo, vector wo) {
    int n = num_elements(xo);
    real median = 0.0;
    int id = 0;
    int ord[n] = sort_indices_asc(xo);
    vector[n] x = xo[ord];
    vector[n] w = cumulative_sum(wo[ord]);
    w = w / w[n];
    for(i in 1:n) {
      if(w[i] >= 0.5) {
        id = i;
        break;
      }
    }
    if (id == 1) {
      median = x[1];
    } else {
      real w1 = w[id - 1];
      real x1 = x[id - 1];
      median = x1 + (x[id] - x1) * (0.5 - w1)/(w[id] - w1);
    }
    return(median);
  }
  /* monotonic effect and sratio pmf copied and modified from brms R package
  *
  *  Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using
  *  Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
  *
  *  Paul-Christian Bürkner (2018). Advanced Bayesian Multilevel Modeling with the R Package
  *  brms. The R Journal, 10(1), 395-411. doi:10.32614/RJ-2018-017
  */
  
  /* compute monotonic effects
  * Args:
  *   scale: a simplex parameter
  *   i: index to sum over the simplex
  * Returns:
  *   a scalar between 0 and 1
  */
  real mo(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return rows(scale) * sum(scale[1:i]);
    }
  }
  
  /* sratio-logit log-PDF for a single response
  * Args:
  *   x: response category (0-indexed)
  *   mu: linear predictor
  *   thres: ordinal thresholds
  * Returns:
  *   a scalar to be added to the log posterior
  */
  real sratio_logit_lpmf(int x, real mu, vector thres) {
    int y = x + 1;
    int ncat = num_elements(thres) + 1;
    vector[ncat] p;
    vector[ncat - 1] q;
    int k = 1;
    while (k <= min(y, ncat - 1)) {
      q[k] = 1 - inv_logit(thres[k] - mu);
      p[k] = 1 - q[k];
      for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];
      k += 1;
    }
    if (y == ncat) {
      p[ncat] = prod(q);
    }
    return log(p[y]);
  }
  
  /* sratio-logit rng, uses lpmf, not optimal but fast enough
  * Args:
  *   mu: linear predictor
  *   thres: ordinal thresholds
  * Returns:
  *   response category (0-indexed)
  */
  int sratio_logit_rng(real mu, vector thres) {
    
    int ncat = num_elements(thres) + 1;
    real u = uniform_rng(0, 1);
    int k = 0;
    real p = exp(sratio_logit_lpmf(k | mu, thres));
    while (p < u && k < (ncat - 1)) {
      k += 1;
      p += exp(sratio_logit_lpmf(k | mu, thres));
    }
    return k;
  }
  
}
data {
  int<lower=1> n;
  // Y, income
  vector<lower=0> [n] income;         
  // X, education levels from 0 to 2 (lowest to highest)
  int<lower=0, upper=2> education[n];
  // Z, grade average from elementary school scaled to [0,1]
  vector<lower=0,upper=1>[n] grade;
  // W, sosioeconomic status of parents from 1 to 3 (lowest to highest)
  int<lower=1, upper=3> ses[n];    
  // I, ITPA scores
  vector[n] itpa;
  // G, 0 = Male, 1 = Female
  int<lower=0, upper=1> gender[n];   
  // Number of Monte Carlo samples per S and G (outer loop)
  int<lower=1> N; 
  // Number of Monte Carlo samples per W (inner loop)
  int<lower=1> M;
}


parameters {
  
  // ordered logistic model for W
  // low sosioeconomic status as reference
  ordered[2] a_w;
  
  // sequential model for X
  // secondary level as reference
  vector[2] a_x; 
  real b_xg;
  real b_xs;
  real b_xz;
  real b_xw;
  simplex[2] simplex_xw;
  
  // ITPA conditioned on SES
  real a_sw;
  real b_sw;
  simplex[2] simplex_sw;
  real<lower=0> sigma_sw;

  // Gender as Bernoulli
  real<lower=0,upper=1> p_gender;
  
  // Income as Gamma
  real a_y; 
  real b_yg;
  real b_ys;
  real b_yz;
  real b_yw;
  simplex[2] simplex_yw;
  real b_yx;
  simplex[2] simplex_yx;
  real<lower=0> shape_y;
  
  // Beta distribution for Z
  real a_z_phi;
  real a_z;
}
model {
  
  for(i in 1:n) {
    
    real mu_education = b_xg * gender[i] + b_xs * itpa[i] + b_xz * grade[i] + 
    b_xw * mo(simplex_xw, ses[i] - 1);
    
    real mu_income = a_y + b_yg * gender[i] + b_ys * itpa[i] + b_yz * grade[i] + 
    b_yx * mo(simplex_yx, education[i]) + b_yw * mo(simplex_yw, ses[i] - 1);
    
    
    
    ses[i] ~ ordered_logistic(0, a_w);
    education[i] ~ sratio_logit(mu_education, a_x);
    income[i] ~ gamma(shape_y, shape_y * exp(-mu_income));
    itpa[i] ~ normal(a_sw + b_sw * mo(simplex_sw, ses[i] - 1), sigma_sw);
  }
  {
    real mu_grade = inv_logit(a_z);
    real phi_grade = exp(a_z_phi);
    grade ~ beta(mu_grade * phi_grade, (1 - mu_grade) * phi_grade);
  }
  gender ~ bernoulli(p_gender);
}


generated quantities {
  
  // causal part
  
  // expected value of Y|do(X=x)
  vector[3] mean_ydox;
  vector[3] median_ydox;
  
  vector[3] mcse_mean;
  vector[3] mcse_median;
  
  // sample from Y|do(X=x)
  vector[3] sample_ydox;
  
  real mu_grade = inv_logit(a_z);
  real phi_grade = exp(a_z_phi);
  // each do_x value
  for(k in 1:3) {  
    
    vector[M * N] sample;
    vector[M * N] weight;
    
    for(i in 1:N) {
       // sample gender and ITPA
      int gender_i = bernoulli_rng(p_gender);
      //sample itpa and ses (and ignore ses i.e. get marginal distribution of itpa)
      int ses_i = ordered_logistic_rng(0, a_w);
      real itpa_i = normal_rng(a_sw + b_sw * mo(simplex_sw, ses_i - 1), sigma_sw);
   
      real trapdoor_grade = beta_rng(mu_grade * phi_grade, (1 - mu_grade) * phi_grade);
      
      vector[M] weight_jk;
      vector[M] sample_jk;
      
      for(j in 1:M) {
        // sample sosioeconomic status of parents
        int ses_j = ordered_logistic_rng(0, a_w);
        
        real mu_education = b_xg * gender_i + b_xs * itpa_i + b_xz * trapdoor_grade +
        b_xw * mo(simplex_xw, ses_j - 1);
        
        real mu_income = a_y + b_yg * gender_i + b_ys * itpa_i + b_yz * trapdoor_grade +
        b_yx * mo(simplex_yx, k - 1) + b_yw * mo(simplex_yw, ses_j - 1);
        
        sample_jk[j] = gamma_rng(shape_y, shape_y * exp(-mu_income));
        weight_jk[j] = sratio_logit_lpmf(k - 1 | mu_education, a_x) +
        normal_lpdf(itpa_i | a_sw + b_sw * mo(simplex_sw, ses_j - 1), sigma_sw);
      }
      
      // normalize the weights
      weight_jk = exp(weight_jk - max(weight_jk));
      // scale with N so that the total sum is 1
      weight_jk =  weight_jk / sum(weight_jk) / N; 
      
      weight[((i-1)*M + 1):(i*M)] = weight_jk;
      sample[((i-1)*M + 1):(i*M)] = sample_jk;
    }
    // compute mean, median and sample one replication
    
    mean_ydox[k] = dot_product(sample, weight);
    median_ydox[k] = weighted_median(sample, weight);
    mcse_mean[k] = sqrt(sum(square(weight .* (sample - mean_ydox[k]))));
    mcse_median[k] = sqrt(sum(square(weight .* (sample - median_ydox[k]))));
    sample_ydox[k] = sample[categorical_rng(weight)];
  }
}
