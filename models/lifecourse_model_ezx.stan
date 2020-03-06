// Model for FSD data, model Z with X

functions {
  
  // monotonic effect and sratio pmf copied and modified from brms R package
  
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
  // W, sosioeconomic status of parents from 0 to 2 (highest to lowest)
  int<lower=0, upper=2> ses[n];       
  // I, ITPA scores
  vector<lower=0>[n] itpa;
  // G, 0 = Male, 1 = Female
  int<lower=-1, upper=1> gender[n];   
  // Number of Monte Carlo samples per I and G (outer loop)
  int<lower=1> N; 
  // Number of Monte Carlo samples per W (inner loop)
  int<lower=1> M;
}


parameters {
  
  // sequential ordinal model for W and X
  // low sosioeconomic status as reference
  vector[2] a_w;
  
  // secondary level as reference
  vector[2] a_x; 
  real b_xg;
  real b_xi;
  real b_xz;
  real b_xw;
  simplex[2] simplex_xw;
  
  // Beta distribution for Z | X, I, G
  real b_zg_phi;
  vector[3] b_zx_phi;
  real a_z;
  real b_zg;
  real b_zi;
  real b_zx;
  simplex[2] simplex_zx;
  
  // ITPA conditioned on W
  real a_iw;
  real b_iw;
  simplex[2] simplex_iw;
  real<lower=0> sigma_iw;
  
  // ITPA unconditional distribution
  real a_i;
  real<lower=0> sigma_i;
  
  // Gender as Bernoulli
  real<lower=0,upper=1> p_gender;
  
  // Income as Gamma
  real a_y; 
  real b_yg;
  real b_yi;
  real b_yz;
  real b_yw;
  simplex[2] simplex_yw;
  real b_yx;
  simplex[2] simplex_yx;
  real<lower=0> shape_y;
  
  
}
model {
  
  for(i in 1:n) {
    
    real mu_education = b_xg * gender[i] + b_xi * itpa[i] + b_xz * grade[i] + 
    b_xw * mo(simplex_xw, ses[i]);
    
    real mu_income = a_y + b_yg * gender[i] + b_yi * itpa[i] + b_yz * grade[i] + 
    b_yx * mo(simplex_yx, education[i]) + b_yw * mo(simplex_yw, ses[i]);
    
       real mu_grade = inv_logit(a_z + b_zg * gender[i] + b_zi * itpa[i] + 
    b_zx * mo(simplex_zx, education[i]));
    
    real phi_grade = exp(b_zg_phi * gender[i] + b_zx_phi[education[i]+1]);
    
    ses[i] ~ sratio_logit(0, a_w);
    
    education[i] ~ sratio_logit(mu_education, a_x);
    
    income[i] ~ gamma(shape_y, shape_y * exp(-mu_income));
    
    grade[i] ~ beta(mu_grade * phi_grade, (1 - mu_grade) * phi_grade);
    
    itpa[i] ~ normal(a_i, sigma_i);
    itpa[i] ~ normal(a_iw + b_iw * mo(simplex_iw, ses[i]), sigma_iw);
    
    gender[i] ~ bernoulli(p_gender);
    
  }
  
  shape_y ~ gamma(0.01, 0.01);
}


generated quantities {
  
  // causal part
  
  // expected value of Y|do(X=x)
  vector[3] eydox;
  vector[3] eydox_mcse;
  
  // sample from Y|do(X=x)
  vector[3] ydox_sample;
  
  
  // each do_x value
  for(k in 1:3) {
    vector[M * N] sample;
    vector[M * N] weight;
    for(i in 1:N) {
      
      // sample gender and ITPA
      real itpa_i = normal_rng(a_i, sigma_i);
      int gender_i = bernoulli_rng(p_gender);
      real trapdoor_grade = 
        inv_logit(a_z + b_zg * gender_i + b_zi * itpa_i + b_zx * mo(simplex_zx, k - 1));

      
      vector[M] weight_jk;
      vector[M] sample_jk;
      
      for(j in 1:M) {
        // sample sosioeconomic status of parents
        int ses_i = sratio_logit_rng(0, a_w);
        
        real mu_education = b_xg * gender_i + b_xi * itpa_i + b_xz * trapdoor_grade +
        b_xw * mo(simplex_xw, ses_i);
        
        real mu_income = a_y + b_yg * gender_i + b_yi * itpa_i + b_yz * trapdoor_grade +
        b_yx * mo(simplex_yx, k - 1) + b_yw * mo(simplex_yw, ses_i);
        
        sample_jk[j] = gamma_rng(shape_y, shape_y * exp(-mu_income));
        weight_jk[j] = sratio_logit_lpmf(k - 1 | mu_education, a_x) +
          normal_lpdf(itpa_i | a_iw + b_iw * mo(simplex_iw, ses_i), sigma_iw);
      }
      
      // normalize the weights
      weight_jk = exp(weight_jk - max(weight_jk));
      // scale with N so that the total sum is 1
      weight_jk =  weight_jk / sum(weight_jk) / N; 
      
      weight[((i-1)*M + 1):(i*M)] = weight_jk;
      sample[((i-1)*M + 1):(i*M)] = sample_jk;
    }
    // compute expected value and sample one replication
    
    eydox[k] = dot_product(sample, weight);
    eydox_mcse[k] = sqrt(sum(square(weight .* (sample - eydox[k]))));
    ydox_sample[k] = sample[categorical_rng(weight)];
  }
}
