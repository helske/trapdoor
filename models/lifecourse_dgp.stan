// Model for FSD data, model Z with X

functions {
  
  // sample median
  
  real median(vector x) {
    int n = num_elements(x);
    int half = n / 2;
    return sort_asc(x)[half];
  }
  // monotonic effect copied from brms R package
  
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

  
}
data {
  int<lower=1> n;
  vector[n] u1;         
  vector[n] u2;         
  vector[n] u3;         
  // Y, income
  vector<lower=0>[n] income;         
  // X, education levels from 0 to 2 (lowest to highest)
  int<lower=0, upper=2> education[n];
  // Z, grade average from elementary school scaled to [0,1]
  vector<lower=0,upper=1>[n] grade;
  // W, sosioeconomic status of parents from 1 to 3 (highest to lowest)
  int<lower=1, upper=3> ses[n];       
  // I, ITPA scores
  vector[n] itpa;
  // G, 0 = Male, 1 = Female
  int<lower=0, upper=1> gender[n];   
  // Number of Monte Carlo samples per I and G (outer loop)
  int<lower=1> N; 
  // Number of Monte Carlo samples per W (inner loop)
  int<lower=1> M;
}

parameters {
  
  // Gender as Bernoulli
  real<lower=0,upper=1> p_gender;
  
  // ITPA
  real a_i;
  real b_i3;
  real<lower=0> sigma_i;
  
  // Income as Gamma
  real a_y; 
  real b_yg;
  real b_yi;
  real b_y1;
  real b_yx;
  simplex[2] simplex_yx;
  real<lower=0> shape_y;
  
  real a_u1;
  real<lower=0> s_u1;
  real a_u3;
  real<lower=0> s_u3;
}
model {
  
  for(i in 1:n) {
    real mu_income = a_y + b_yg * gender[i] + b_yi * itpa[i] +
    b_yx * mo(simplex_yx, education[i]) + b_y1 * u1[i];
    
    income[i] ~ gamma(10000 * shape_y, 10000 * shape_y * exp(-mu_income));
  }
  itpa ~ normal(a_i + b_i3 * u3, sigma_i);
  gender ~ bernoulli(p_gender);
  u1 ~ normal(a_u1, s_u1);
  u3 ~ normal(a_u3, s_u3);
}


generated quantities {
  
  // expected value of Y|do(X=x)
  vector[3] mean_ydox;
  vector[3] median_ydox;
  // each do_x value
  for(k in 1:3) {
    vector[N] sample;
    for(i in 1:N) {
      // sample gender, ITPA, u1 and u3
      real u1_i = normal_rng(a_u1, s_u1);
      real u3_i = normal_rng(a_u3, s_u3);
      real itpa_i = normal_rng(a_i + b_i3 * u3_i, sigma_i);
      int gender_i = bernoulli_rng(p_gender);
      
      real mu_income = a_y + b_yg * gender_i + b_yi * itpa_i +
      b_yx * mo(simplex_yx, k - 1) + b_y1 * u1_i;
      
      sample[i] = gamma_rng(shape_y, shape_y * exp(-mu_income));
    }
    mean_ydox[k] = mean(sample);
    median_ydox[k] = median(sample);
  }
}
