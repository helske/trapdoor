data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] z;
  vector[n] w;
  int<lower=1> N;
  real do_x;
}

parameters {
  real a_w;
  real<lower=0> s_w;
  
  real a_x;
  real b_xz;
  real b_xw;
  real<lower=0> s_x;
  
  real a_y;
  real b_yx;
  real b_yz;
  real b_yw;
  real<lower=0> s_y;
  
  real a_z;
  real b_zx;
  real<lower=0> s_z;
}

model {
  w ~ normal(a_w, s_w);
  x ~ normal(a_x + b_xz * z + b_xw * w, s_x);
  y ~ normal(a_y + b_yx * x + b_yz * z + b_yw * w, s_y);
  z ~ normal(a_z + b_zx * x, s_z);
}


generated quantities {
  real d = (b_yw * b_xw * s_w^2)/(b_xw^2 * s_w^2 + s_x^2);
  real ay = a_y + (b_yw * s_x^2 * a_w)/(b_xw^2 * s_w^2 + s_x^2) - d * a_x;
  real bx = b_yx + d;
  real bz = b_yz - b_xz * d;
  
  real trapdoor = a_z + b_zx * do_x;
  real mean_analytical = ay + bx * do_x + bz * trapdoor;
  
  real mean_montecarlo;
  real yrep;
  real mcse;
  
  vector[N] samples;
  vector[N] weights;
  for(i in 1:N) {
    real wi = normal_rng(a_w, s_w);
    samples[i] = normal_rng(a_y + b_yx * do_x + b_yz * trapdoor + b_yw * wi, s_y);
    weights[i] = normal_lpdf(do_x | a_x + b_xz * trapdoor + b_xw * wi, s_x);
  }
  
  weights = exp(weights - max(weights));
  weights = weights / sum(weights);
  mean_montecarlo = dot_product(weights, samples);
  yrep = samples[categorical_rng(weights)];
  mcse = sqrt(sum(square(weights .* (samples - mean_montecarlo))));
  
}
