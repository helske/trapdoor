data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] z;
  vector[n] w;
  int<lower=1> Kx;
  int<lower=1> Kz;
  vector[Kx] do_x;
  vector[Kz] fixed_z;
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
}

model {
  w ~ normal(a_w, s_w);
  x ~ normal(a_x + b_xz * z + b_xw * w, s_x);
  y ~ normal(a_y + b_yx * x + b_yz * z + b_yw * w, s_y);
}


generated quantities {
  matrix[Kx, Kz] mean_ydox;
  real d = (b_yw * b_xw * s_w^2)/(b_xw^2 * s_w^2 + s_x^2);
  real ay = a_y + (b_yw * s_x^2 * a_w)/(b_xw^2 * s_w^2 + s_x^2) - d * a_x;
  real bx = b_yx + d;
  real bz = b_yz - b_xz * d;
  for(i in 1:Kz) {
    for(j in 1:Kx) {
      mean_ydox[j, i] = ay + bx * do_x[j] + bz * fixed_z[i];
    }
  }
}
