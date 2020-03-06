data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
  vector[n] u;
  int<lower=1> Kx;
  vector[Kx] do_x;
}

parameters {
  real a_u;
  real<lower=0> s_u;

  real a_y;
  real b_yx;
  real b_yu;
  real<lower=0> s_y;
}

model {
  u ~ normal(a_u, s_u);
  y ~ normal(a_y + b_yx * x + b_yu * u, s_y);
}
generated quantities {
  vector[Kx] mean_ydox = a_y + b_yx * do_x + b_yu * a_u;
}
