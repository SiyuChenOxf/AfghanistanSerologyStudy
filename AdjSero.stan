data {
  int<lower=0> N;  
  int<lower=0,upper=1> y[N]; 
  real<lower=0, upper=1> kse;
  real<lower=0, upper=1> ksp;
}
parameters {
 
      real<lower=0, upper=1> p;
}

model {
  p ~ beta(1,1);
  y ~ bernoulli(kse*p+(1-ksp)*(1-p));
}