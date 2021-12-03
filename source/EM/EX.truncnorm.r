# Helper functions for the E step
EX.truncnorm = function(mu, sig, a, b){
  a.ratio = (a-mu)/sig
  b.ratio = (b-mu)/sig
  num = (dnorm(a.ratio) - dnorm(b.ratio))*sig
  den = pnorm(b.ratio) - pnorm(a.ratio)
  EX = mu + num / den
  num2 = ( a.ratio * dnorm(a.ratio) - b.ratio * dnorm(b.ratio))
  EX2 = (1 + num2 / den) * sig^2 + 2 * mu * EX - mu^2
  cbind(EX, EX2)
}

EX.truncnorm.rightcens = function(mu, sig, a, toler = 1e-10){
  a.ratio = (a-mu)/sig
  Q = ifelse(pnorm(a.ratio) == 1, toler, 1 - pnorm(a.ratio) )
  dpnorm.a.ratio = dnorm(a.ratio) / Q
  EX = mu + sig * dpnorm.a.ratio
  EX2 = mu^2 + sig^2 + sig * (a + mu) * dpnorm.a.ratio
  cbind(EX, EX2)
}

EX.norm = function(mu, sig){
  EX = mu
  EX2 = sig^2 + EX^2
  cbind(EX, EX2)
}