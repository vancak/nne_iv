#################################
##### SANDWICH MATRIX LOGIT #####
#################################
library(pracma)

y = 1
a = 1 
z = 1

qvec <- function(
  #y,a, z,
  b0, b1, b2, b3,
  psi0, psi1,
  pz,
  p0, p1, pb,
  g0, g1, gb) {
  out <- c(
    (y - plogis(b0 + b1*a + b2*z + b3*z*a)),
    (y - plogis(b0 + b1*a + b2*z + b3*z*a))*a,
    (y - plogis(b0 + b1*a + b2*z + b3*z*a))*z, 
    (y - plogis(b0 + b1*a + b2*z + b3*z*a))*z*a,
    (z - pz) * plogis(b0 + b1*a + b2*z + b3*z*a + psi0*(1 - a)), 
    (z - pz) * plogis(b0 + b1*a + b2*z + b3*z*a - psi1*a),
    z - pz,
    ( plogis(b0 + b1*a + b2*z + b3*z*a + psi0*(1 - a)) - plogis(b0 + b1*a + b2*z + b3*z*a) )  - p0 * (1 - a),
    ( plogis(b0 + b1*a + b2*z + b3*z*a) - plogis(b0 + b1*a + b2*z + b3*z*a - psi1*a) ) - p1 * a,
    ( plogis(b0 + b1*a + b2*z + b3*a*z + psi0*(1 - a)) - plogis(b0 + b1*a + b2*z + b3*a*z - psi1*a) ) - pb,
    1/p0 - g0,
    1/p1 - g1,
    1/pb - gb) 

  return(out)
}


qvec2 <- function(x) {
  out <- qvec(#y = 1, a = 1, z = 1, 
    b0   = x[1], b1   = x[2], b2 = x[3], b3 = x[4],
    psi0 = x[5], psi1 = x[6],
    pz   = x[7],
    p0   = x[8], p1  = x[9], pb = x[10],
    g0   = x[11], g1  = x[12], gb = x[13])
  return(out)
}
