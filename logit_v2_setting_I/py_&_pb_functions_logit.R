##################################
##### MARGINAL PROBABILITIES #####
##################################

##### marginal p_y #####
py_fun <- function(beta, gamma, pz = 0.5) {
 py <-  (1 - pz) * ( plogis( beta[1] )           * ( 1 - plogis(gamma[1]) ) +
                     plogis( beta[1] + beta[2] ) * plogis(gamma[1]) ) +
        pz      * (  plogis( beta[1] + beta[3] ) * ( 1 - plogis(gamma[1] + gamma[2]) ) +
                     plogis( beta[1] + beta[2] + beta[3] + beta[4] ) * plogis(gamma[1] + gamma[2]) )
 return(py)
}

py_fun(beta  = c(1, 1, 1, 1), 
       gamma = c(-0.83, 3))

###### marginal Pb ######
pb_fun <- function(beta, gamma, psi, pz) {
  pb <- (1 - pz) *  ( ( ( plogis( beta[1] + beta[2] ) - plogis( beta[1] + beta[2] - psi[2]) ) * plogis( gamma[1] ) +
                    ( plogis( beta[1] + psi[1] ) - plogis( beta[1] ) ) * (1 - plogis( gamma[1] ) ) ) ) +
        pz       *    ( ( plogis( beta[1] + beta[2] + beta[3] + beta[4] ) - plogis( beta[1] + beta[2] + beta[3] + beta[4] - psi[2]) ) * plogis( gamma[1] + gamma[2] ) +
                    ( plogis( beta[1] + beta[3] + psi[1] ) - plogis( beta[1] + beta[3] ) ) * ( 1 - plogis( gamma[1] + gamma[2] ) ) )
  return(pb) 
}

# pb <- pb_fun(
#        beta  = c(1, 1, 1, 1), 
#        gamma = c(-0.83, 3), 
#        psi   = c(1, 1.5)
# )
# 
# 
# NNT = 1 / pb
# NNT

###### marginal Pb(u) ######
pbu_fun <- function(beta, gamma, psi, pz.u = 0.5) {
  pbu <- pz.u      * ( plogis(beta[1] + beta[3] + psi[1]) - plogis(beta[1] + beta[3]) ) +
        (1 - pz.u) * ( plogis(beta[1] + psi[1]) - plogis(beta[1]) )
  return(pbu) 
}

# pbu <- pbu_fun(
#   beta  = c(1, 1, 1, 1), 
#   gamma = c(-0.83, 3), 
#   psi   = c(1, 1.5)
# )
# 
# NNE = 1 / pbu
# NNE

###### marginal Pb(e) ######
pbe_fun <- function(beta, gamma, psi, pz.e = 0.5) {
  pbe = pz.e      * ( plogis(beta[1] + beta[2] + beta[3] + beta[4]) - plogis(beta[1] + beta[2] + beta[3] + beta[4] - psi[2]) ) +
       (1 - pz.e) * ( plogis(beta[1] + beta[2]) - plogis(beta[1] + beta[2] - psi[2]) )
  return(pbe) 
}

# pbe <- pbe_fun(
#   beta  = c(1, 1, 1, 1), 
#   gamma = c(-0.83, 3), 
#   psi   = c(1, 1.5)
# )
# 
# EIN = 1 / pbe
# EIN

# for I_u psi = (psi_u, psi_e) - valid IV
au_fun =  function(beta, gamma, psi) {
 au =  ( plogis( beta[1] + beta[3] ) * ( 1 - plogis( gamma[1] + gamma[2] ) )  +
         plogis( beta[1] + beta[2] + beta[3] + beta[4] - psi[2] ) *  plogis( gamma[1] + gamma[2] ) )  -
       ( plogis( beta[1] ) * ( 1 -  plogis( gamma[1] ) )  +
         plogis( beta[1] + beta[2] - psi[2]) * plogis( gamma[1] )  )
 return(au)
}

# au <- au_fun(beta   = c(1, 1, 1, 1), 
#              gamma  = c(-0.83, 3), 
#              psi    = c(1, 1.5))
# 
# au

# for I_e - valid IV

ae_fun <- function(beta, gamma, psi) {
ae =  ( plogis( beta[1] + beta[3] + psi[1] ) * ( 1 - plogis( gamma[1] + gamma[2] ) )  +
        plogis( beta[1] + beta[2] + beta[3] + beta[4] ) * plogis( gamma[1] + gamma[2] ) )  -
      ( plogis( beta[1] + psi[1] ) * ( 1 -  plogis( gamma[1] ) )  +
        plogis( beta[1] + beta[2] ) * plogis( gamma[1] )  )
 return(ae)
}

# ae <- ae_fun(beta = c(1, 1, 1, 1), 
#              gamma =  c(-0.83, 3), 
#              psi    = c(1, 1.5))
# 
# ae


##############################
##### SOLVING W.R.T BETA #####

fun_psi2 <- function(beta   = c(-1.537871, 1, 4.294861, -4.60807), 
                     gamma  =  c(-0.83, 3),       # for P(X=1) = 0.6
                     py     = 0.3,               # for P(Y=1) = 0.3
                     psi    = c(1, 1.5),          # psi = (psi_u, psi_e)
                     pb     = 1/3, 
                     pz     = 0.5
) {
  y    <- numeric()
  y[1] <- au_fun(beta, gamma, psi)
  y[2] <- ae_fun(beta, gamma, psi)
  y[3] <- py - py_fun(beta, gamma, pz)
  y[4] <- pb - pb_fun(beta, gamma, psi, pz)
  
  return(y)
}

# fun_psi2()

# psi = c(0.5, 1.5), NNT = 3
fun_x2 = function(x) { fun_psi2(beta = c(x[1], x[2], x[3], x[4]), 
                                gamma =  c(-0.83, 3),       # for P(X=1) = 0.6
                                py     = 0.3,               # for P(Y=1) = 0.3
                                psi    = c(1, 1.5),         # psi = (psi_u, psi_e)
                                pb     = 1/3, 
                                pz     = 0.5
) }

# library(nleqslv)

nleqslv(x  = c(0, 0, 0, 0), 
        fn = fun_x2)

betNNT3 <- nleqslv(x  = c(0, 0, 0, 0), 
                   fn = fun_x2)$x
betNNT3

# psi = c(1, 1.5), NNT = 5
# PROBLEMS WITH PSI_E CONVERGENCE
fun_x2 = function(x) { fun_psi2(beta = c(x[1], x[2], x[3], x[4]), 
                                gamma  = c(-0.83, 3),       # for P(X=1) = 0.6
                                py     = 0.3,               # for P(Y=1) = 0.3
                                psi    = c(1, 1.5),         # psi = (psi_u, psi_e)
                                pb     = 1/4.63, 
                                pz     = 0.5
) }
