py_fun = function(beta, gamma, pz = 0.5) {
 py =   (1 - pz) *  ( pnorm( beta[1] )           * ( 1 - plogis( gamma[1] ) ) +
                      pnorm( beta[1] + beta[2] ) * plogis( gamma[1] ) ) +
             pz  * ( pnorm( beta[1] + beta[3] ) * ( 1 - plogis(gamma[1] + gamma[2]) ) +
                     pnorm( beta[1] + beta[2] + beta[3] + beta[4] ) * plogis( gamma[1] + gamma[2] ) )
 return(py)
}

# py_fun(beta  = c(1, 1, 1, 1), 
#        gamma = c(-0.83, 3))

###### FOR P_b ######
pb_fun <- function(beta, gamma, psi, pz = 0.5) {
  pb <- (1 - pz) *  ( ( ( pnorm( beta[1] + beta[2] ) - pnorm( beta[1] + beta[2] - psi[2]) ) * plogis( gamma[1] ) +
                          ( pnorm( beta[1] + psi[1] ) - pnorm( beta[1] ) ) * (1 - plogis( gamma[1] ) ) ) ) +
    pz       *    ( ( pnorm( beta[1] + beta[2] + beta[3] + beta[4] ) - pnorm( beta[1] + beta[2] + beta[3] + beta[4] - psi[2]) ) * plogis( gamma[1] + gamma[2] ) +
                      ( pnorm( beta[1] + beta[3] + psi[1] ) - pnorm( beta[1] + beta[3] ) ) * ( 1 - plogis( gamma[1] + gamma[2] ) ) )
  return(pb) 
}

# pb <- pb_fun(
#        beta  = betNNT4, 
#        gamma = c(-0.83, 3), 
#        psi   = c(0.5, 1.5)
#        )
# 
# 1 / pb

###### FOR P_b_u ######
pbu_fun <- function(beta, gamma, psi, pz.u = 0.5) {
  pbu <- pz.u       * ( ( pnorm(beta[1] + beta[3] + psi[1]) - pnorm(beta[1] + beta[3]) ) ) +
         (1 - pz.u) * (   pnorm(beta[1] + psi[1]) - pnorm(beta[1]) )
  
  return(pbu) 
}

pbu <- pbu_fun(
  beta  = c(1, 1, 1, 1), 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5), 
  pz.u  = 0.5
)


NNE = 1 / pbu
NNE

###### FOR P_b_e ######
pbe_fun <- function(beta, gamma, psi, pz.e = 0.5) {
  pbe <- (1 - pz.e) * ( ( pnorm(beta[1] + beta[2]) - pnorm(beta[1] + beta[2] - psi[2]) ) ) +
          pz.e      * (   pnorm(beta[1] + beta[2] + beta[3] + beta[4]) - pnorm(beta[1] + beta[2] + beta[3] + beta[4] - psi[2]) )
  return(pbe) 
}

# pbe <- pbe_fun(
#   beta  = betNNT4, 
#   gamma = c(-0.83, 3), 
#   psi   = c(1, 1.5), 
#   pz.e  = 0.5
# )
# 
# EIN = 1 / pbe
# EIN

# for I_u psi = (psi_u, psi_e)
au_fun =  function(beta, gamma, psi) {
 au =  ( pnorm( beta[1] + beta[3] ) * ( 1 - plogis( gamma[1] + gamma[2] ) )  +
               pnorm( beta[1] + beta[2] + beta[3] + beta[4] - psi[2] ) *  plogis( gamma[1] + gamma[2] ) )  -
       ( pnorm( beta[1] ) * ( 1 -  plogis( gamma[1] ) )  +
           pnorm( beta[1] + beta[2] - psi[2]) * plogis( gamma[1] )  )
 return(au)
}

# au <- au_fun(beta = betNNT4, 
#              gamma =  c(-0.83, 3), 
#              psi    = c(1, 1.5))
# 
# au

# for I_e

ae_fun <- function(beta, gamma, psi) {
ae =  ( pnorm( beta[1] + beta[3] + psi[1] ) * ( 1 - plogis( gamma[1] + gamma[2] ) )  +
             pnorm( beta[1] + beta[2] + beta[3] + beta[4] ) * plogis( gamma[1] + gamma[2] ) )  -
      ( pnorm( beta[1] + psi[1] ) * ( 1 -  plogis( gamma[1] ) )  +
             pnorm( beta[1] + beta[2] ) * plogis( gamma[1] )  )
 return(ae)
}

# ae <- ae_fun(beta = betNNT4, 
#              gamma =  c(-0.83, 3), 
#              psi    = c(1, 1.5))
# 
# ae


######
### SOLVING W.R.T BETA ###
fun_psi2 <- function(beta = c(-1.537871, 1, 4.294861, -4.60807), 
                     gamma =  c(-0.83, 3),       # for P(X=1) = 0.6
                     py     = 0.3,               # for P(Y=1) = 0.3
                     psi    = c(0.5, 1.5),          # psi = (psi_u, psi_e)
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
                                pb     = 1/3 
) }

# library(nleqslv)

nleqslv(x  = c(0, 0, 0, 0), 
        fn = fun_x2)

betNNT3 <- nleqslv(x  = c(0, 0, 0, 0), 
               fn = fun_x2)$x
betNNT3

### NNT ###
pb <- pb_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5)
)

NNT = 1 / pb
NNT

### NNE ###
pbu <- pbu_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5)
)


NNE = 1 / pbu
NNE

### EIN ###
pbe <- pbe_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5)
)

EIN = 1 / pbe
EIN

# # psi = c(1, 1.5), NNT = 5
# # PROBLEMS WITH PSI_E CONVERGENCE
# fun_x2 = function(x) { fun_psi2(beta = c(x[1], x[2], x[3], x[4]), 
#                                 gamma =  c(-0.83, 3),       # for P(X=1) = 0.6
#                                 py     = 0.3,               # for P(Y=1) = 0.3
#                                 psi    = c(1, 1.5),         # psi = (psi_u, psi_e)
#                                 pb     = 1/4 
# ) }
# 
# # library(nleqslv)
# 
# betNNT4 <- nleqslv(x = c(0, 0, 0, 0), 
#                    fn = fun_x2)$x
# betNNT4 
# 
# # psi = c(0.5, 2), NNT = 2.8
# fun_x2 = function(x) { fun_psi2(beta = c(x[1], x[2], x[3], x[4]), 
#                                 gamma  =  c(-0.83, 3),      # for P(X=1) = 0.6
#                                 py     = 0.3,               # for P(Y=1) = 0.3
#                                 psi    = c(1, 2),         # psi = (psi_u, psi_e)
#                                 pb     = 1/3 
# ) }
# 
# # library(nleqslv)
# 
# nleqslv(x = c(0, 0, 0, 0), 
#         fn = fun_x2)
# 
# 
# # psi = c(0.5, 2), NNT = 10
# # PROBLEMS WITH PSI_E CONVERGENCE
# fun_x2 = function(x) { fun_psi2(beta   = c(x[1], x[2], x[3], x[4]), 
#                                 gamma  = c(-0.83, 3),       # for P(X=1) = 0.6
#                                 py     = 0.7,               # for P(Y=1) = 0.3
#                                 psi    = c(1, 3),         # psi = (psi_u, psi_e)
#                                 pb     = 1/3
# ) }
# 
# # library(nleqslv)
# 
# nleqslv(x = c(0, 0, 0, 0), 
#         fn = fun_x2)
# 
# ##### CHECK #####
# # bb <- nleqslv(x = c(0, 0, 0, 0), 
# #               fn = fun_x2)$x
# # 
# # fun_psi2(beta = bb, 
# #          gamma =  c(-0.83, 3),       # for P(X=1) = 0.6
# #          py     = 0.3,               # for P(Y=1) = 0.3
# #          psi    = c(1, 3),         # psi = (psi_u, psi_e)
# #          pb     = 1/4 )