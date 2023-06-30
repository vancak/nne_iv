#######################
##### SIMULATIONS #####
#######################

library(nleqslv)
source("psi_functions_probit.R")
source("dh_functions_probit.R")
source("py_&_pb_functions_probit.R")
source("sandwich_matrix_probit_clean_v.R")

##### ASSOCIATION MODEL COEFFICIENTS #####
nleqslv(x  = c(0, 0, 0, 0), 
        fn = fun_x2)

betNNT3 <- nleqslv(x  = c(0, 0, 0, 0), 
                   fn = fun_x2)$x
betNNT3

### NNT ###
pb <- pb_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5), 
  pz    = 0.5
)

NNT = 1 / pb
NNT

### NNE ###
gamma = c(-0.83, 3)
pz_u  <- (1 - plogis(gamma[1] + gamma[2]))/2 * ( (1 - plogis(gamma[1] + gamma[2]))/2
                                                +
                                                 (1 - plogis(gamma[1]))/2 ) ^ (-1)

pbu <- pbu_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5), 
  pz.u  = pz_u
)


NNE = 1 / pbu
NNE

### EIN ###
pz_e  <- plogis(gamma[1] + gamma[2])/2 * ( plogis(gamma[1] + gamma[2])/2
                                                 +
                                             plogis(gamma[1])/2 ) ^ (-1)

pbe <- pbe_fun(
  beta  = betNNT3, 
  gamma = c(-0.83, 3), 
  psi   = c(1, 1.5), 
  pz.e  = pz_e
)

EIN = 1 / pbe
EIN

#######################
##### SIMULATIONS #####
#######################
n = c(2000, 4000, 8000, 16000)
#n = 4000
m = 1000
# simulation 

set.seed(123)

CI_NNT  <- matrix(NA, ncol = 2, nrow = m)
CI_EIN  <- matrix(NA, ncol = 2, nrow = m)
CI_NNE  <- matrix(NA, ncol = 2, nrow = m)

NNT_mat <- matrix(NA, ncol = 2, nrow = m)
EIN_mat <- matrix(NA, ncol = 2, nrow = m)
NNE_mat <- matrix(NA, ncol = 2, nrow = m)

psi_mat <- matrix(NA, ncol = 2, nrow = m)

NNT_psi   <- data.frame(matrix(NA, ncol = length(n), nrow = m))
NNT_undj  <- data.frame(matrix(NA, ncol = length(n), nrow = m))
EIN_psi   <- data.frame(matrix(NA, ncol = length(n), nrow = m))
EIN_undj  <- data.frame(matrix(NA, ncol = length(n), nrow = m))
NNE_psi   <- data.frame(matrix(NA, ncol = length(n), nrow = m))
NNE_undj  <- data.frame(matrix(NA, ncol = length(n), nrow = m))

NNT_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))
EIN_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))
NNE_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))

names(NNT_psi)  <- n
names(NNT_undj) <- n
names(EIN_psi)  <- n
names(EIN_undj) <- n
names(NNE_psi)  <- n
names(NNE_undj) <- n


k = 1

for (j in n) {
for (i in 1:m) {
 
  # set the marginal P(Z=1)
  Z    <- rbinom(j, 1, 0.5)
  
  # mean(z)
  gamz <- 3
  
  # set the marginal P(X=1) = p_x11
  # gam0 = fun_gam0(p_x11 = 0.6, gamzz = gamz )
  
  A    <- rbinom(j, 1, prob = plogis( -0.83 + 3 * Z ) )  
  
#  coef_gam <- coef(glm(A ~ Z, family = "binomial"))
  
    beta = betNNT3                                                ### NNT = 3
  # beta = c(-1.20761098, -0.07916555,  3.95015511, -3.46592022)  ### NNT = 5
  # beta <- c(-1.983679,  1.138180,  3.332185, -2.965986)         ### NNT = 2.8
  
  Y = rbinom(j, 1, prob = pnorm(beta[1] + beta[2]*A + beta[3]*Z + beta[4]*I(A*Z) ))
   
  coef_bet <- coef(glm(Y ~ A + Z + I(A*Z), family = binomial(link = "probit")))
  
  psi_e    <- nleqslv(x = 0, fn = dh_u_beta, beta = coef_bet)$x
  
  psi_u    <- nleqslv(x = 0, fn = dh_e_beta, beta = coef_bet)$x

  coef_psi <- c(psi_u, psi_e)
  
  # check 
  #
  # nleqslv(x = 0, fn = dh_u_beta, beta = coef_bet)
  # nleqslv(x = 0, fn = dh_e_beta, beta = coef_bet)
  
  # psi_mat[i,] <- c(psi_u, ifelse(abs(psi_e) < 10, psi_e, 2))
  
  p_z <- mean(Z)
  
  psi_mat[i,] <- c(psi_u, psi_e)
  
  pb <- mean( pnorm(coef_bet[1] + coef_bet[2]*A + coef_bet[3]*Z + coef_bet[4]*A*Z + coef_psi[1]*(1 - A)) - 
                pnorm(coef_bet[1] + coef_bet[2]*A + coef_bet[3]*Z + coef_bet[4]*A*Z - coef_psi[2]*A) )
  
  NNT_mat[i,] = c( 1 / (mean(Y[A == 1]) - mean(Y[A == 0])), 
                   1 / pb )
                     
  # pb_fun(beta = coef_bet, gamma = coef_gam, psi = c(psi_u, psi_e), pz = p_z) 

  # coef_red   <- coef(glm(Y ~ A, family = "binomial"))
  
  p_ein1     <- sum( A * pnorm(coef_bet[1] + coef_bet[2] * 1 + coef_bet[3] * Z + coef_bet[4] * I(1*Z)) ) / sum(A)
  p_ein0     <- sum( A * pnorm(coef_bet[1] + coef_bet[2] * 0 + coef_bet[3] * Z + coef_bet[4] * 0) ) / sum(A)
  
  p_ze       <- mean(Z[A == 1])
  
  pbe         <- pbe_fun(beta = coef_bet, gamma = coef_gam, psi = c(psi_u, psi_e), pz.e = p_ze) 

  EIN_mat[i,] <- c( 1 / (p_ein1 - p_ein0),
                    1 / pbe )

  p_nne0     <- sum( (1 - A) * pnorm(coef_bet[1] + coef_bet[2] * 0 + coef_bet[3] * Z + coef_bet[4] * 0)) / sum(1 - A)
  p_nne1     <- sum( (1 - A) * pnorm(coef_bet[1] + coef_bet[2] * 1 + coef_bet[3] * Z + coef_bet[4] * Z * 1)) / sum(1 - A)
  
  p_zu       <- mean(Z[A == 0])
  
  pbu        <- pbu_fun(beta = coef_bet, gamma = coef_gam, psi = c(psi_u, psi_e), pz.u = p_zu)

  NNE_mat[i,] = c( 1 / (p_nne1 - p_nne0), 
                   1 / pbu )
  ####################################  
  ##### CIs SANDWICH NNT NNE EIN #####
  ####################################
  
  ##### THE A BREAD MATRIX #####
  bread_mat <- matrix(0, 13, 13)
  
#  pb_new <- pb  
  
  for (b in 1:length(Y)) {
    y <-  Y[b]
    a <-  A[b]
    z <-  Z[b]
    
    jac_mat <- -jacobian(f = qvec2, x0 = c(coef_bet,
                                           coef_psi,
                                           p_z,
                                           pbu, pbe, pb, 
                                           1/pbu, 1/pbe, 1/pb))
    bread_mat <- bread_mat + jac_mat
  }
  
  bread_mat <- 1/length(Y) * bread_mat
  #bread_mat
  
  if (kappa(bread_mat) > 10^15) { 
    
    print("Singular") 
    
    i <- i + 1
    
  }  else {
  inv_a <- solve(bread_mat)
  
  ##### THE B MEAT MATRIX #####
  meat_mat <- matrix(0, 13, 13)
  
  for (b in 1:length(Y)) {
    y <- Y[b]
    a <- A[b]
    z <- Z[b]
    
    out_prod <-  outer(qvec2(x = c(coef_bet,
                                   coef_psi,
                                   p_z,
                                   pbu, pbe, pb, 
                                   1/pbu, 1/pbe, 1/pb)), 
                       qvec2(x = c(coef_bet,
                                   coef_psi,
                                   p_z,
                                   pbu, pbe, pb, 
                                   1/pbu, 1/pbe, 1/pb)) )
    meat_mat <- meat_mat + out_prod
  }
  
  
  meat_mat <- 1/length(Y) * meat_mat
  # meat_mat
  
  ##### THE SANDWICH inv(A) %*% B %*% t(inv(A)) MATRIX ##### 
  
  sand_mat <- 1 / length(Y) * inv_a %*% meat_mat %*% t(inv_a) 
  
  CI_NNE[i,1] <- NNE_mat[i,2] - 1.96 * sqrt(sand_mat[11, 11]); CI_NNE[i,2]  <- NNE_mat[i,2] + 1.96 * sqrt(sand_mat[11, 11])
  CI_EIN[i,1] <- EIN_mat[i,2] - 1.96 * sqrt(sand_mat[12, 12]); CI_EIN[i,2]  <- EIN_mat[i,2] + 1.96 * sqrt(sand_mat[12, 12])
  CI_NNT[i,1] <- NNT_mat[i,2] - 1.96 * sqrt(sand_mat[13, 13]); CI_NNT[i,2]  <- NNT_mat[i,2] + 1.96 * sqrt(sand_mat[13, 13]) 
  
  print(i)
}
}
  
  NNT_psi[,k]   <- NNT_mat[,2]
  NNT_undj[,k]  <- NNT_mat[,1]
  EIN_psi[,k]   <- EIN_mat[,2]
  EIN_undj[,k]  <- EIN_mat[,1]
  NNE_psi[,k]   <- NNE_mat[,2]
  NNE_undj[,k]  <- NNE_mat[,1]
  
  NNT_psi_CI[,(2*k - 1)]  <- CI_NNT[,1]; NNT_psi_CI[,2*k]  <- CI_NNT[,2]
  NNE_psi_CI[,(2*k - 1)]  <- CI_NNE[,1]; NNE_psi_CI[,2*k]  <- CI_NNE[,2]
  EIN_psi_CI[,(2*k - 1)]  <- CI_EIN[,1]; EIN_psi_CI[,2*k]  <- CI_EIN[,2]
  
  ##### forwarding the loop  
  k = k + 1
  
  print(k)
}

COVER_NNE <- c()
COVER_EIN <- c()
COVER_NNT <- c()

for (k in c(1, 3, 5, 7)) {
#  k = 1
  COVER_NNE[(k + 1)/2] <- mean(ifelse(NNE_psi_CI[,k] < NNE & NNE_psi_CI[,k + 1] > NNE, 1, 0), na.rm = T)
  COVER_EIN[(k + 1)/2] <- mean(ifelse(EIN_psi_CI[,k] < EIN & EIN_psi_CI[,k + 1] > EIN, 1, 0), na.rm = T)
  COVER_NNT[(k + 1)/2] <- mean(ifelse(NNT_psi_CI[,k] < NNT & NNT_psi_CI[,k + 1] > NNT, 1, 0), na.rm = T)
}

colMeans(psi_mat)
EIN
colMeans(EIN_psi)
colMeans(EIN_undj)
NNE
colMeans(NNE_psi)
colMeans(NNE_undj)
NNT
colMeans(NNT_psi)
colMeans(NNT_undj)

apply(psi_mat,  2, sd)
apply(NNT_mat,  2, sd)
apply(NNE_psi,  2, sd)
apply(EIN_psi,  2, sd)

write.csv(NNE_psi_CI, "NNE_psi_CI_PROBIT_m1000.csv", row.names = F)
write.csv(EIN_psi_CI, "EIN_psi_CI_PROBIT_m1000.csv", row.names = F)
write.csv(NNT_psi_CI, "NNT_psi_CI_PROBIT_m1000.csv", row.names = F)

