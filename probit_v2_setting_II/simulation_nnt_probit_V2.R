#######################
##### SIMULATIONS #####
#######################
library(nleqslv)
library(pracma)
source("dh_functions_probit_V2.R")
source("py_&_pb_functions_probit_V2.R")
source("sandwich_matrix_probit_clean_v_V2.R")

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
  psi   = c(0.5, 1), 
  pz    = 0.7
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
  psi   = c(0.5, 1), 
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
  psi   = c(0.5, 1), 
  pz.e  = pz_e
)

EIN = 1 / pbe
EIN

#######################
##### SIMULATIONS #####
#######################
n = c(500, 1000, 2000, 4000)
#n = 4000
m = 1000
# simulation 

set.seed(1984)

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
IV_str    <- data.frame(matrix(NA, ncol = length(n), nrow = m))

NNT_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))
EIN_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))
NNE_psi_CI   <- data.frame(matrix(NA, ncol = length(n) * 2, nrow = m))

names(NNT_psi)  <- n
names(NNT_undj) <- n
names(EIN_psi)  <- n
names(EIN_undj) <- n
names(NNE_psi)  <- n
names(NNE_undj) <- n
names(IV_str)   <- n


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
  IV_str[i,k]    <- summary(glm(A ~ Z, family = "binomial"))$coeff[,"z value"][2]
  
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
  
  if ( kappa(bread_mat) < 10 ^ 12 ) {

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
  } else { 
    CI_NNE[i,1] <- NA; CI_NNE[i,2]  <- NA
    CI_EIN[i,1] <- NA; CI_EIN[i,2]  <- NA
    CI_NNT[i,1] <- NA; CI_NNT[i,2]  <- NA 
    
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

write.csv(psi_mat, "psi_mat.csv", row.names = F)

EIN
write.csv(EIN_psi, "EIN_psi.csv", row.names = F)
write.csv(EIN_undj, "EIN_undj.csv", row.names = F)

apply(EIN_psi[, -ncol(EIN_psi)], 2, FUN = function(x) {sum(x < 1 | x > 1000, na.rm = T)})

NNE
write.csv(NNE_psi, "NNE_psi.csv", row.names = F)
write.csv(NNE_undj, "NNE_undj.csv", row.names = F)

apply(NNE_psi[, -ncol(NNE_psi)], 2, FUN = function(x) {sum(x < 1 | x > 1000, na.rm = T)})

NNT
write.csv(NNT_psi, "NNT_psi.csv", row.names = F)
write.csv(NNT_undj, "NNT_undj.csv", row.names = F)

apply(NNT_psi[, -ncol(NNT_psi)], 2, FUN = function(x) {sum(x < 1 | x > 1000, na.rm = T)})


colMeans(psi_mat, na.rm = T)
EIN
colMeans(EIN_psi[EIN_psi > 1,-ncol(EIN_psi)], na.rm = T)
colMeans(EIN_undj)
NNE
colMeans(NNE_psi[NNE_psi > 1,-ncol(NNE_psi)], na.rm = T)
colMeans(NNE_undj)
NNT
colMeans(NNT_psi[NNT_psi > 1,-ncol(NNT_psi)], na.rm = T)
colMeans(NNT_undj)

# apply(psi_mat,  2, sd)
# sd_NNT <- apply(NNT_mat,  2, sd)

### SD vectors 
sd_NNT <- apply(NNT_psi[NNT_psi < 1000 & NNT_psi > 0,-ncol(NNT_psi)],  2, sd, na.rm = T)
sd_NNE <- apply(NNE_psi[NNE_psi < 1000 & NNE_psi > 0,-ncol(NNE_psi)],  2, sd, na.rm = T)
sd_EIN <- apply(EIN_psi[EIN_psi < 1000 & EIN_psi > 0,-ncol(EIN_psi)],  2, sd, na.rm = T)

sd_mat <- data.frame(rbind(sd_NNT, sd_NNE, sd_EIN))

### AVERAGE BIAS
av_bias <- function(x, y) {
  mean( abs( x - y ), na.rm = T) 
}

avb_NNT <- apply(NNT_psi[NNT_psi < 1000 & NNT_psi > 0,-ncol(NNT_psi)], 2, FUN = av_bias, y = NNT)
avb_NNE <- apply(NNE_psi[NNE_psi < 1000 & NNE_psi > 0,-ncol(NNE_psi)], 2, FUN = av_bias, y = NNE)
avb_EIN <- apply(EIN_psi[EIN_psi < 1000 & EIN_psi > 0,-ncol(EIN_psi)], 2, FUN = av_bias, y = EIN)

avb_mat <- data.frame(rbind(avb_NNT, avb_NNE, avb_EIN))

### count irelevnat point estimators
nrow( NNT_psi[NNT_psi > 1000 | NNT_psi < 0, -ncol(NNT_psi)])
nrow( NNE_psi[NNE_psi > 1000 | NNE_psi < 0, -ncol(NNT_psi)])
nrow( EIN_psi[EIN_psi > 1000 | EIN_psi < 0, -ncol(NNT_psi)])

colMeans(IV_str)

write.csv(IV_str, "IV_strength.csv", row.names = T)

### CSV FILES ###
## SE
write.csv(sd_mat, "sd_PROBIT_m1000.csv", row.names = T)

## AVERAGE BIAS
write.csv(avb_mat, "avb_PROBIT_m1000.csv", row.names = T)

## CIs
write.csv(NNE_psi_CI, "NNE_psi_CI_PROBIT_m1000.csv", row.names = F)
write.csv(EIN_psi_CI, "EIN_psi_CI_PROBIT_m1000.csv", row.names = F)
write.csv(NNT_psi_CI, "NNT_psi_CI_PROBIT_m1000.csv", row.names = F)

### COVERAGE MATRIX - CREATION ###
NNE_psi_CI <- read.csv("NNE_psi_CI_PROBIT_m1000.csv")
EIN_psi_CI <- read.csv("EIN_psi_CI_PROBIT_m1000.csv")
NNT_psi_CI <- read.csv("NNT_psi_CI_PROBIT_m1000.csv")

# COUNT NON-INVERTIBLE COV MATRICES
NNE_sing <- apply(NNE_psi_CI[, c(1, 3, 5, 7)], 2, FUN = function(x) {sum(is.na(x))})
EIN_sing <- apply(EIN_psi_CI[, c(1, 3, 5, 7)], 2, FUN = function(x) {sum(is.na(x))})
NNT_sing <- apply(NNT_psi_CI[, c(1, 3, 5, 7)], 2, FUN = function(x) {sum(is.na(x))})

write.csv(data.frame( rbind(NNE_sing,
                            EIN_sing,
                            NNT_sing)), "singular_matrices.csv", row.names = F)

# COUNT PRACTICALLY INFINITE CIs
NNT_inf_ci <- apply(NNT_psi_CI[, c(2, 4, 6, 8)], 2, FUN = function(x) {sum(x > 10^3 , na.rm = T)})
EIN_inf_ci <- apply(EIN_psi_CI[, c(2, 4, 6, 8)], 2, FUN = function(x) {sum(x > 10^3 , na.rm = T)})
NNE_inf_ci <- apply(NNE_psi_CI[, c(2, 4, 6, 8)], 2, FUN = function(x) {sum(x > 10^3 , na.rm = T)})

write.csv(data.frame( rbind(NNE_inf_ci,
                            EIN_inf_ci,
                            NNT_inf_ci)), "infinite_cis_by_measure.csv", row.names = F)


t(data.frame( rbind(EIN_inf_ci,
                    NNE_inf_ci,
                    NNT_inf_ci)))

COVER_NNE <- c()
COVER_EIN <- c()
COVER_NNT <- c()

for (k in c(1, 3, 5, 7)) {
  #  k=1
  COVER_NNE[(k + 1)/2] <- mean(ifelse(NNE_psi_CI[,k] < NNE & NNE_psi_CI[,k + 1] > NNE, 1, 0), na.rm = T)
  COVER_EIN[(k + 1)/2] <- mean(ifelse(EIN_psi_CI[,k] < EIN & EIN_psi_CI[,k + 1] > EIN, 1, 0), na.rm = T)
  COVER_NNT[(k + 1)/2] <- mean(ifelse(NNT_psi_CI[,k] < NNT & NNT_psi_CI[,k + 1] > NNT, 1, 0), na.rm = T)
}

### WRITING THE LONG FILE ###
### NNT ###
library(reshape2)

NNT_mat2 <- as.data.frame(NNT_mat)

colnames(NNT_mat2) <- c("NNT_UN", "NNT_IV")

NNT_mat2$id <- 1:m

NNT_psi$id <- 1:m

bbb2 <- melt(NNT_psi,
             id.var        = "id", 
             variable.name = "n", 
             value.name    = "NNT") 

NNT_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(NNT_undj,
               id.var        = "id", 
               variable.name = "n", 
               value.name    = "NNT") 

# apply(NNT_psi, 2, mean)

#bbb_un$NNT <- ifelse(bbb_un$NNT > 20, NA, bbb_un$NNT)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

# bbb_both2$NNT <- ifelse(bbb_both2$NNT > 15, NA, bbb_both2$NNT)
# bbb_both2$NNT <- ifelse(bbb_both2$NNT < -25, NA, bbb_both2$NNT)

write.csv(bbb_both2, "NNTm1000_PROBIT_5_3.csv", row.names = F)

### EIN ###
EIN_mat2 <- as.data.frame(EIN_mat)

colnames(EIN_mat2) <- c("EIN_UN", "EIN_IV")

EIN_mat2$id <- 1:m

EIN_psi$id  <- 1:m

bbb2 <- melt(EIN_psi,
             id.var        = "id",
             variable.name = "n",
             value.name    = "EIN")

EIN_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(EIN_undj,
               id.var        = "id",
               variable.name = "n",
               value.name    = "EIN")

# apply(NNT_psi, 2, mean)

# bbb_un$EIN <- ifelse(bbb_un$EIN > 20, NA, bbb_un$EIN)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

# bbb_both2$EIN <- ifelse(bbb_both2$EIN > 10, NA, bbb_both2$EIN)
# bbb_both2$EIN <- ifelse(bbb_both2$EIN < -6, NA, bbb_both2$EIN)

write.csv(bbb_both2, "EINm1000_PROBIT_4_496.csv", row.names = F)

### NNE ###

NNE_mat2 <- as.data.frame(NNE_mat)

colnames(NNE_mat2) <- c("NNE_UN", "NNE_IV")

NNE_mat2$id <- 1:m

NNE_psi$id  <- 1:m

bbb2 <- melt(NNE_psi,
             id.var        = "id", 
             variable.name = "n", 
             value.name    = "NNE") 

NNE_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(NNE_undj,
               id.var        = "id", 
               variable.name = "n", 
               value.name    = "NNE") 

# apply(NNT_psi, 2, mean)

# bbb_un$NNE <- ifelse(bbb_un$NNE > 20, NA, bbb_un$NNE)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

#write.csv(bbb_both2, "NNE_NNTbet4.csv", row.names = F)


# bbb_both2$NNE <- ifelse(bbb_both2$NNE > 10, NA, bbb_both2$NNE)
# bbb_both2$NNE <- ifelse(bbb_both2$NNE < -7, NA, bbb_both2$NNE)

write.csv(bbb_both2, "NNEm1000_PROBIT_7_248.csv", row.names = F)