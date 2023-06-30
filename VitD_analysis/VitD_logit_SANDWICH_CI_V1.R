###############################
##### VitD LOGIT SANDWICH #####
###############################

library(ivtools)
library(nleqslv)
source("psi_functions_logit.R")
source("dh_functions_logit.R")
source("py_&_pb_functions_logit.R")
source("sandwich_matrix_logit_clean_v.R")
# set.seed(12345)

data("VitD")

g <- function(x) {
  out <- ifelse( x > 0, 1/x, Inf)
  return(out)
}

#hist(VitD$vitd, breaks = 100)
#mean(VitD$vitd)
#sd(VitD$vitd)

VitD$expose <- ifelse(VitD$vitd >=  30, 1, 0)

summary(VitD$vitd)
mean(VitD$expose)

A <- VitD$expose

Z <- VitD$filaggrin

coef_gam <- coef(glm(expose ~ filaggrin, data = VitD, family = "binomial"))

# summary(glm(expose ~ filaggrin, data = VitD, family = "binomial"))

#######################
##### LOGIT MODEL #####
#######################

coef_bet <- coef(glm(I(1-death) ~ expose + filaggrin + I(expose*filaggrin), data = VitD, family = "binomial"))

psi_e = nleqslv(x = 1.5, fn = dh_u_beta, beta = coef_bet)$x

psi_u = nleqslv(x = 5, fn = dh_e_beta, beta = coef_bet)$x

coef_psi <- c(psi_u, psi_e)

# psi_mat[i,] <- c(psi_u, ifelse(abs(psi_e) < 10, psi_e, 2))

Y <- 1 - VitD$death

##### NNT-IV & NNT UNASJUSTED #####
pb <- pb_fun(beta   = coef_bet,
              gamma = coef_gam,
              psi   = coef_psi,
              pz    = mean(Z))

NNT <- g(pb)
NNT

pb_new <- mean(
  ( plogis(coef_bet[1] + coef_bet[2]*A + coef_bet[3]*Z + coef_bet[4]*A*Z + coef_psi[1]*(1 - A)) - 
      plogis(coef_bet[1] + coef_bet[2]*A + coef_bet[3]*Z + coef_bet[4]*A*Z - coef_psi[2]*A) ) 
  )

NNT_mat     = c( g( (mean(Y[A == 1]) - mean(Y[A == 0]))), 
                 g( pb_new ) )

NNT_mat

pb <- pb_new
##### NNE-IV #####
p_zu       <- mean(Z[A == 0])

pbu <- pbu_fun(
  beta  = coef_bet, 
#  gamma = coef_gam, 
  psi   = coef_psi, 
  pz    = p_zu
)

NNE = g(pbu)
NNE

##### NNE UNADJUSTED #####
p_nne0     <- sum( (1 - A) * plogis(coef_bet[1] + coef_bet[2] * 0 + coef_bet[3] * Z + coef_bet[4] * 0)) / sum(1 - A)
p_nne1     <- sum( (1 - A) * plogis(coef_bet[1] + coef_bet[2] * 1 + coef_bet[3] * Z + coef_bet[4] * Z * 1)) / sum(1 - A)

NNE_undj = g(p_nne1 - p_nne0) 
NNE_undj

##### EIN-IV #####
p_ze       <- mean(Z[A == 1])

pbe <- pbe_fun(
  beta  = coef_bet, 
#  gamma = coef_gam, 
  psi   = coef_psi, 
  pz.e  = p_ze
)

EIN = g(pbe)
EIN

##### UNADJUSTED EIN #####
p_ein1     <- sum( A * plogis(coef_bet[1] + coef_bet[2] * 1 + coef_bet[3] * Z + coef_bet[4] * I(1*Z))) / sum(A)
p_ein0     <- sum( A * plogis(coef_bet[1] + coef_bet[2] * 0 + coef_bet[3] * Z + coef_bet[4] * 0)) / sum(A)

EIN_undj   <-  g(p_ein1 - p_ein0)
EIN_undj


#######################
##### CI SANDWICH #####
#######################
pmz <- mean(Z)
library(pracma)

### REALITY CHECK
grad_vec <- rep(0, 13)

for (i in 1:length(Y)) {
  y <-  Y[i]
  a <-  A[i]
  z <-  Z[i]
  
  grad_tmp <- qvec2(x = c(coef_bet,
                          coef_psi,
                          pmz,
                          pbu, pbe, pb, 
                          1/pbu, 1/pbe, 1/pb))
  grad_vec <- grad_vec + grad_tmp
}

sum(grad_vec)

##### THE A BREAD MATRIX #####
bread_mat <- matrix(0, 13, 13)

for (i in 1:length(Y)) {
  y <-  Y[i]
  a <-  A[i]
  z <-  Z[i]
  
  jac_mat <- -jacobian(f = qvec2, x0 = c(coef_bet,
                                         coef_psi,
                                         pmz,
                                         pbu, pbe, pb_new, 
                                         1/pbu, 1/pbe, 1/pb_new))
  bread_mat <- bread_mat + jac_mat
}

bread_mat <- 1/length(Y) * bread_mat
#bread_mat

inv_a <- solve(bread_mat)

##### THE B MEAT MATRIX #####
meat_mat <- matrix(0, 13, 13)

for (i in 1:length(Y)) {
  y <- Y[i]
  a <- A[i]
  z <- Z[i]
  
 out_prod <-  outer(qvec2(x = c(coef_bet,
                                coef_psi,
                                pmz,
                                pbu, pbe, pb_new, 
                                1/pbu, 1/pbe, 1/pb_new)), 
                          qvec2(x = c(coef_bet,
                                      coef_psi,
                                      pmz,
                                      pbu, pbe, pb_new, 
                                      1/pbu, 1/pbe, 1/pb_new)) )
  meat_mat <- meat_mat + out_prod
}

meat_mat <- 1/length(Y) * meat_mat
#meat_mat

##### THE SANDWICH inv(A) %*% B %*% t(inv(A)) MATRIX ##### 

sand_mat <- 1 / length(Y) * inv_a %*% meat_mat %*% t(inv_a) 

NNE - 1.96 * sqrt(sand_mat[11, 11]); NNE + 1.96 * sqrt(sand_mat[11, 11])
EIN - 1.96 * sqrt(sand_mat[12, 12]); EIN + 1.96 * sqrt(sand_mat[12, 12])
NNT - 1.96 * sqrt(sand_mat[13, 13]); NNT + 1.96 * sqrt(sand_mat[13, 13])
