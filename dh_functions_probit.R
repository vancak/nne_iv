##############################################
### D * h functions - probit causal model  ###
##############################################

dh_u_probit = function(psi_e, beta) { 
              mu_zz = mean(Z)
              sum( ( Z - mu_zz ) *   pnorm(beta[1]    +   beta[2]*A + 
                                           beta[3]*Z  +   beta[4]*A*Z - 
                                           psi_e * A) )  
} 


dh_e_probit = function(psi_u, beta) { 
              mu_zz = mean(Z)
              sum( ( Z - mu_zz ) *   pnorm(beta[1]    +   beta[2]*A   + 
                                           beta[3]*Z +   beta[4]*A*Z + 
                                           psi_u * (1 - A)) )  
} 

# fixing beta vector, solving for psi exposed 
dh_u_beta = function(x, beta) {
  dh_u_probit(psi_e = x, beta)
} 

# fixing  beta vector, solving for psi unexposed
dh_e_beta = function(x, beta) {
  dh_e_probit(psi_u = x, beta)
} 

