#############################################
### D * h functions - logit causal model  ###
#############################################

dh_u = function(psi_e, beta) { 
  mu_zz = mean(Z)
  sum( ( Z - mu_zz ) *   plogis( beta[1]   +   beta[2]*A + 
                                 beta[3]*Z +   beta[4]*A*Z - 
                                 psi_e * A) )  
} 


dh_e = function(psi_u, beta) { 
  mu_zz = mean(Z)
  sum( ( Z - mu_zz ) *   plogis( beta[1]   +   beta[2]*A   + 
                                 beta[3]*Z +   beta[4]*A*Z + 
                                 psi_u * (1 - A)) )  
} 

# fixing beta vector, solving for psi exposed 
dh_u_beta = function(x, beta) {
  dh_u(psi_e = x, beta)
} 

# fixing  beta vector, solving for psi unexposed
dh_e_beta = function(x, beta) {
  dh_e(psi_u = x, beta)
} 
