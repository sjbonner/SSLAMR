model{
  ## Likelihood
  for (i in 1:n1) {
    mu[index1[i]] <- beta0 * width[index1[i]] + 
      inprod(X[i,slots[i,1:nslots[i]]], beta[slots[i,1:nslots[i]]])
  }
 
  for (i in 1:n2) {
    mu[index2[i]] <- beta0 * width[index2[i]]
  }
  
  for(i in 1:(n1 + n2)){
    y[i] ~ dpois(mu[i])
  }
  
  ## Prior distributions
  for (k in 1:K) {
    beta_tmp[k] ~ dt(0, 1/beta_tmp_sd^2, beta_tmp_k)T(0,)
    gamma[k] ~ dbern(gamma_p)
    beta[k] <- beta_tmp[k] * gamma[k]
  }
  
  beta0 ~ dnorm(beta0_mu,1/beta0_sd^2)T(0,)
}