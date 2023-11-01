model{
  ## Likelihood
  for (i in 1:n) {
    mu[i] <- beta0 * width[i] + 
      inprod(X[i,slots[i,1:nslots[i]]], beta[slots[i,1:nslots[i]]])
    y[i] ~ dpois(mu[i])
  }
  
  ## Prior distributions
  for (k in 1:K) {
    beta_tmp[k] ~ dnorm(beta_tmp_mu,1/beta_tmp_sd^2)T(0,)
    gamma[k] ~ dbern(gamma_p)
    beta[k] <- beta_tmp[k] * gamma[k]
  }
  
  beta0 ~ dnorm(beta0_mu,1/beta0_sd^2)T(0,)
  
  beta_tmp_sd ~ dt(0, beta_tmp_sd_tau, beta_tmp_sd_k)T(0,)
}