prior_simple <- function(pars = NULL){

  ## beta -- half-t
  ## Identify parameters based on quantile matching
  
  ## Set default quantiles
  if(is.null(pars$q_beta))
    q_beta <- c(5, 1000)

  if(is.null(pars$p_beta))  
    p_beta <- c(.5, .99)
  
  ## Set default max df
  if(is.null(pars$max_df))
    max_df <- 100
  
  ## Convert probabilities back to full-t
  pstar <- 1-(1-p_beta)/2
  
  ## Compute sigma over range of df
  beta_df <- tibble(df = 1:max_df) |> 
    rowwise() |> 
    mutate(sigma = sum(q_beta^2)/sum(q_beta * qt(pstar , df)),
           d = sum((q_beta - sigma * qt(pstar, df))^2)) |>
    ungroup() |> 
    filter(rank(d) == 1)
  
  ## Extract parameters
  beta_tmp_mu <- 0
  beta_tmp_sd <- pull(beta_df,"sigma")[1]
  beta_tmp_k <- pull(beta_df,"df")[1]
  
  ## gamma -- Bernoulli
  
  ## Set default
  if(is.null(pars$gamma_p))
    gamma_p <- .5
  
  ## beta0 -- half-normal
  ## Identify parameters based on quantile matching
  
  ## Set default quantiles
  q_beta0 = 10
  p_beta0 = .99
  
  ## Convert probabilities back to full-normal
  pstar <- 1 - (1- p_beta0)/2
  beta0_mu <- 0
  beta0_sd <- sum(q_beta0^2)/sum(q_beta0 * qnorm(pstar))
  
  ## Return list of prior parameters
  list(beta_tmp_mu = beta_tmp_mu,
       beta_tmp_sd = beta_tmp_sd,
       beta_tmp_k = beta_tmp_k,
       gamma_p = gamma_p,
       beta0_mu = beta0_mu,
       beta0_sd = beta0_sd)
}

prior_hierarchical <- function(pars){
  ## beta_tmp -- half-t
  
  ## Set default quantiles
  if(is.null(pars$q_beta))
    q_beta <- rbind(c(5, 1000),c(50,10000))
  
  if(is.null(pars$p_beta))  
    p_beta <- c(.5, .99)
  
  ## Set default max df
  if(is.null(pars$max_df))
    max_df <- 100
  
  ## Convert probabilities back to full-t
  pstar <- 1-(1-p_beta)/2
  
  ## Compute sigma median over range of df
  beta_df_1 <- tibble(df = 1:max_df) |> 
    rowwise() |> 
    mutate(sigma1 = sum(q_beta[1,]^2)/sum(q_beta[1,] * qt(pstar , df)),
           d1 = sum((q_beta[1,] - sigma1 * qt(pstar, df))^2))
  
  ## Compute sigma maximum over range of df
  beta_df_2 <- tibble(df = 1:max_df) |> 
    rowwise() |> 
    mutate(sigma2 = sum(q_beta[2,]^2)/sum(q_beta[2,] * qt(pstar , df)),
           d2 = sum((q_beta[2,] - sigma2 * qt(pstar, df))^2)) 
  
  ## Identify best combination of sigma and df over both
  beta_df <- full_join(beta_df_1, beta_df_2, by = "df") |> 
    mutate(d = (sqrt(d1) + sqrt(d2))^2) |> 
    ungroup() |> 
    filter(rank(d) == 1)
  
  ## Identify parameters for prior on sigma
  q_sigma <- c(beta_df$sigma1[1], beta_df$sigma2[1])
  
  sigma_df <- tibble(df = 1:max_df) |> 
    rowwise() |> 
    mutate(tau = sum(q_sigma^2)/sum(q_sigma * qt(pstar , df)),
           d = sum((q_sigma - tau * qt(pstar, df))^2)) |>
    ungroup() |> 
    filter(rank(d) == 1)
  
  ## Extract parameters
  beta_tmp_mu <- 0
  beta_tmp_k <- beta_df$df[1]
  beta_tmp_sd_tau <- sigma_df$tau[1]
  beta_tmp_sd_k <- sigma_df$df[1]

  ## gamma -- bernoulli
  if(is.null(pars$gamma))
    gamma_p <- .5
  else
    gamma_p <- pars$gamma$p
  
  ## beta_0 -- half-normal
  if(is.null(pars$beta0)){
    beta0_mu <- 0
    beta0_sd <- 100
  }
  
  ## Return list of prior parameters
  list(beta_tmp_mu = beta_tmp_mu,
       beta_tmp_k = beta_tmp_k,
       beta_tmp_sd_tau = beta_tmp_sd_tau,
       beta_tmp_sd_k = beta_tmp_sd_k,
       gamma_p = gamma_p,
       beta0_mu = beta0_mu,
       beta0_sd = beta0_sd)
}