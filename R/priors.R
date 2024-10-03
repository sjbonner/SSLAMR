prior_simple <- function(pars = NULL, design){

  ## beta -- half-t
  ## Identify parameters based on quantile matching
  
  ## Set default quantiles
  if(is.null(pars$p_beta))  
    p_beta <- c(.5, .99)
  
  if(is.null(pars$q_beta))
    q_beta <- c(1000 * qt(.75,1)/qt(.995,1), 1000)

  ## Set default max df
  if(is.null(pars$max_df))
    max_df <- 30
  
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
    gamma_p <- rep(.5, ncol(design))
  else
    gamma_p <- ifelse(is.na(pars$gamma$p), .5, pars$gamma$p)
  
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

prior_hierarchical <- function(pars = NULL, design){
  ## beta_tmp -- half-t
  
  ## Set default quantiles
  if(is.null(pars$p_beta))  
    p_beta <- c(.5, .99)
  
  if(is.null(pars$q_beta))
    q_beta <- c(1000 * qt(.75,1)/qt(.995,1), 1000)
  
  ## Set default min and max df
  if(is.null(pars$min_df))
    min_df <- 1
  
  if(is.null(pars$max_df))
    max_df <- 30
  
  ## Marginal CDF of beta
  cdf_beta <- function(b, tau, df){
    
    integrand <- function(sigma, tau, df){
      (sigma > 0) * (2*pnorm(b/sigma) - 1) * (2/tau * dt(sigma/tau, df = df))
    }
    
    integrate(integrand, 0, Inf, tau = tau, df = df)
  }
  
  ## Marginal quantile function of beta
  inv_cdf_beta <- function(p, tau, df){
    target <- function(b, tau, df, p){
      cdf_beta(b, tau, df)$value - p
    }
    
    sapply(p, function(p){
      uniroot(target,c(0,1e6), tau = tau, df = df, p = p)$root
    })
  }
  
  ## Solve for tau
  tau_solver <- function(p, q, df){
    target <- function(tau, p, q, df){
      sum((inv_cdf_beta(p, tau, df) - q)^2)
    }
    
    tau_init <- mean(q/inv_cdf_beta(p, 1, df))
    
    optimise(target, c(0.1, 10000), p = p, q = q, df = df)$minimum
  }

  ## Identify parameters for prior on sigma
  sigma_df <- tibble(df = 1:max_df) |> 
    rowwise() |> 
    mutate(tau = tau_solver(p_beta, q_beta, df),
           d = sum((q_beta - inv_cdf_beta(p_beta, tau, df))^2)) |>
    ungroup() |> 
    filter(rank(d) == 1)
  
  ## Extract parameters
  beta_tmp_mu <- 0
  beta_tmp_sd_tau <- sigma_df$tau[1]
  beta_tmp_sd_k <- sigma_df$df[1]

  ## gamma -- bernoulli
  if(is.null(pars$gamma))
    gamma_p <- rep(.5, ncol(design))
  else
    gamma_p <- ifelse(is.na(pars$gamma$p), .5, pars$gamma$p)
  
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
       beta_tmp_sd_tau = beta_tmp_sd_tau,
       beta_tmp_sd_k = beta_tmp_sd_k,
       gamma_p = gamma_p,
       beta0_mu = beta0_mu,
       beta0_sd = beta0_sd)
}