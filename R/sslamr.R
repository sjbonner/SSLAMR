sslamr_inits <- function(chain,
                         design,
                         width,
                         counts){
  n <- nrow(design)
  K <- ncol(design)
  
  if(chain == 1){
    ## Everything present in equal amounts
    beta0 <- 0.1
    beta_tmp <- rep(mean(counts),K)
    gamma <- rep(1, K)
  }
  else if(chain == 2){
    ## Nothing present
    beta0 <- mean(counts/width)
    beta_tmp <- rep(1,K)
    gamma <- rep(0,K)
  }
  else{
    # Random mixture containing 10% of the predictors
    
    # Select 10 percent of the possible predictors
    tmp1 <- sort(sample(1:K, min(K, max(10,floor(.1 * K)))))
    
    # Fit a Poisson GLM using standard methods
    tmp <- glm(counts ~ width + design[,tmp1] - 1, 
               family = poisson(link = "identity"),
               start = rep(mean(counts),length(tmp1) + 1))
    
    # Set the initial value for the intercept
    beta0 <- tmp$coefficients["width"]
    
    # Extract the estimates for the remaining coefficients
    b <- tmp$coefficients[names(tmp$coefficients) != "width"]
    
    # Identify which coefficients were positive
    tmp2 <- which(b >= 0)
    
    # Set the initial values of gamma and beta for those predictors
    # with positive estimates
    gamma <- rep(0,K)
    gamma[tmp1[tmp2]] <- 1
    
    beta_tmp <- rep(1,K)
    beta_tmp[tmp1[tmp2]] <- b[tmp2]
  }
  
  # Return list
  list(gamma = gamma,
       beta0 = beta0,
       beta_tmp = beta_tmp)
}

# Model running function ####

# this functions takes the name of the file and runs the below model
# designed to be looped across a list of files
# returns the list of MCMC samples from JAGS, which can be moved to
# other summary functions

sslamr_sample <- function(data, 
                          n.chains = 1,
                          n.adapt = 1000, 
                          n.burnin = 1000,
                          n.sampling = 10000,
                          prior_par = NULL,
                          model = "hierarchical"){
  
  # Argument checks
  if(!model %in% c("hierarchical","log_hierarchical","simple"))
    stop("The argument `model` must be one of 'log_hierarchical', 'hierarchical' or 'simple'.")
  
  # Set default prior parameters
  if(is.null(prior_par$beta0))
    prior_par$beta0 <- list(mu = 0, sd = 100)
  if(is.null(prior_par$beta_tmp))
    prior_par$beta_tmp <- list(mu = 0)
  if(is.null(prior_par$beta_tmp_sd))
    prior_par$beta_tmp_sd <- list(k = 4, tau = (qt(.75,4)/100)^2)
  if(is.null(prior_par$gamma))
    prior_par$gamma <- list(p = .5)
  
  # Extract data values
  design <- data$data %>%
    select(-Interval, -Isotopes, - Mass, -Lower, 
           -Upper, -Width, - Peaks, -Count) %>%
    as.matrix()
  
  n <- nrow(design) # Number of groups
  K <- ncol(design) # Number of candidate isotopes
  
  width <- data$data %>% # Extract interval widths
    pull("Width")
  
  counts <- data$data %>% # Extract counts
    pull("Count") 
  
  ## Identify non-zero entries of the design matrix
  slot_list <- which(design > 0, arr.ind = TRUE) %>%
    as_tibble() %>%
    arrange(row)
  
  nslots <- slot_list %>%
    group_by(row) %>%
    summarize(n = n()) %>%
    pull("n")
  
  slot_mat <- slot_list %>%
    group_by(row) %>%
    mutate(i = 1:n()) %>%
    pivot_wider(names_from = i, values_from = col) %>%
    ungroup() %>%
    select(-row) %>%
    as.matrix()
  
  index1 <- unique(pull(slot_list,"row"))
  
  index2 <- (1:nrow(design))[-index1]
    
  jags_data <- list(n1 = length(index1),
                    n2 = length(index2),
                    index1 = index1,
                    index2 = index2,
                    K = K,
                    width = width,
                    X = design[index1,],
                    slots = slot_mat,
                    nslots = nslots,
                    y = counts)
  
  # Add prior parameters to data
  if(model %in% c("hierarchical","log_hierarchical"))
    jags_data <- c(jags_data,
                   beta0_mu = prior_par$beta0$mu,
                   beta0_sd = prior_par$beta0$sd,
                   beta_tmp_mu = prior_par$beta_tmp$mu,
                   beta_tmp_sd_k = prior_par$beta_tmp_sd$k,
                   beta_tmp_sd_tau = prior_par$beta_tmp_sd$tau,
                   gamma_p = prior_par$gamma$p)
  else if(model == "simple")
    jags_data <- c(jags_data,
                   beta0_mu = prior_par$beta0$mu,
                   beta0_sd = prior_par$beta0$sd,
                   beta_tmp_k = prior_par$beta_tmp$k,
                   beta_tmp_sd = prior_par$beta_tmp$sd,
                   gamma_p = prior_par$gamma$p)
  
  # Set initial values
  jags_inits <- lapply(1:n.chains, sslamr_inits, design = design, width = width, counts = counts)
  
  # model
  if(model == "hierarchical")
    model_file <- system.file("JAGS/spike_and_slab_jags_hierarchical.R",package="SSLAMR")
  else if(model == "log_hierarchical")
    model_file <- system.file("JAGS/spike_and_slab_jags_log_hierarchical.R",package="SSLAMR")
  else if(model == "simple")
    model_file <- system.file("JAGS/spike_and_slab_jags_simple.R",package="SSLAMR")
  else
    stop("Unknown model ",model,".")
  
  # Initialize model and run adapting phase
  jags_model <- jags.model(model_file,
                           data = jags_data,
                           inits = jags_inits,
                           n.chains = n.chains,
                           n.adapt = n.adapt)
  
  # Define variables to monitor
  monitor <- c("beta","gamma", "mu", "beta0","beta_tmp_sd")
  
  # Burnin phase
  burnin <- coda.samples(jags_model,
                          variable.names = monitor,
                          n.iter = n.burnin)
  
  # Sampling phase
  samples <- coda.samples(jags_model,
                          variable.names = monitor,
                          n.iter = n.sampling)
    
  return(list(inits = jags_inits,
              burnin = burnin,
              samples = samples))
}

# All in one wrapper
#' Title
#'
#' @param spectrum
#' @param candidates 
#' @param n.adapt 
#' @param n.chains 
#' @param verbose 
#' @param min_abundance 
#' @param epsilon 
#' @param n.burnin 
#' @param n.sampling 
#' @param run_model 
#' @param prior_par 
#' @param isoinfo 
#' @param isotope_data 
#' @param xlsx_out 
#' @param adducts 
#' @param replace_isoinfo 
#' @param group_candidates If TRUE then group candidates with the same chemical formula.
#' @param max_isotopes 
#' @param reload_results 
#' @param skip_isotopes 
#' @param min_mass_charge 
#' @param max_mass_charge 
#'
#' @importFrom writexl write_xlsx
#' @importFrom tictoc tic toc
#' @return
#' @export
#'
#' @examples
sslamr <- function(spectrum = NULL,
                   candidates = NULL,
                   isotope_data = NULL,
                   adducts = NULL,
                   isoinfo = NULL,
                   replace_isoinfo = FALSE,
                   min_abundance = .001,
                   max_isotopes = Inf,
                   skip_isotopes = 0,
                   group_candidates = TRUE,
                   epsilon = .05,
                   min_mass_charge = NULL,
                   max_mass_charge = NULL,
                   n.chains = 3,
                   n.adapt = 1000,
                   n.burnin = 1000,
                   n.sampling = 10000,
                   prior_par = NULL,
                   model = "hierarchical",
                   run_model = TRUE,
                   reload_results = FALSE,
                   xlsx_out = NULL,
                   verbose = TRUE){
  
  if(!reload_results){
    # Check input
    if (is.null(spectrum))
      stop("You must supply the spectrum data.")
    
    if(is.null(candidates) & is.null(isotope_data))
      stop("You must supply one of either the list of candidates or the formatted isotope data.")
    
    if((!is.null(candidates) & !is.null(isotope_data)))
      stop("You cannot supply both the list of candidates and the formatted isotope data.")
    
    if(!is.null(isotope_data) & !is.null(adducts))
      stop("You cannot supply a list of adducts with already formatted isotope data.")
    
    # Read input files
    if(verbose) message("Loading input data...")
    
    # 1) Spectrum data
    if(is.character(spectrum)){
      if(verbose) message("   Loading spectrum data...")
      spectrum <- read_any(spectrum)
    }  
    
    # 2) Candidate data
    if(is.character(candidates)){
      if(verbose) message("   Loading candidates data (",candidates,")")
      candidates <- read_any(candidates)
    }
    
    # 3) Isotope data
    if(is.character(isotope_data)){
      if(verbose) message("   Loading isotope data (",isotope_data,")")
      isotope_data <- read_any(isotope_data)
    }
    
    # 4) Adduct data
    if(is.character(adducts)){
      if(verbose) message("   Loading adducts data (",adducts,")")
      adducts <- read_any(adducts)
    }
    
    # 5) Element isotope information
    if(is.character(isoinfo)){
      if(verbose) message("    Loading element isotope information (",isoinfo,")")
      isoinfo <- read_any(isoinfo)
    }
    
    # Add element isotope information to ecipex's existing list
    if(!replace_isoinfo)
      isoinfo <- isoinfo %>%
        bind_rows(ecipex::nistiso)
    
    # Modify candidates via adducts
    if(!is.null(adducts))
      candidates <- modify_adducts(candidates , adducts)
    
    # Group candidates with identical chemical formulas
    if(group_candidates)
      candidates <- candidates %>%
        group_by(Formula,Charge) %>%
        summarize(Name = paste(Name,collapse = "/"))
    
    # define design matrix
    if(verbose) message("Processing data...")
    
    data <- sslamr_data(spectrum, 
                        candidates = candidates,
                        isotope_data = isotope_data,
                        min_abundance = min_abundance,
                        max_isotopes = max_isotopes,
                        skip_isotopes = skip_isotopes,
                        epsilon = epsilon, 
                        min_mass_charge = min_mass_charge,
                        max_mass_charge = max_mass_charge,
                        isoinfo = isoinfo,
                        verbose = verbose)
    
    # spike slab sampling
    if(run_model){
      if(verbose) message("Running sampler...")
      ss.model <- sslamr_sample(data,
                                n.adapt = n.adapt,
                                n.chains = n.chains,
                                n.burnin = n.burnin,
                                n.sampling = n.sampling,
                                prior_par = prior_par,
                                model = model)
      
      # Convert burnin to data frame
      b.df <- get_samples_df(ss.model$burnin)
      
      # Convert samples to a data frame
      s.df <- get_samples_df(ss.model$samples)
      
      # results
      if(verbose) message("Summarizing results...")
      
      
      if(verbose) message("    Summarizing beta and gamma")
      if(verbose) tic() 
      bg.summ <- beta.gamma_summ(s.df, data$candidates)
      if(verbose) toc()
      
      if(verbose) message("    Computing standard deviations")
      if(verbose) tic() 
      beta_sd.summ <- beta_sd_summ(s.df)
      if(verbose) toc()
      
      if(verbose) message("    Summarizing intercept")
      if(verbose) tic() 
      int.summ <- intercept_summ(s.df)
      if(verbose) toc()
      
      if(verbose) message("    Summarizing fitted values")
      if(verbose) tic() 
      fit.summ <- fitted_summ(s.df, data$data)
      if(verbose) toc()
      
      # Convergence diagnostics
      if(verbose) message("Computing convergence diagnostics...")
      if(verbose) tic() 
      mu_burnin <- lapply(ss.model$burnin,function(mcmc){
        index <- grep("mu",colnames(ss.model$burnin[[1]]))
        mcmc[,index]
      })
      
      convergence <- coda::gelman.diag(mu_burnin,
                                       multivariate = FALSE)
      
      # effective sizes
      mu_samples <- lapply(ss.model$samples,function(mcmc){
        index <- grep("mu",colnames(ss.model$samples[[1]]))
        mcmc[,index]
      })
      
      eff.size <- coda::effectiveSize(mu_samples)
      
      fit.summ <- as_tibble(cbind(fit.summ, EffectiveSize=eff.size))
      if(verbose) toc()
    }

  # package results
  parameters <- tibble(min_abundance = min_abundance,
                       epsilon = epsilon,
                       min_mass_charge = min_mass_charge,
                       max_mass_charge = max_mass_charge,
                       n_chains = n.chains,
                       n_adapt = n.adapt,
                       n_burnin = n.burnin,
                       n_sampling  = n.sampling) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value")

  if(run_model){
    results <- list(data=data,
                    convergence = convergence,
                    burnin = b.df,
                    samples = s.df,
                    coefficients=bg.summ,
                    beta_sd = beta_sd.summ,
                    intercept = int.summ,
                    fitted = fit.summ,
                    parameters = parameters)
  }
  else {
    results <- list(data = data,
                   parameters = parameters)
  }
    
  if(!is.null(xlsx_out)){
    ## Save formatted output to Excel spreadsheet
      if(run_model){
        xlsx_output <- c(results$data,
                         list(coefficients = bg.summ,
                              intercept = int.summ,
                              beta_sd = beta_sd.summ,
                              fitted = fit.summ,
                              parameters = parameters))
      }
      else{
        xlsx_output <- c(results$data,
                         list(parameters = parameters))
      }
      
      write_xlsx(xlsx_output,xlsx_out)
    }
  }
  
  else{
    # Reload existing results
    
    # Check for output file
    if(is.null(xlsx_out))
      stop("Please supply the name of an output file.")
    
    if(!file.exists(xlsx_out))
      stop("The file ",xlsx_out,"does not exist.")
    
    results <- list(candidates = read_xlsx(xlsx_out,sheet = "candidates"),
                    intervals = read_xlsx(xlsx_out,sheet = "intervals"),
                    spectrum = read_xlsx(xlsx_out,sheet = "spectrum"),
                    data = read_xlsx(xlsx_out,sheet = "data"),
                    coefficients = read_xlsx(xlsx_out,sheet = "coefficients"),
                    intercept = read_xlsx(xlsx_out,sheet = "intercept"),
                    beta_sd = read_xlsx(xlsx_out,sheet = "beta_sd"),
                    fitted = read_xlsx(xlsx_out,sheet = "fitted"),
                    parameters = read_xlsx(xlsx_out,sheet = "parameters"))
  }
  
  # Return output
  return(results)
}
  