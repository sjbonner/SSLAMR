sslamr_inits <- function(chain,
                         design,
                         width,
                         counts,
                         prior_par,
                         seed = NULL){
  n <- nrow(design)
  K <- ncol(design)
  
  if(!is.null(seed))
    set.seed(seed)
  
  if(chain == 1){
    ## Everything present in equal amounts
    beta0 <- 0.1
    beta_tmp <- rep(mean(counts),K)
    gamma_tmp <- rep(1, K)
  }
  else if(chain == 2){
    ## As little present as possible
    beta0 <- mean(counts/width)
    beta_tmp <- rep(1,K)
    gamma_tmp <- ifelse(prior_par$gamma_p == 1, 1, 0)
  }
  else{
    # Random mixture containing 10% of the predictors plus any that have a prior of 1
    
    # Select 10 percent of the possible predictors
    tmp1 <- sample(1:K, min(K, max(10,floor(.1 * K)))) |> 
      c(which(prior_par$gamma_p == 1)) |> 
      unique() |> 
      sort()
    
    # Fit a Poisson GLM using standard methods
    # Warnings are suppressed since this is likely to produce boundary
    # estimates for some parameters.
    withr::with_options(new = list(warn = -1),
                        {tmp <- glm(counts ~ width + design[,tmp1] - 1, 
                                    family = poisson(link = "identity"),
                                    start = rep(mean(counts),length(tmp1) + 1))})
    
    # Set the initial value for the intercept
    beta0 <- tmp$coefficients["width"]
    
    # Extract the estimates for the remaining coefficients
    b <- tmp$coefficients[names(tmp$coefficients) != "width"]
    
    # Identify which coefficients were positive
    tmp2 <- which(b >= 0) |> 
      c(which(prior_par$gamma_p == 1)) |> 
      unique()
    
    # Set the initial values of gamma and beta for those predictors
    # with positive estimates
    gamma_tmp <- rep(0,K)
    gamma_tmp[tmp1[tmp2]] <- 1
    
    beta_tmp <- rep(1,K)
    beta_tmp[tmp1[tmp2]] <- pmax(b[tmp2],0)
  }
  
  # Return list
  list(gamma_tmp = gamma_tmp,
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
                          n.thin = 1,
                          prior_par = NULL,
                          model = "hierarchical",
                          inits_seed = NULL){
  
  # Argument checks
  if(!model %in% c("hierarchical","log_hierarchical","simple"))
    stop("The argument `model` must be one of 'log_hierarchical', 'hierarchical' or 'simple'.")
  
  # Extract data values
  design <- data$design |> 
    select(-Interval) |> 
    as.matrix()
  
  n <- nrow(design) # Number of groups
  K <- ncol(design) # Number of candidate isotopes
  
  # Identify parents of remaining candidates
  names <- colnames(design)
  
  parent_name <- tibble(ID = names) %>%
    left_join(unique(select(data$candidates, ID, Parent)), by = "ID") %>%
    pull("Parent")
  
  parent <- sapply(parent_name, function(name) which(names == name), simplify = TRUE)
  
  # Extract interval widths
  width <- data$intervals %>% 
    pull("Width")
  
  # Extract counts
  counts <- data$counts$Count  
    
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
                    y = counts,
                    parent = parent)

  # Set prior parameters
  
  if(model == "hierarchical"){
    prior_par <- prior_hierarchical(prior_par, design)
  }
  else if(model == "simple")
    prior_par <- prior_simple(prior_par, design)
  else{
    stop("Model ", model, "is not yet operational.\n")
  }
  
  # Add prior parameters to data
  jags_data <- c(jags_data,prior_par)
 
  # Set initial values
  jags_inits <- lapply(1:n.chains, 
                       sslamr_inits, 
                       design = design, 
                       width = width, 
                       counts = counts,
                       prior_par = prior_par,
                       seed = inits_seed)
  
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
  
  # Common parameters
  monitor <- c("beta","gamma", "mu", "beta0")
  
  # Hierarchical model
  if(model == "hierarchical") 
    monitor <- c(monitor, "beta_tmp_sd")
  
  # Burnin phase
  burnin <- coda.samples(jags_model,
                         variable.names = monitor,
                         n.iter = n.burnin)
  
  # Sampling phase
  samples <- coda.samples(jags_model,
                          variable.names = monitor,
                          n.iter = n.sampling,
                          thin = n.thin)
    
  return(list(inits = jags_inits,
              burnin = burnin,
              samples = samples))
}

#' Spike-and-Slab Analysis for Mass Spectrometry in R
#'
#' @description Implements the Bayesian Poisson generalized linear model of
#'   Bonner, Bonner, Walker and Cao for deconvolution of isotope patterns in
#'   spectra obtained from mass spectrometry using the spike-and-slab prior to
#'   simultaneously identify and quantify compounds of interest within a sample.
#'
#' @section Input Data:
#'
#'   The inputs `spectrum`, `candidates`, `isoinfo`, and `adducts` may either be
#'   specified as strings or data frames. If a string is provided, then it is
#'   interpreted as the path to a file containing the relevant information and
#'   is loaded with [read_any()]. Otherwise, it is assumed to be data frame that
#'   already contains the requisite information in the correct format.
#'
#'   The following list specifies the format for these arguments (whether
#'   contained in an external file or passed as an existing data frame):
#'
#' * `spectrum` two columns: `Mass/Charge` and `Count`.
#' * `candidates` three columns: `Name`, `Formula`, and `Charge`.
#' * `isoinfo ` four columns `element`, `mass`, `abundance`, and `nucleons`.
#'   Abundances for the same element must sum to 1.
#' * `adducts` may either be a string or a data frame. If a string then it is
#'   interpreted as the path to a file containing the observed spectrum and is
#'   read with [read_any()]. Otherwise, it is assumed to be a data frame. The
#'   data must contain three columns: `Name`, `Formula`, and `Action`. Action
#'   may either be `-` indicating that the fragment may be removed from the
#'   candidates or `+` indicating that the adduct may be added to the
#'   candidates.
#'
#' @section Isotope Data:
#'
#'   The `isotope_data` may either be `NULL` (default), a string, or a data
#'   frame. If `NULL` then the relative abundances for the isotope of the
#'   candidates are computed with the [ecipex()] function. If a string then it
#'   is interpreted as the path to a file containing the processed isotope data
#'   and is read with [read_any()]. Otherwise, it is assumed to be a data frame.
#'   The data must contain five columns: `ID`, `Charge`, `Isotope`, `Mass`, and
#'   `Abund`.
#'
#' @section Models:
#'
#'   The `model` argument determines the structure of the prior (the model) for
#'   the abundances of candidates present in the sample. If `simple` then
#'   abundances are assigned independent prior distributions with fixed
#'   parameters. If `hierarchical` then the abundances are modelled as draws
#'   from a distribution whose parameters are determined by the data and
#'   assigned further prior distributions.
#'
#' @section Prescreening:
#'
#'   Prescreening allows some candidates to be removed from the analysis before
#'   fitting the model. Let \eqn{y_i} denote the count for the i-th entry in the
#'   observed spectrum and \eqn{x_{ij}} the relative contribution of the j-th
#'   candidate to this value. If `prescreen_weight` is `FALSE` then the
#'   candidate is include in modelling only if the sum of the counts in the
#'   spectrum associated with the candidate, \deqn{\sum_{i=1}^n y_i 1(x_{ij}
#'   >0),} is greater than `prescreen`. If `prescreen_weight` is `TRUE` then 
#'   counts are weighted by the relative abundance of the isotopes, scaled so 
#'   that the most common isotope has a weight of 1. The
#'   candidate is include in modelling only if the weighted sum,
#'   \deqn{\sum_{i=1}^n y_i \frac{x_{ij}}{\max_{j}x_{ij}}} is greater than `prescreen`. Setting
#'   `prescreen=0` will include all candidates in the model.
#'
#' @section Rounding:
#'
#'   Intensities in the observed spectrum (`Count`) are assumed to be integer
#'   counts. If non-integer values are detected then these may be rounded in
#'   different ways. Options are:
#'   * `nearest` round to the nearest integer in the usual way.
#'   * `floor` round down to the next highest integer.
#'   * `ceiling` round up to the next lowest integer.
#'   * `bernoulli` if \eqn{d} is the decimal part of the count then round up
#'   with probability \eqn{d} and down with probability \eqn{1-p}.
#'
#'   The advantage of `bernoulli` is that the expected value of the sum of the
#'   rounded values is equal to the sum of the original values.
#'
#' @param spectrum Observed spectrum. May be a string specifying the path to the
#'   data file or an existing data frame. See Input Details for further
#'   information. (character or dataframe)
#' @param candidates List of candidate compounds. May be a string specifying the
#'   path to the data file or an existing data frame. See Input Details for
#'   further information. (character or dataframe)
#' @param n.adapt Number of iterations in MCMC sampler's adapting phase.
#'   (numeric)
#' @param n.chains Number of parallel chains in MCMC sampler. (numeric)
#' @param verbose If TRUE then print progress messages. (boolean)
#' @param min_abundance Minimum proportion retained within candidate isotope
#'   patterns. (numeric)
#' @param epsilon Tolerance of binning algorithm. (numeric)
#' @param n.burnin Number of iterations in MCMC sampler's burn-in phase.
#'   (numeric)
#' @param n.sampling Number of iterations in the MCMC sampler's sampling phase.
#'   (numeric)
#' @param run_model If TRUE then run MCMC sampler. If FALSE then prepare input
#'   without runing sampler. Boolean.
#' @param prior_par List of lists specifying the parameters of the prior
#'   distributions. List.
#' @param isoinfo Custom element isotope information that either replaces or
#'   adds to the definitions in [ecipex()] (see replace_isoinfo). May be a
#'   string specifying the path to the data file or an existing data frame. See
#'   Input Details for further information. (character or dataframe)
#' @param isotope_data Preprocessed candidate isotope information. May be a
#'   string specifying the path to the data file or an existing data frame. See
#'   Isotope Data for further information. (character or dataframe)
#' @param xlsx_out Path including name of output file. (character)
#' @param adducts List of adducts. May be a string specifying the path to the
#'   data file or an existing data frame. See Input Details for further
#'   information. (character or dataframe)
#' @param replace_isoinfo If TRUE then the information in `isoinfo` replaces the
#'   data on the ratios of the elemental isotopes used by [ecipex()]. This
#'   overrides the default isotope ratios. Otherwise, the information in
#'   `isoinfo` is appended to the data used by [ecipex()]. (boolean)
#' @param group_formula If TRUE then group candidates with the same chemical
#'   formula. (boolean)
#' @param max_isotopes Maximum number of isotopes retained in the pattern of any
#'   candidate. (numeric)
#' @param reload_results If TRUE then MCMC sampler is not run and output is
#'   retrieved from the file specified by `xlsx_out`. Otherwise any existing
#'   file at this location will be overwritten. (boolean)
#' @param skip_isotopes Number of isotopes at the start of the isotope pattern
#'   that are ignored. Defaults to 0 implying that all isotopes with relative
#'   abundance greater than `min_abundance` up to `max_isotopes` will be
#'   retained in the patter of any candidate. (numeric)
#' @param min_mass_charge Lower bound of the mass charge window for the
#'   analysis. Any observations in the spectrum or isotopes of candidates with
#'   mass/charge ratio below this value are ignored. (numeric)
#' @param max_mass_charge Upper bound of the mass charge window for the
#'   analysis. Any observations in the spectrum or isotopes of candidates with
#'   mass/charge ratio above this value are ignored. (numeric)
#' @param binning Specifies how the peaks in the spectrum and the entries in the
#'   candidate isotope patterns are associated to construct the design matrix.
#'   If TRUE then observations in the spectrum and elements of the isotope
#'   pattern are binned. Otherwise, elements of each isotope pattern are
#'   associated with the nearest observation in the spectrum up to a distance of
#'   `epsilon`. (boolean)
#' @param model Specifies the form of the model for the abundances of candidates
#'   present in the sample. Options are "simple" or "hierarchical". See Models
#'   for further information. (character)
#' @param prescreen Minimum count for prescreening candidate. See Prescreening
#'   for further information. (numeric)
#' @param rounding Method for rounding non-integer counts. Options are
#'   "nearest", "floor", "ceiling", or "bernoulli". See Rounding for further
#'   information. (character)
#' @param prescreen_weight If TRUE then counts are weighted be the relative
#'   isotope abundance during the prescreening phase. See Prescreening for
#'   further information. (boolean)
#' @param mixtures If TRUE then summarize information about the sampled 
#'   mixtures. (boolean)
#' @param prescreen_prior Prescreening threshold for prior inclusion probability. 
#' Any candidate with prior inclusion probability greater than or equal to this 
#' value will be retained in the analysis regardless of the counts in the associated 
#' intervals. (numeric)
#' @param n.thin Thinning parameter for MCMC sampler. (integer)
#' @param group_pattern If TRUE then candidates with similar isotope pattens are 
#' grouped when computing posterior summary statistics for presence and 
#' abundance. (boolean)
#' @param pattern_tol Tolerance for grouping candidates with similar isotope patterns. 
#' If `group_pattern` is true then candidates with a maximum difference of their 
#' isotope patterns less than `pattern_tol` are grouped when computing posterior 
#' summary statistics for the presence and abundance. (numeric)
#'
#' @importFrom writexl write_xlsx
#' @importFrom tictoc tic toc
#' @return A list of objects produced by processing the data and fitting the
#'   model. Objects marked (ro) are only included if `run` is `TRUE`.
#'   * `data`
#'   * `convergence` (ro)
#'   * `burnin` (ro)
#'   * `samples` (ro)
#'   * `coefficients` (ro)
#'   * `intercept` (ro)
#'   * `fitted` (ro)
#'   * `parameters` (ro)
#' @export
#'
#' @examples
#' #TBD
sslamr <- function(spectrum = NULL,
                   candidates = NULL,
                   isotope_data = NULL,
                   adducts = NULL,
                   isoinfo = NULL,
                   replace_isoinfo = FALSE,
                   min_abundance = .001,
                   max_isotopes = Inf,
                   skip_isotopes = 0,
                   group_formula = NULL,
                   group_pattern= TRUE,
                   pattern_tol = .05,
                   binning = TRUE,
                   epsilon = .05,
                   min_mass_charge = NULL,
                   max_mass_charge = NULL,
                   prescreen = 0,
                   prescreen_prior = 1,
                   prescreen_weight = FALSE,
                   rounding = "nearest",
                   n.chains = 3,
                   n.adapt = 1000,
                   n.burnin = 1000,
                   n.sampling = 10000,
                   n.thin = 1,
                   prior_par = NULL,
                   model = "hierarchical",
                   inits_seed = NULL,
                   run_model = TRUE,
                   mixtures = FALSE,
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
  
    # Save input arguments for writing to output
    parameters <- mget(ls(environment(), sorted=F)) |> 
      c(match.call(expand.dots=F)$...) |>
      lapply(\(x) ifelse(is.null(x), "", x)) |> 
      as_tibble() |> 
      mutate(across(everything(), as.character)) |>
      pivot_longer(everything(), names_to = "Argument", values_to = "Value")
    
    # Add SSLAMR version
    parameters <- tibble(Argument = "SSLAMR Version", 
                         Value = as.character(packageVersion("SSLAMR"))) |> 
      bind_rows(parameters)
    
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
      if(verbose) message("   Loading element isotope information (",isoinfo,")")
      isoinfo <- read_any(isoinfo)
    }
    
    # Add element isotope information to ecipex's existing list
    if(!replace_isoinfo)
      isoinfo <- isoinfo %>%
        bind_rows(ecipex::nistiso)
    
    # Add Prior column to candidate list if not specified
    if(!("Prior" %in% colnames(candidates))){
      candidates <- candidates %>%
        add_column(Prior = NA)
    }
    
    
    # Modify candidates via adducts
    if(!is.null(adducts)){
      candidates <- candidates %>%
        modify_adducts(adducts)
    }
    else{
      candidates <- candidates %>%
        mutate(Parent = Name)
    }
    
    if(verbose) message("Processing data...")
    tic()
    
    # Group candidates with identical chemical formulas
    if(!is.null(group_formula)){
      # Group candidates has been set by user. Check setting is applicable.
      if(group_formula & !is.null(adducts))
        stop("Error: Candidates cannot be grouped when adducts are included. Please set group_formula to FALSE.\n")

    }
    else{
      # Group candidates not set by user but not applicable.
      if(!is.null(adducts)){
        if(verbose) message("   Adducts provided. Not grouping candidates...")
        group_formula <- FALSE
      }
      else if(any(!is.na(candidates$Prior))){
        if(verbose) message("   Adducts provided. Not grouping candidates...")
        group_formula <- FALSE
      }
      else{
        group_formula <- TRUE
      }
    }
    
    if(group_formula){
      # Split candidates list
      candidates1 <- candidates %>%
        filter(!is.na(Prior))
      
      # Group candidates
      candidates2 <- candidates %>%
        filter(is.na(Prior)) %>%
        group_by(Formula,Charge) %>%
        summarize(Name = paste(Name,collapse = "/"),
                  .groups = "drop")
      
      # Rejoin list
      candidates <- candidates1 %>%
        bind_rows(candidates2) %>%
        arrange(Name)
    }
    
    # Pull prior information from candidates data
    if(has_name(candidates,"Prior")){
      gamma_p <- pull(candidates,"Prior")
      names(gamma_p) <- pull(candidates,"Name")
      
      prior_par$gamma <- list(p = gamma_p)
      
      candidates <- candidates %>%
        select(-Prior)
    }
    
    # define design matrix
    data <- sslamr_data(spectrum, 
                        candidates = candidates,
                        isotope_data = isotope_data,
                        min_abundance = min_abundance,
                        max_isotopes = max_isotopes,
                        skip_isotopes = skip_isotopes,
                        binning = binning,
                        epsilon = epsilon, 
                        min_mass_charge = min_mass_charge,
                        max_mass_charge = max_mass_charge,
                        prior_par = prior_par,
                        prescreen = prescreen,
                        prescreen_prior = prescreen_prior,
                        prescreen_weight = prescreen_weight,
                        rounding = rounding,
                        pattern_tol = pattern_tol,
                        isoinfo = isoinfo,
                        verbose = verbose)
    toc_out <- toc(quiet = !verbose)
    timing <- tibble(Stage = "Processing data",
                     Time = toc_out$toc - toc_out$tic)
    
    # Modify prior parameters to account for pre-screening
    if(!is.null(prior_par$gamma$p)){
      prior_par$gamma$p <- prior_par$gamma$p[colnames(data$design[-1])]
    }

    # MCMC sampling
    if(run_model){
      if(verbose) message("Running sampler...")
      tic()
      ss.model <- sslamr_sample(data,
                                n.adapt = n.adapt,
                                n.chains = n.chains,
                                n.burnin = n.burnin,
                                n.sampling = n.sampling,
                                n.thin = n.thin,
                                prior_par = prior_par,
                                model = model,
                                inits_seed = inits_seed)
      toc_out <- toc(quiet = !verbose)
      timing <- timing |> 
        add_row(Stage = "Running sampler",
                Time = toc_out$toc - toc_out$tic)

      # Convert burnin to data frame
      b.df <- get_samples_df(ss.model$burnin)
      
      # Convert samples to a data frame
      s.df <- get_samples_df(ss.model$samples) %>%
        pivot_longer(-c(Chain,Iteration),names_to = "Parameter",values_to = "Value") %>%
        mutate(ID = ifelse(Parameter %in% c("beta0","beta_tmp_sd"), NA, 
                       as.integer(str_extract(Parameter,"[0-9]+"))),
               Parameter = ifelse(Parameter %in% c("beta0","beta_tmp_sd"), Parameter,
                                  str_extract(Parameter,"[a-z]+"))) %>%
        mutate(Name = ifelse(Parameter %in% c("beta","gamma"),colnames(data$design)[ID+1],NA))
      
      # results
      if(verbose) message("Summarizing results...")
      
      # Identify groups by pattern
      if(group_pattern){
        if(verbose)
          message("   Identyfing similar isotope patterns...")
        
        # Remove candidates which never appear
        not_present <- s.df |> 
          filter(Parameter == "gamma") |> 
          group_by(Name,ID) |> 
          summarize(Any = any(Value > 0),
                    .groups = "drop") |>
          filter(!Any)
          
        # Group remaining candidates
        groups <- group_by_pattern(data$design[-c(1,pull(not_present,"ID") + 1)], 
                                   tol = pattern_tol)
        
        # Recombine
        data$groups <- groups |>
          select(-Group_ID) |> 
          bind_rows(tibble(Name = not_present |> pull("Name"),
                           Group_Name = not_present |> pull("Name"))) |> 
          arrange(Group_Name) |> 
          group_by(Group_Name) |>
          mutate(Group_ID = cur_group_id())
      }
      else{
        # Assign each candidate to its own group
        data$groups <- tibble(Name = colnames(design)) |> 
          mutate(Group_ID = 1:n(),
                 Group_Name = Name)
      }
      
      if(verbose) message("    Summarizing beta and gamma")
      tic() 
      bg.summ <- beta.gamma_summ(s.df, 
                                 design = data$design, 
                                 groups = data$groups)
      
      toc_out <- toc(quiet = !verbose)
      timing <- timing |> 
        add_row(Stage = "Summarizing beta and gamma",
                Time = toc_out$toc - toc_out$tic)
      
      if(mixtures){
        if(verbose) message("    Summarizing mixtures... be patient...")
        tic() 
        mixtures.summ <- mixtures_summ(
          samp_df = s.df, 
          design = data$design,
          groups = data$groups)
        
        toc_out <- toc(quiet = !verbose)
        timing <- timing |> 
          add_row(Stage = "Summarizing mixtures",
                  Time = toc_out$toc - toc_out$tic)
      }
      
      if(model == "hierarchical"){
        if(verbose) message("    Computing standard deviations")
        tic() 
        beta_sd.summ <- beta_sd_summ(s.df)
        toc_out <- toc(quiet = !verbose)
        timing <- timing |>
          add_row(Stage = "Computing standard deviations",
                  Time = toc_out$toc - toc_out$tic)
      }
      
      if(verbose) message("    Summarizing intercept")
      tic() 
      int.summ <- intercept_summ(s.df)
      toc_out <- toc(quiet = !verbose)
      timing <- timing |>
        add_row(Stage = "Summarizing intercept",
                Time = toc_out$toc - toc_out$tic)
      
      if(verbose) message("    Summarizing fitted values")
      tic() 
      fit.summ <- fitted_summ(s.df, data$counts)
      toc_out <- toc(quiet = !verbose)
      timing <- timing |>
        add_row(Stage = "Summarizing fitted values",
                Time = toc_out$toc - toc_out$tic)
      
      # Convergence diagnostics
      if(verbose) message("Computing convergence diagnostics...")
      tic() 
      
      index <- grep("mu",colnames(ss.model$samples[[1]]))
      
      mu_burnin <- ss.model$burnin[,index] |> 
        window(start = n.adapt + floor(n.burnin/2))
      
      convergence <- coda::gelman.diag(mu_burnin,
                                       multivariate = FALSE)
      
      # effective sizes
      mu_samples <- lapply(ss.model$samples,function(mcmc){
        index <- grep("mu",colnames(ss.model$samples[[1]]))
        mcmc[,index]
      })
      
      eff.size <- coda::effectiveSize(mu_samples)
      
      fit.summ <- as_tibble(cbind(fit.summ, EffectiveSize=eff.size))
      toc_out <- toc(quiet = !verbose)
      timing <- timing |>
        add_row(Stage = "Computing convergence diagnostics",
                Time = toc_out$toc - toc_out$tic)
    }

  # package results
  if(run_model){
    results <- list(data=data,
                    convergence = as_tibble(convergence$psrf, rownames = "Parameter"),
                    burnin = b.df,
                    samples = s.df,
                    coefficients=bg.summ,
                    intercept = int.summ,
                    fitted = fit.summ,
                    parameters = parameters,
                    timing = timing)
    
    if(mixtures)
      results$mixtures <- mixtures.summ
    
    if(model == "hierarchical"){
      results$beta_sd <- beta_sd.summ
      
    }
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
                              fitted = fit.summ,
                              parameters = parameters,
                              convergence = as_tibble(convergence$psrf, rownames = "Parameter"),
                              timing = timing))
        
        if(mixtures)
          xlsx_output$mixtures <- mixtures.summ
        
        if(model == "hierarchical")
          xlsx_output$beta_sd <- beta_sd.summ
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
    
    # Retrieve base output objects
    results <- list(data = list(candidates = read_xlsx(xlsx_out,sheet = "candidates"),
                                intervals = read_xlsx(xlsx_out,sheet = "intervals"),
                                spectrum = read_xlsx(xlsx_out,sheet = "spectrum"),
                                design = read_xlsx(xlsx_out,sheet = "design"),
                                counts = read_xlsx(xlsx_out,sheet = "counts")),
                    coefficients = read_xlsx(xlsx_out,sheet = "coefficients"),
                    intercept = read_xlsx(xlsx_out,sheet = "intercept"),
                    fitted = read_xlsx(xlsx_out,sheet = "fitted"),
                    parameters = read_xlsx(xlsx_out,sheet = "parameters"),
                    convergence = read_xlsx(xlsx_out,sheet = "convergence"),
                    timing = read_xlsx(xlsx_out,sheet = "timing"))
    
    # Retrieve optional output objects
    if("beta_sd" %in% excel_sheets(xlsx_out))
      results$beta_sd <- read_xlsx(xlsx_out,sheet = "beta_sd")
    
    if("mixtures" %in% excel_sheets(xlsx_out))
      results$mixtures <- read_xlsx(xlsx_out,sheet = "mixtures")
  }
  
  # Return output
  return(results)
}
  