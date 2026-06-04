sslamr_summarize <- function(model,
                             data,
                             group_pattern,
                             pattern_tol){
  # Convert burnin to data frame
  b.df <- get_samples_df(model$burnin)
  
  # Convert samples to a data frame
  s.df <- get_samples_df(model$samples) %>%
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
  coefficient.summ <- coefficient_summ(s.df, 
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
  
  index <- grep("mu",colnames(model$samples[[1]]))
  
  mu_burnin <- model$burnin[,index] |> 
    window(start = n.adapt + floor(n.burnin/2))
  
  convergence <- coda::gelman.diag(mu_burnin,
                                   multivariate = FALSE)
  
  # effective sizes
  mu_samples <- lapply(model$samples,function(mcmc){
    index <- grep("mu",colnames(model$samples[[1]]))
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
                  coefficients=coefficient.summ,
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
                     list(coefficients = coefficient.summ,
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
get_samples_df <- function(samples){
  sample_df <- lapply(1:length(samples),function(k){
    as_tibble(samples[[k]]) %>%
      add_column(Chain = k, Iteration = 1:nrow(.), .before = 1)
  }) %>%
    bind_rows()
  return(sample_df)
}

beta_sd_summ <- function(samp_df){
  
  samples_beta_sd <- samp_df %>%
    filter(Parameter == "Value") %>%
    select(Chain, Iteration, Value)
  
  summ_beta_sd <- samples_beta_sd %>%
    summarize(Mean = mean(Value, na.rm = TRUE),
              Median = median(Value, na.rm = TRUE),
              SD = sd(Value, na.rm = TRUE),
              Q2.5 = quantile(Value,.025,na.rm = TRUE),
              Q25 = quantile(Value,.25,na.rm = TRUE),
              Q75 = quantile(Value,.75,na.rm = TRUE),
              Q97.5 = quantile(Value,.975,na.rm = TRUE))
  
  return(summ_beta_sd)
}

intercept_summ <- function(samp_df){
  samples_beta0 <- samp_df %>%
    filter(Parameter == "beta0") %>%
    select(Chain, Iteration, Value)
  
  summ_beta0 <- samples_beta0 %>%
    summarize(Mean = mean(Value, na.rm = TRUE),
              Median = median(Value, na.rm = TRUE),
              SD = sd(Value, na.rm = TRUE),
              Q2.5 = quantile(Value,.025,na.rm = TRUE),
              Q25 = quantile(Value,.25,na.rm = TRUE),
              Q75 = quantile(Value,.75,na.rm = TRUE),
              Q97.5 = quantile(Value,.975,na.rm = TRUE))
  
  return(summ_beta0)
}

coefficient_summ <- function(samp_df, design, groups){
  
  ## Extract names of candidates
  candidate_names <- colnames(design)[-1]
  
  ## Transform from array to tibble
  samples_beta_gamma <- samp_df %>%
    filter(Parameter %in% c("beta","gamma")) |>
    pivot_wider(names_from = Parameter, values_from = Value)

  ## Remove candidates with probability of presence of zero
  samples_beta_gamma <- samples_beta_gamma |> 
    group_by(ID) |> 
    mutate(P = mean(gamma)) |> 
    ungroup()
  
  summ_beta_gamma_1 <- samples_beta_gamma |> 
    filter(P == 0) |> 
    group_by(ID) |> 
    summarize(ID = first(ID),
              P_present = 0,
              P_absent = 1 ,
              Mean = NA,
              Median = NA,
              SD = NA,
              Q2.5 = NA,
              Q25 = NA,
              Q75 = NA,
              Q97.5 = NA,
              .groups = "drop") |> 
    mutate(Name = candidate_names[ID])
  
  samples_beta_gamma_2 <- samples_beta_gamma |> 
    filter(P > 0) |>
    left_join(groups, by = "Name") |> 
    group_by(Chain, Iteration, Group_ID) |> 
    reframe(beta = sum(beta, na.rm = TRUE),
            gamma = 1*any(gamma == 1),
            Name = Group_Name)
  
  summ_beta_gamma_2 <- samples_beta_gamma_2 %>%
    group_by(Name) %>%
    mutate(beta = ifelse(gamma > 0, beta, NA)) |> 
    summarize(P_present = mean(gamma),
              P_absent = 1 - P_present,
              Mean = mean(beta, na.rm = TRUE),
              Median = median(beta, na.rm = TRUE),
              SD = sd(beta, na.rm = TRUE),
              Q2.5 = quantile(beta,.025,na.rm = TRUE),
              Q25 = quantile(beta,.25,na.rm = TRUE),
              Q75 = quantile(beta,.75,na.rm = TRUE),
              Q97.5 = quantile(beta,.975,na.rm = TRUE))
  
  ## Combine summary tables 
  summ_beta_gamma_2 |> 
    bind_rows(summ_beta_gamma_1) |> 
    select(-ID) |> 
    arrange(Name)
} 

mixtures_summ <- function(samp_df, 
                          design, 
                          groups,
                          min_prob = .0001){
  
  ## Retrieve candidate names from design matrix
  candidate_names <- colnames(design)[-1]
  
  ## Process MCMC output
  samp_df <- samp_df |> 
    filter(Parameter %in% c("beta","gamma")) |>
    pivot_wider(names_from = Parameter, values_from = Value) |> 
    filter(gamma > 0)
  
  ## Group results
  group_df <- samp_df |> 
    left_join(groups, by = "Name") |> 
    group_by(Chain, Iteration, Group_Name) |> 
    summarize(Name = paste(Name[gamma == 1], collapse = "/"),
              Abundance = sum(beta),
              .groups = "drop") 
  
  ## Identify mixtures
  group_df <- group_df |> 
    group_by(Name) |> 
    mutate(Group_ID = cur_group_id()) |> 
    ungroup() |> 
    group_by(Chain,Iteration) |> 
    mutate(Mixture_Raw = paste(sort(Group_ID), collapse = "/"))
  
  ## Compute occurrences and proportions for each mixture
  mixtures <- group_df |>  
    group_by(Mixture_Raw) |>
    summarize(n=n(),.groups = "drop") |> 
    mutate(p = n/sum(n)) |> 
    arrange(desc(p)) |> 
    rowid_to_column(var = "Mixture")
  
  ## Label mixtures in grouped output
  group_df <- group_df |> 
    full_join(mixtures, by = c("Mixture_Raw")) |> 
    select(-Mixture_Raw,-n)

  ## Compute summaries
  group_df |> 
    filter(p > min_prob) |> 
    group_by(Mixture, p, Name) |> 
    summarise(Mean = mean(Abundance),
              Median = median(Abundance),
              SD = sd(Abundance),
              Q2.5 = quantile(Abundance, .025),
              Q25 = quantile(Abundance, .25),
              Q75	= quantile(Abundance, .75),
              Q97.5 = quantile(Abundance, .975),
              .groups = "drop")
  

}

fitted_summ <- function(samp_df, counts){
  
  samples_mu <- samp_df %>%
    filter(Parameter == "mu") |> 
    select(Chain, Iteration, Interval = ID, Fitted = Value) |> 
    full_join(select(counts,Interval,Count), by = "Interval") %>%
    group_by(Chain, Iteration) %>%
    mutate(Residual = Count - Fitted,
           Pearson = Residual / sqrt(Fitted)) %>%
    ungroup()
  
  mu_summ <- samples_mu %>%
    group_by(Interval, Count) %>%
    summarize(MeanFitted = mean(Fitted),
              Q2.5Fitted = quantile(Fitted,.025),
              Q25Fitted = quantile(Fitted,.25),
              Q75Fitted = quantile(Fitted,.75),
              Q97.5Fitted = quantile(Fitted,.975),
              MeanResidual=mean(Residual),
              MeanPearson=mean(Pearson)) %>%
    ungroup()
  
  return(mu_summ)
}


