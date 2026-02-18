greedy_fit_single <- function(name, x, count){
  
  fit1 <- glm(count ~ x,
              family = poisson(link = "identity"),
              start = c(mean(count),0))
  
  # Extract output values
  # 1) Coefficients
  coeffs <- fit1$coefficients
  
  # 2) Deviance
  #deviance <- fit$nulldev * fit$dev.ratio
  deviance <- fit1$null.deviance - fit1$deviance
  
  # 3) P-value
  p_val <- pchisq(deviance, 1, lower.tail = FALSE)
  
  tibble(Abundance = coeffs[2], 
         Deviance = deviance, 
         `p-value` = p_val)
}

greedy_fit <- function(data,
                       xlsx_out = NULL, 
                       critical = .05){

  # Group candidates
  if(!is.null(data$groups)){
    design1 <- data$design |> 
      pivot_longer(-Interval,
                   names_to = "Name",
                   values_to = "Value") |> 
      left_join(data$groups, by = "Name") |> 
      group_by(Group_Name, Interval) |> 
      summarize(Value = mean(Value),
                .groups = "drop") |> 
      rename(Name = Group_Name)
  }
  else{
    design1 <- data$design |> 
      pivot_longer(-Interval,
                   names_to = "Name",
                   values_to = "Value")
  }
  
  # Extract positive values
  design1 <- design1 |> 
    filter(Value > 0) 
  
  # Extract observed counts
  Count <- data$counts |> 
    pull("Count")
  
  output <- tibble()
  
  continue <- TRUE

  while(continue){
     # Remove any candidates associated with all empty intervals
    design1 <- design1 |> 
      mutate(Counts = Count[Interval]) |> 
      group_by(Name) |> 
      mutate(Total = sum(Counts)) |> 
      filter(Total > 0) |> 
      select(-Total)
    
    # Fit models separately for each candidate
    fits <- design1 |> 
      mutate(Counts = Count[Interval]) |> 
      group_by(Name) |> 
      summarize(greedy_fit_single(first(Name), Value, Count[Interval]),
                .groups = "drop") |> 
      filter(Abundance > 0) |>
      arrange(`p-value`)

    if(fits$`p-value`[1] < critical){
      
      # Save fit information
      output <- output |> 
        bind_rows(fits[1,])
      
      # Subtract from counts
      name <- fits |> 
        head(1) |> 
        pull("Name")
      
      abund <- fits |> 
        head(1) |> 
        pull("Abundance")
      
      tmp <- design1 |> 
        filter(Name == name)
      
      Count[tmp$Interval] <- (Count[tmp$Interval] - tmp$Value * abund) %>%
        round() %>%
        pmax(0)

      # Update design matrix  
      design1 <- design1 |> 
        filter(Name != name)
    }
    else{
      continue <- FALSE
    }
  }
    
  # Save output
  write_xlsx(output, 
             path = xlsx_out)
}

greedy_fit_2 <- function(spectrum = NULL,
                         candidates = NULL,
                         isotope_data = NULL,
                         adducts = NULL,
                         isoinfo = NULL,
                         replace_isoinfo = FALSE,
                         min_abundance = .001,
                         group_formula = NULL,
                         group_pattern = TRUE,
                         pattern_tol = .05,
                         binning = TRUE,
                         epsilon = .05,
                         rounding = "nearest",
                         min_mass_charge = NULL,
                         max_mass_charge = NULL,
                         prescreen = 0,
                         prescreen_prior = 1,
                         prescreen_weight = FALSE,
                         metric = c("p-value","abundance"),
                         critical = .05,
                         xlsx_out = NULL,
                         verbose = TRUE){
  
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
  
  # Determine ranking metric
  metric <- match.arg(metric)
  
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
  
  # Define design matrix
  data <- sslamr_data(spectrum, 
                      candidates = candidates,
                      isotope_data = isotope_data,
                      min_abundance = min_abundance,
                      max_isotopes = Inf,
                      skip_isotopes = 0,
                      binning = binning,
                      epsilon = epsilon, 
                      min_mass_charge = min_mass_charge,
                      max_mass_charge = max_mass_charge,
                      prior_par = NULL,
                      prescreen = prescreen,
                      prescreen_prior = prescreen_prior,
                      prescreen_weight = prescreen_weight,
                      rounding = rounding,
                      pattern_tol = pattern_tol,
                      isoinfo = isoinfo,
                      verbose = verbose)
  
  # Convert design matrix to data frame
  design <- data$design |> 
    pivot_longer(-Interval,
                 names_to = "Name",
                 values_to = "Proportion")
  
  # Extract rows with positive isotope proportions
  design <- design |> 
    filter(Proportion > 0) 
  
  # Add parent information
  design <- design |> 
    left_join(data$candidates |> select(Name = ID, Parent) |> distinct(), by = "Name")
  
  # Extract observed counts
  Count <- data$counts |> 
    pull("Count")
  
  # Initialize output table
  results <- tibble()
  
  continue <- TRUE
  
  if(verbose){
    message("Fitting model ...\n")
    pb <- txtProgressBar(max = length(unique(design$Name)), style = 1)
    k <- 1
  }
  
  while(continue){
    if(verbose)
      setTxtProgressBar(pb, k)
    
    # Remove candidates whose parents are not yet in the model
    if(nrow(results) == 0){
      design1 <- design |> 
        filter(Parent == Name)
    }
    else{
      design1 <- design |> 
        mutate(Keep = (Parent == Name) + (Parent %in% results$Name)) |>
        filter(Keep > 0)
    }
   
    # Remove any candidates associated with all empty intervals
    design1 <- design1 |> 
      mutate(Counts = Count[Interval]) |> 
      group_by(Name) |> 
      mutate(Total = sum(Counts)) |> 
      filter(Total > 0) |> 
      select(-Total)
    
    # Fit models separately for each candidate
    fits <- design1 |> 
      mutate(Counts = Count[Interval]) |> 
      group_by(Name) |> 
      mutate(Total = sum(Count)) |> 
      summarize(greedy_fit_single(first(Name), Proportion, Count[Interval]),
                .groups = "drop") |> 
      filter(Abundance > 0, `p-value` < critical)
    
    # Arrange by chosen metric
    if(metric == "p-value"){
      fits <- fits |> 
        arrange(`p-value`)
    }
    else if(metric == "abundance"){
      fits <- fits |> 
        arrange(desc(Abundance))
    }
    
    if(nrow(fits) > 0){
      
      # Save fit information
      results <- results |> 
        bind_rows(fits[1,])
      
      # Subtract from counts
      name <- fits |> 
        head(1) |> 
        pull("Name")
      
      abund <- fits |> 
        head(1) |> 
        pull("Abundance")
      
      tmp <- design1 |> 
        filter(Name == name)
      
      Count[tmp$Interval] <- (Count[tmp$Interval] - tmp$Proportion * abund) %>%
        round() %>%
        pmax(0)
      
      # Update design matrix  
      design <- design |> 
        filter(Name != name)
      
      if(verbose)
        k <- k + 1
    }
    else{
      continue <- FALSE
    }
  }
  
  output <- list(original = results)
  
  if(group_pattern){
    
    # Identify pattern groups
    pattern_groups <- data$design |> 
      group_by_pattern()
    
    # Restrict to candidates in results
    test <- results |> 
      left_join(pattern_groups, by = "Name") |> 
      select(Group_ID, Group_Name) |> 
      arrange(Group_ID) |> 
      unique()
    
    # Split group names to identify potential candidates
    test1 <- test |>
      separate_longer_delim(Group_Name, delim = "/") |> 
      rename(Name = Group_Name) |> 
      left_join(select(candidates, Name, Parent), by = "Name")
    
    # Keep only candidates whose parents are also in results
    test2 <- test1 |> 
      mutate(Keep = (Parent %in% pull(test1,Name))) |> 
      filter(Keep) |> 
      select(Group_ID, Name)
    
    # Rebuild group names
    test3 <- test2 |> 
      group_by(Group_ID) |> 
      arrange(Name) |> 
      mutate(Group_Name = paste(Name, collapse = "/")) |> 
      ungroup() |> 
      select(Name, Group_Name)
    
    # Compile group output
    group_results <- results |> 
      left_join(test3, by = "Name") |> 
      select(Group_Name, Abundance) |> 
      group_by(Group_Name) |>
      summarize(Abundance = sum(Abundance),
                .groups = "drop")
    
    output$grouped <- group_results
  }  
  
  # Save output
  write_xlsx(output,
             path = xlsx_out)
  
  return(list(data = list(design),
              output = output))
}