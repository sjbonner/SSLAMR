get_iso_info <- function(name, 
                         formula, 
                         charge, 
                         min_abundance = -Inf, 
                         max_isotopes = Inf,
                         skip_isotopes = 0,
                         isoinfo = NULL){
  
  if(is.null(isoinfo)){
    ## Retrieve information from ecipex
    dist <- ecipex(formula, sortby = "mass",
                   groupby = "mass", gross = TRUE)[[1]][,-3] 
  }
  else{
    dist <- ecipex(formula, sortby = "mass",isoinfo = isoinfo,
                   groupby = "mass", gross = TRUE)[[1]][,-3]
  }
  ## Format table
  dist %>% 
    as_tibble() %>% 
    filter(abundance > min_abundance) %>%
    rowid_to_column(var = "Isotope") %>%
    filter(Isotope > skip_isotopes, Isotope <= max_isotopes + skip_isotopes) %>%
    add_column(ID=name, Charge = charge, .before=1) %>%
    rename(Mass = centroidMass,
           Abund = abundance)
}

# takes the grouped counts from the previous function and randomly
# rounds it based on the remainder. this accounts for the uncertainty
# of rounding counts. intended for replication
round_bernoulli <- function(counts, ran.seed = unclass(Sys.time())){
  set.seed(ran.seed)
  probs <- counts %% 1
  decision <- sapply(probs, function(x){rbinom(1, size=1, prob=x)})
  new_counts <- counts - probs + decision
  return(new_counts)
}

# this function gets the lipids info from a data frame containing
# name in the first column and formula in the second column.

candidate_info <- function(candidates, 
                           min_abundance = .001, 
                           max_isotopes = Inf,
                           skip_isotopes = 0,
                           min_mass_charge = -Inf,
                           max_mass_charge = Inf,
                           isoinfo = NULL){
  candidates %>%
    rowwise() %>%
    do(get_iso_info(.$Name, .$Formula, .$Charge, 
                    min_abundance = min_abundance, 
                    max_isotopes = max_isotopes,
                    skip_isotopes = skip_isotopes,
                    isoinfo = isoinfo)) %>%
    ungroup() %>%
    mutate(Mass = (Mass - Charge * 0.000548579909)/abs(Charge)) %>%
    filter(Mass >= min_mass_charge, Mass <= max_mass_charge)
}

# lipid bins: defines the boundaries of the bins for creating
# the design matrix
create_bins <- function(isotope_df, epsilon = .05, min_mass_charge = 0, max_mass_charge = 1000){
  
  bins1 <- tibble(Mass = sort(unique(pull(isotope_df,"Mass")))) %>%
    mutate(Delta = Mass - lag(Mass, default = 0),
           Increment = Delta > 2 * epsilon,
           Group = cumsum(Increment)) %>%
    group_by(Group) %>%
    summarize(Isotopes = n(),
              Mean = mean(Mass),
              Lower = min(Mass) - epsilon,
              Upper = max(Mass) + epsilon,
              .groups = "drop") %>%
    select(-Group)
  
  bins2 <- tibble(Lower = c(min_mass_charge,pull(bins1,"Upper")),
                  Upper = c(pull(bins1, Lower), max_mass_charge),
                  Isotopes = 0,
                  Mean = NA)
  
  bind_rows(bins1,bins2) %>%
    arrange(Lower) %>%
    rowid_to_column("Group") %>%
    mutate(Width = Upper - Lower)
}

# this function takes the output from lipids_info and 
# creates the design matrix
build_design <- function(candidate_bins, epsilon=0.05){
  
  design <- candidate_bins %>%
    complete(ID,Group,fill = list(Isotope = NA, Abund = 0)) %>%
    arrange(ID) %>%
    dplyr::select(-Isotope,-Mass) %>%
    group_by(ID,Group) %>%
    summarize(Abund = sum(Abund)) %>%
    ungroup() %>%
    pivot_wider(names_from = ID, values_from = Abund)
  
  return(design)
}

assign_to_bins <- function(peaks, bins, MC = "Mass", lower = "Lower", upper = "Upper"){
  # Assign quantities to bins
  peaks %>%
    mutate(Group = findInterval(.data[[MC]], sort(c(-Inf,bins[[lower]],bins[[upper]],Inf))),
           Group = ifelse(Group %% 2 == 0, Group/2, NA))
}

summarize_counts <- function(MS, bins, epsilon=0.05, ran.seed=unclass(Sys.time())){
  
  # Round counts and sum within bins
  MS <- MS %>%
    mutate(Count = round_bernoulli(Count, ran.seed = ran.seed)) %>%
    group_by(Group) %>%
    summarize(Peaks = n(),
              Count = sum(Count))
  
  # Complete table with empty bins
  bins %>%
    select(Group) %>%
    left_join(MS, by = "Group") %>%
    replace_na(list(Peaks = 0, Count = 0))
}

#' Data processing
#'
#' @param MS Data frame of observed peaks with columns: Mass/Charge (numeric) and Intensity (numeric).
#' @param candidates Data from of candidate molecules with columns: Name (character) and Formula (character).
#' @param epsilon Resolution (default = .05)
#' @param ran.seed Random seed for rounding counts (integer)
#' @param verbose If true then print messages tracking steps (boolean, default = TRUE)
#' @param min_abundance 
#' @param isoinfo 
#'
#' @return List with 5 components: candidates -- an augmented data frame with the candidate information, 
#' bins -- data frame with bin information, MS -- augmented data frame with peak data, design -- design matrix, counts -- vector of summed counts for each bin 
#' @export
#'
sslamr_data <- function(MS, 
                        candidates = NULL,
                        isotope_data = NULL,
                        min_abundance = .001,
                        max_isotopes = Inf,
                        skip_isotopes = 0,
                        epsilon=0.05, 
                        min_mass_charge = NULL,
                        max_mass_charge = NULL,
                        ran.seed=unclass(Sys.time()),
                        isoinfo = NULL,
                        verbose = FALSE){
  
  # Set default range for mass/charge
  if(is.null(min_mass_charge)){
    min_mass_charge <- 0
  }
  if(is.null(max_mass_charge)){
    max_mass_charge <- max(MS$`Mass/Charge` + 2 * epsilon)
  }
  
  # Process candidate molecules
  
  # 1) Retrieve candidate isotope patterns
  if(is.null(isotope_data)) {
    if(verbose)
      message("  Retrieving candidate isotope data...")
    
    isotope_data <- candidate_info(candidates, 
                                   min_abundance = min_abundance,
                                   max_isotopes = max_isotopes,
                                   skip_isotopes = skip_isotopes,
                                   min_mass_charge = min_mass_charge,
                                   max_mass_charge = max_mass_charge,
                                   isoinfo = isoinfo)
  }
  
  # 2) Define bins
  if(verbose)
    message("  Defining bins...")
  bins <- create_bins(isotope_data, 
                      epsilon, 
                      min_mass_charge = min_mass_charge,
                      max_mass_charge = max_mass_charge)
  
  # 3) Assign isotopes to bins
  candidate_bins <- isotope_data %>%
    assign_to_bins(bins,"Mass")
  
  # 4) Construct design matrix 
  if(verbose)
    message("  Constructing design matrix...")
  design <- build_design(candidate_bins, epsilon=epsilon)
  
  # Process spectrum
  
  # 1) Assign peaks to bins
  if(verbose)
    message("  Assigning peaks to bins...")
  
  ccf <- 1 / (MS %>% filter(Intensity > 0) %>% pull("Intensity") %>% min())
  
  MS <- MS %>% 
    filter(Intensity > 0) %>% # Remove peaks with zero intensity
    filter(`Mass/Charge` >= min_mass_charge, `Mass/Charge` <= max_mass_charge) %>% # Restrict to analysis window
    mutate(Count = Intensity * ccf) %>% # Convert intensity to count
    assign_to_bins(bins, "Mass/Charge")
  
  # 2) Summarize counts within bins
  counts <- summarize_counts(MS, bins, ran.seed = ran.seed)
  
  data <- bins %>% 
    full_join(counts, by = "Group") %>%
    full_join(design, by = "Group") %>%
    arrange(Group) %>%
    mutate(across(!c(Group, Isotopes, Mean, Lower, Upper, Width, Peaks, Count), ~replace_na(.x, 0)))
  
  # Return complete output
  list(candidates = candidate_bins,
       bins = bins,
       MS = MS,
       data = data)
}

prescreen <- function(candidates, MS, prescreen){
  if(verbose)
    message("  Prescreening...")
  
  nin <- nrow(candidates)
  
  # Compute average count in bins with positive abundance for each candidate
  tmp <- design %>% 
    pivot_longer(-Group,names_to = "Name",values_to = "Abundance") %>% 
    left_join(counts,by = "Group") %>% 
    group_by(Name) %>% 
    mutate(Count_Avg = mean(Count[Abundance > 0])) 
  
  # Retain candidates with average counts greater than prescreen
  design <- tmp %>% 
    filter(Count_Avg >= prescreen) %>%
    select(Group, Name, Abundance) %>%
    pivot_wider(names_from = Name, values_from = Abundance)
  
  nout <- ncol(design) - 1
  
  if(verbose)
    message("    ",nout," of ", nin, " candidates retained.")
}


