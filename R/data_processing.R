
# This function takes the name and lipid attained from
# get_PC_formula gets its isotope distribution from ecipex
# Could use improvement
get_iso_info <- function(name, 
                         formula, 
                         charge, 
                         min_abundance = -Inf, 
                         max_isotopes = Inf,
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
    filter(Isotope <= max_isotopes) %>%
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
                           isoinfo = NULL){
  candidates %>%
    rowwise() %>%
    do(get_iso_info(.$Name, .$Formula, .$Charge, 
                    min_abundance = min_abundance, 
                    max_isotopes = max_isotopes,
                    isoinfo = isoinfo)) %>%
    ungroup() %>%
    mutate(Mass = (Mass - Charge * 0.000548579909)/abs(Charge))
}

# lipid bins: defines the boundaries of the bins for creating
# the design matrix
create_bins <- function(isotope_df, epsilon = .05){
  isotope_df %>%
    arrange(Mass) %>%
    mutate(diff = Mass-lag(Mass),
           step = ifelse(is.na(diff), 0, diff > epsilon),
           Group = cumsum(step) + 1) %>%
    group_by(Group) %>%
    mutate(GroupMass = round(mean(Mass)/epsilon)*epsilon,
           Lower = min(Mass) - epsilon/2,
           Upper = max(Mass) + epsilon/2) %>%
    dplyr::select(-diff, -step) %>%
    ungroup() %>%
    arrange(ID,Isotope)
}


# this function takes the output from lipids_info and 
# creates the design matrix
build_design <- function(lip_df, epsilon=0.05){
  
  design <- lip_df %>%
    complete(ID,nesting(GroupMass,Group),fill = list(Isotope = NA, Abund = 0)) %>%
    arrange(ID) %>%
    dplyr::select(-Isotope,-Mass) %>%
    group_by(ID,Group) %>%
    summarize(Abund = sum(Abund)) %>%
    ungroup() %>%
    pivot_wider(names_from = ID, values_from = Abund)
  
  return(design)
}

peaks_to_bins <- function(peaks, bins, MC = "Mass", lower = "Lower", upper = "Upper"){
  # Assign peaks to bins
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
                        epsilon=0.05, 
                        ran.seed=unclass(Sys.time()),
                        isoinfo = NULL,
                        verbose = FALSE){
  
  # Process candidate molecules
  
  # 1) Retrieve candidate isotope patterns
  if(is.null(isotope_data)) {
    if(verbose)
      message("  Retrieving candidate isotope data...")
    
    isotope_data <- candidate_info(candidates, 
                                   min_abundance = min_abundance,
                                   max_isotopes = max_isotopes,
                                   isoinfo = isoinfo)
  }
  
  # 2) Assign candidates to bins
  if(verbose)
    message("  Binning candidate isotopes...")
  candidate_bins <- create_bins(isotope_data, epsilon)
  
  # 3) Summarize bin information
  bins <- candidate_bins %>%
    group_by(Group) %>%
    summarize(Isotopes = n(), 
              Mass = GroupMass[1],
              Lower = Lower[1],
              Upper = Upper[1],
              Width = Upper[1] - Lower[1]) %>%
    ungroup()
  
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
    mutate(Count = Intensity * ccf) %>% # Convert intensity to count
    peaks_to_bins(bins, "Mass/Charge")
  
  # 2) Summarize counts within bins
  counts <- summarize_counts(MS, bins, ran.seed = ran.seed)
  
  data <- bins %>% 
    full_join(counts, by = "Group") %>%
    full_join(design, by = "Group") %>%
    arrange(Group)
  
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


