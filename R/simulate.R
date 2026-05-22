#' Mass Spectrum Simulation
#'
#' @param candidates List of candidate compounds included in the mixture. 
#' @param isoinfo Custom element isotope information that either replaces or
#'   adds to the definitions in [ecipex()].
#' @param mc_cal_error Mass/charge calibration error. This error is assumed to
#'   be systematic and is added directly to the mass/charge values for all peaks.
#' @param mc_sd Standard deviation of the mass/charge observation error. 
#' @param mc_min_diff Minimum discernable mass/charge difference. Peaks separated
#'   by less than this value are combined.
#' @param intercept Mean count of the noise per amu.
#' @param mc_min Lower bound of the mass/charge window.
#' @param mc_max Upper bound of the mass/charge window.
#' @param vif Variance inflation factor, the ratio of the mean abundance to the
#'   variance. If `vif=1`, the default, then counts are generate from Poisson distributions. 
#'   Otherwise, counts are generated from negative binomials with vif defining
#'   the ratio of the mean and variance. 
#' @param saturation If finite, then counts are truncated to be no greater than
#'   this value to mimic saturation of the sensor. Defaults to `Inf` so that
#'   no saturation occurs.
#'
#' @returns A tibble with the simulated mass spectrum containing columns for 
#'    the mass/charge ratio and intensity. 
#' @export
#'
#' @examples
#' # TBD
sslamr_simulate <- function(candidates, 
                            isoinfo = NULL,
                            mc_cal_error,
                            mc_sd,
                            mc_min_diff,
                            intercept,
                            mc_min,
                            mc_max,
                            vif = 1,
                            saturation = Inf){
  
  # Check arguments
  if(vif < 1)
    stop("The variance inflation factor (vif) must be at least 1.\n")
  
  if(saturation < 1)
    stop("The maximum count (saturation) must be at least 1.\n")
  
  # Generate isotope information
  isotopes <- candidate_info(candidates, isoinfo = isoinfo)
  
  # Generate counts
  if(vif > 1){
    # Simulate from negative binomial
    spectrum <- candidates %>%
      select(ID = Name, Abundance) %>%
      full_join(isotopes, by = "ID") %>%
      group_by(ID) %>%
      mutate(mu = Abundance * Abund/sum(Abund),
             beta = vif - 1,
             lambda = rgamma(n(), shape = mu/beta, scale = beta),
             Count = rpois(n(), lambda)) %>%
      ungroup() %>%
      filter(Count > 0)
  }
  else{
    # Simulate from Poisson
    spectrum <- candidates %>%
      select(ID = Name, Abundance) %>%
      full_join(isotopes, by = "ID") %>%
      group_by(ID) %>%
      mutate(mu = Abundance * Abund/sum(Abund),
             Count = rpois(n(), mu)) %>%
      ungroup() %>%
      filter(Count > 0)
  }  
  # Add MC error
  spectrum <- spectrum %>%
    mutate(Mass = Mass + mc_cal_error + rnorm(n(), 0, mc_sd))
  
  # Add noise
  n_val <- rpois(1, (mc_max - mc_min) * intercept)
  
  spectrum <- spectrum %>%
    bind_rows(tibble(ID = paste0("Dummy",as.character(1:n_val)),
                     Mass = sort(runif(n_val, mc_min, mc_max)),
                     Count = 1,
                     Isotope = 1))
  
  # Group nearby observations
  bins <- spectrum %>%
    create_bins_1(epsilon = mc_min_diff) 
  
  spectrum <- spectrum %>%
    assign_to_bins(bins) %>%
    group_by(Interval) %>%
    summarize(Mass = round(sum(Mass * Count)/sum(Count),5),
              Count = pmin(sum(Count), saturation)) %>%
    arrange(Mass)
  
  # Return spectrum
  spectrum %>%
    select(`Mass/Charge` = Mass, Intensity = Count) %>%
    arrange(`Mass/Charge`)
  
}