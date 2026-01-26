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


