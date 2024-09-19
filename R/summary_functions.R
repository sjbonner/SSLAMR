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
    select(Chain, Iteration, beta_tmp_sd)
  
  summ_beta_sd <- samples_beta_sd %>%
    summarize(Mean = mean(beta_tmp_sd, na.rm = TRUE),
              Median = median(beta_tmp_sd, na.rm = TRUE),
              SD = sd(beta_tmp_sd, na.rm = TRUE),
              Q2.5 = quantile(beta_tmp_sd,.025,na.rm = TRUE),
              Q25 = quantile(beta_tmp_sd,.25,na.rm = TRUE),
              Q75 = quantile(beta_tmp_sd,.75,na.rm = TRUE),
              Q97.5 = quantile(beta_tmp_sd,.975,na.rm = TRUE))
  
  return(summ_beta_sd)
}

intercept_summ <- function(samp_df){
  samples_beta0 <- samp_df %>%
    select(Chain, Iteration, beta0)
  
  summ_beta0 <- samples_beta0 %>%
    summarize(Mean = mean(beta0, na.rm = TRUE),
              Median = median(beta0, na.rm = TRUE),
              SD = sd(beta0, na.rm = TRUE),
              Q2.5 = quantile(beta0,.025,na.rm = TRUE),
              Q25 = quantile(beta0,.25,na.rm = TRUE),
              Q75 = quantile(beta0,.75,na.rm = TRUE),
              Q97.5 = quantile(beta0,.975,na.rm = TRUE))
  
  return(summ_beta0)
}

beta.gamma_summ <- function(samp_df, design){
  candidate_names <- colnames(design)[-1]
    
  samples_beta_gamma <- samp_df %>%
    select(Chain, Iteration, starts_with("beta["),starts_with("gamma["))%>%
    pivot_longer(c(starts_with("beta"),starts_with("gamma")),
                 names_to = "Parameter",values_to = "Value") %>%
    mutate(Name = as.integer(str_extract(Parameter,"[0-9]+")),
           Parameter = str_extract(Parameter,"[a-z]+")) %>%
    pivot_wider(names_from = Parameter,values_from = Value) %>%
    mutate(beta = ifelse(gamma, beta, NA),
           Name = candidate_names[Name])
  
  summ_beta_gamma <- samples_beta_gamma %>%
    group_by(Name) %>%
    summarize(P_present = mean(gamma),
              P_absent = 1 - P_present,
              Mean = mean(beta, na.rm = TRUE),
              Median = median(beta, na.rm = TRUE),
              SD = sd(beta, na.rm = TRUE),
              Q2.5 = quantile(beta,.025,na.rm = TRUE),
              Q25 = quantile(beta,.25,na.rm = TRUE),
              Q75 = quantile(beta,.75,na.rm = TRUE),
              Q97.5 = quantile(beta,.975,na.rm = TRUE))
  
  return(summ_beta_gamma)
} 

mixtures_summ <- function(samp_df, design){
  
  ## Retrieve candidate names from design matrix
  candidate_names <- colnames(design)[-1]
  
  ## Extract presence indicators
  gamma <- samp_df |> 
    select(starts_with("gamma")) |>
    rowid_to_column(var = "Iteration") |> 
    pivot_longer(starts_with("gamma"),names_to = "Parameter",values_to = "Presence")|> 
    mutate(Name = candidate_names[as.integer(str_extract(Parameter,"[0-9]+"))]) |> 
    select(-Parameter) 
  
  ## Extract abundance values
  beta <- samp_df |> 
    select(starts_with("beta[")) |>
    rowid_to_column(var = "Iteration") |> 
    pivot_longer(starts_with("beta["),names_to = "Parameter",values_to = "Abundance")|> 
    mutate(Name = candidate_names[as.integer(str_extract(Parameter,"[0-9]+"))]) |> 
    select(-Parameter) 
  
  ## Label mixtures using dplyr's default sorting
  mixtures1 <- gamma |> 
    pivot_wider(names_from = Name, values_from = Presence) |> 
    group_by(across(-c(Iteration))) |> 
    mutate(Label = cur_group_id()) |> 
    ungroup() |> 
    select(Iteration, Label)
  
  gamma <- mixtures1 |> 
    full_join(gamma, by = "Iteration")
  
  ## Compute occurrences and proportions for each unique mixture
  mixtures1 <- mixtures1 |>  
    group_by(across(-c(Iteration))) |>
    summarize(n=n(),.groups = "drop") |> 
    mutate(p = n/sum(n))
  
  ## Relabel in descending order of occurrence
  mixtures1 <- mixtures1 |>
    arrange(desc(p)) |> 
    rowid_to_column(var = "Mixture")
  
  ## Combine proportions, presence, and abundance
  mixtures2 <- mixtures1 |> 
    full_join(gamma, by = "Label") |> 
    left_join(beta, by = c("Iteration","Name")) |> 
    select(-Iteration) |> 
    select(-Label)
  
  ## Compute summaries
  mixtures2 |> 
    filter(Presence > 0) |> 
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
    select(Chain, Iteration, starts_with("mu")) %>%
    pivot_longer(starts_with("mu"),names_to = "Parameter", values_to = "Fitted") %>%
    mutate(Interval = as.integer(str_extract(Parameter,"[0-9]+"))) %>%
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

candidate_summary <- function(scan.obj, candidate){
  lipsum <- filter(scan.obj$bg.summ, Name==candidate) %>%
    mutate(Time=scan.obj$scan.time)
  
  return(lipsum)
}

mass_summary <- function(scan.obj, massID){
  lipsum <- filter(scan.obj$fit.summ, MassID==massID) %>%
    mutate(Time=scan.obj$scan.time)
  
  return(lipsum)
}

r_summ <- function(samp_df){
  samples_r <- samp_df %>%
    select(Chain, Iteration, r)
  
  summ_r <- samples_r %>%
    summarize(Mean = mean(r, na.rm = TRUE),
              Median = median(r, na.rm = TRUE),
              SD = sd(r, na.rm = TRUE),
              Q2.5 = quantile(r,.025,na.rm = TRUE),
              Q25 = quantile(r,.25,na.rm = TRUE),
              Q75 = quantile(r,.75,na.rm = TRUE),
              Q97.5 = quantile(r,.975,na.rm = TRUE))
  
  return(summ_r)
}

intercept_summary <- function(scan.obj){
  return(data.frame(scan.obj$int.summ, Time=scan.obj$scan.time))
}

