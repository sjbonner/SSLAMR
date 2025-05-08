fit_single <- function(name, x, count){
  
  #browser()
  
  fit <- glmnet(cbind(1,x),
                count,
                family = poisson(link = "identity"),
                #intercept = FALSE,
                lower.limits = 0, 
                lambda = 0)
  
  fit1 <- glm(count ~ x,
              family = poisson(link = "identity"),
              start = c(mean(count),0))
  
  # Extract output values
  # 1) Coefficients
  #coeffs <- fit$beta %>%
  #  as.vector()
  
  coeffs <- fit1$coefficients
  
  # 2) Deviance
  #deviance <- fit$nulldev * fit$dev.ratio
  deviance <- fit1$null.deviance - fit1$deviance
  
  # 3) P-value
  p_val <- pchisq(deviance, 1, lower.tail = FALSE)
  
  tibble(abundance = coeffs[2], 
         deviance = deviance, 
         p = p_val)
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
      mutate(Total = sum(Count)) |> 
      summarize(fit_single(first(Name), Value, Count[Interval]),
                .groups = "drop") |> 
      filter(p < critical, abundance > 0) |>
      arrange(p)

    if(any(fits$p < critical)){
      
      # Save fit information
      output <- output |> 
        bind_rows(fits[1,])
      
      # Subtract from counts
      name <- fits |> 
        head(1) |> 
        pull("Name")
      abund <- fits |> 
        head(1) |> 
        pull("abundance")
      
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
    
  # Return
  output
}