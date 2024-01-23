greedy_fit <- function(results, xlsx_out = NULL){
  
  # Extract mass ordered list of candidates
  candidates <- results$data$candidates %>%
    group_by(ID) %>%
    mutate(min_mass = min(Mass),
           row_num = row_number()) %>%
    filter(row_num == 1) %>%
    ungroup() %>%
    arrange(min_mass, ID) %>%
    pull("ID")
  
  
  # Extract observed counts
  Count <- results$data$data$Count
  coeffs <- matrix(nrow = length(candidates), ncol = 4)
  
  # Initialize progress bar
  pb <- txtProgressBar(0, length(candidates), style = 3)
  
  for(k in 1:length(candidates)){
    
    # Update progress bar
    setTxtProgressBar(pb, k)
    
    # Identify non-zero rows of design matrix for next candidate
    ind <- which(results$data$data[,candidates[k]] > 0)

    # Extract entries of design matrix    
    x <- results$data$data[ind,] %>%
      pull(candidates[k])
    
    X <- cbind(1,x)
    
    if(any(Count[ind] > 0)){
      # Fit Poisson model with non-negative constraint
      fit <- glmnet(X, 
                    Count[ind],
                    family = poisson(link = "identity"), 
                    lower.limits = 0, 
                    lambda = 0)
      
      # Extract output values
      # 1) Coefficients
      coeffs[k,1:2] <- fit$beta %>%
        as.vector()
      
      # 2) Deviance
      coeffs[k, 3] <- fit$nulldev * fit$dev.ratio
      
      # 3) P-value
      coeffs[k, 4] <- pchisq(coeffs[k, 3], 1, lower.tail = FALSE)
      
      if(coeffs[k,2] > 0 & coeffs[k,4] < .05){
        # Remove candidate contribution
        Count[ind] <- (Count[ind] - x * coeffs[k,2]) %>%
          round() %>%
          pmax(0)
      }
    }
    else{
      coeffs[k,1:2] <- 0
    }
  }
  
  # Combine results
  results <- tibble(Candidate = candidates, 
                    Beta = coeffs[,2],
                    Chi_square = coeffs[,3],
                    p_value = coeffs[,4])
  
  if(!is.null(xlsx_out))
    write_xlsx(results, xlsx_out)
  
  results
}