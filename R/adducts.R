formula_to_tibble <- function(data, name = "Formula"){
  
  data %>% 
    group_by_all() %>%
    ## Split formula into element/counts pairs
    summarize(split = unlist(str_extract_all(Formula, "[A-Z][a-z]*[0-9]*+")), .groups = "drop") %>%
    ## Extract element and count from pairs
    mutate(element = str_extract(split,"[A-Z][a-z]?+"),
           count = as.numeric(str_extract(split, "[0-9]*$"))) %>%
    ## Replace missing counts with 1
    replace_na(list(count = 1)) %>%
    select(-split)
}

modify_adducts <- function(candidates, adducts){
  ## Convert candidate formulas to tibble
  candidates1 <- candidates %>%
    select(-Charge) %>%
    formula_to_tibble()
  
  ## Convert adduct formulas to tibble
  adducts1 <- adducts %>%
    select(-Action) %>%
    formula_to_tibble()
  
  ## Identify unique elements among all candidates and adducts
  all_element <- tibble(element = unique(c(unique(candidates1$element), 
                                         unique(adducts1$element))))
  
  ## Complete tibbles with all elements 
  candidates2 <- candidates %>% 
    crossing(all_element) %>%
    full_join(candidates1, by = c("Name","Formula","element")) %>%
    replace_na(list(count = 0))
  
  adducts2 <- adducts %>% 
    crossing(all_element) %>%
    full_join(adducts1, by = c("Name","Formula","element")) %>%
    replace_na(list(count = 0))
  
  ## Join tibbles
  output <- full_join(candidates2,adducts2, by = "element", suffix = c(".cand",".adduct")) %>%
    mutate(Parent = Name.cand,
           Name = paste(Name.cand,Action,Name.adduct),
           count = count.cand + 
             (((Action == "+") - (Action == "-")) * count.adduct)) %>%
    select(Parent, Name, Charge, element,count) %>%
    filter(count != 0)
  
  ## Remove any impossibilities
  output <- output %>%
    group_by(Name, Charge) %>%
    mutate(valid = all(count > 0)) %>%
    filter(valid)
  
  ## Reconstitute formulas
  output <- output %>%
    group_by(Parent, Name, Charge) %>%
    summarize(Formula = paste0(c(rbind(element,count)),collapse = ""))
  
  ## Append to list of candidate
  candidates %>%
    mutate(Parent = Name) |> 
    bind_rows(output)
}

