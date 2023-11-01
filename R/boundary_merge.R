# Load packages
library(tidyverse)

# Define boundaries
n_bins <- 10

bounds <- sort(100 * runif(2 *n_bins))

bins <- tibble(ID=1:n_bins,
               Lower = bounds[2*(1:n_bins) - 1],
               Upper = bounds[2*(1:n_bins)])

# Simulate data
m <- 100

peaks <- tibble(ID = 1:m,
                Mass = sort(runif(m, 0,100)))

# Identify bins
peaks1 <- peaks %>%
  mutate(Lower_cut = findInterval(Mass, c(bins$Lower,Inf)),
         Upper_cut = findInterval(Mass, c(-Inf, bins$Upper)),
         Bin = ifelse(Lower_cut == Upper_cut, Lower_cut, NA))


# Identify bins 2
peaks2 <- peaks %>%
  mutate(Bin = findInterval(Mass, sort(c(-Inf,bins$Lower,bins$Upper,Inf))),
         Bin = ifelse(Bin %% 2 == 0, Bin/2, NA))

