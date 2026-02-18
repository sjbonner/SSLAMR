# Data processing

Data processing

## Usage

``` r
sslamr_data(
  spectrum,
  candidates = NULL,
  isotope_data = NULL,
  min_abundance = 0.001,
  max_isotopes = Inf,
  skip_isotopes = 0,
  epsilon = 0.05,
  binning = TRUE,
  min_mass_charge = NULL,
  max_mass_charge = NULL,
  prior_par = NULL,
  prescreen = NULL,
  prescreen_prior = NULL,
  prescreen_weight = FALSE,
  rounding = "nearest",
  pattern_tol = 0.05,
  ran.seed = unclass(Sys.time()),
  isoinfo = NULL,
  verbose = FALSE
)
```

## Arguments

- spectrum:

  Data frame of observed peaks with columns: Mass/Charge (numeric) and
  Intensity (numeric).

- candidates:

  Data from of candidate molecules with columns: Name (character) and
  Formula (character).

- isotope_data:

- min_abundance:

- max_isotopes:

- skip_isotopes:

- epsilon:

  Resolution (default = .05)

- min_mass_charge:

- max_mass_charge:

- ran.seed:

  Random seed for rounding counts (integer)

- isoinfo:

- verbose:

  If true then print messages tracking steps (boolean, default = TRUE)

- round:

## Value

List with 5 components: candidates – an augmented data frame with the
candidate information, intervals – data frame with interval information,
spectrum – augmented data frame with peak data, design – design matrix,
counts – vector of summed counts for each bin
