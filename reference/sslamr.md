# Spike-and-Slab Analysis for Mass Spectrometry in R

Implements the Bayesian Poisson generalized linear model of Bonner,
Bonner, Walker and Cao for deconvolution of isotope patterns in spectra
obtained from mass spectrometry using the spike-and-slab prior to
simultaneously identify and quantify compounds of interest within a
sample.

## Usage

``` r
sslamr(
  spectrum = NULL,
  candidates = NULL,
  isotope_data = NULL,
  adducts = NULL,
  isoinfo = NULL,
  replace_isoinfo = FALSE,
  min_abundance = 0.001,
  max_isotopes = Inf,
  skip_isotopes = 0,
  group_formula = NULL,
  group_pattern = TRUE,
  pattern_tol = 0.05,
  binning = TRUE,
  epsilon = 0.05,
  min_mass_charge = NULL,
  max_mass_charge = NULL,
  prescreen = 0,
  prescreen_prior = 1,
  prescreen_weight = FALSE,
  rounding = "nearest",
  n.chains = 3,
  n.adapt = 1000,
  n.burnin = 1000,
  n.sampling = 10000,
  n.thin = 1,
  prior_par = NULL,
  model = "hierarchical",
  inits_seed = NULL,
  run_model = TRUE,
  mixtures = FALSE,
  reload_results = FALSE,
  xlsx_out = NULL,
  verbose = TRUE
)
```

## Arguments

- spectrum:

  Observed spectrum. May be a string specifying the path to the data
  file or an existing data frame. See Input Details for further
  information. (character or dataframe)

- candidates:

  List of candidate compounds. May be a string specifying the path to
  the data file or an existing data frame. See Input Details for further
  information. (character or dataframe)

- isotope_data:

  Preprocessed candidate isotope information. May be a string specifying
  the path to the data file or an existing data frame. See Isotope Data
  for further information. (character or dataframe)

- adducts:

  List of adducts. May be a string specifying the path to the data file
  or an existing data frame. See Input Details for further information.
  (character or dataframe)

- isoinfo:

  Custom element isotope information that either replaces or adds to the
  definitions in `ecipex()` (see replace_isoinfo). May be a string
  specifying the path to the data file or an existing data frame. See
  Input Details for further information. (character or dataframe)

- replace_isoinfo:

  If TRUE then the information in `isoinfo` replaces the data on the
  ratios of the elemental isotopes used by `ecipex()`. This overrides
  the default isotope ratios. Otherwise, the information in `isoinfo` is
  appended to the data used by `ecipex()`. (boolean)

- min_abundance:

  Minimum proportion retained within candidate isotope patterns.
  (numeric)

- max_isotopes:

  Maximum number of isotopes retained in the pattern of any candidate.
  (numeric)

- skip_isotopes:

  Number of isotopes at the start of the isotope pattern that are
  ignored. Defaults to 0 implying that all isotopes with relative
  abundance greater than `min_abundance` up to `max_isotopes` will be
  retained in the patter of any candidate. (numeric)

- group_formula:

  If TRUE then group candidates with the same chemical formula.
  (boolean)

- group_pattern:

  If TRUE then candidates with similar isotope pattens are grouped when
  computing posterior summary statistics for presence and abundance.
  (boolean)

- pattern_tol:

  Tolerance for grouping candidates with similar isotope patterns. If
  `group_pattern` is true then candidates with a maximum difference of
  their isotope patterns less than `pattern_tol` are grouped when
  computing posterior summary statistics for the presence and abundance.
  (numeric)

- binning:

  Specifies how the peaks in the spectrum and the entries in the
  candidate isotope patterns are associated to construct the design
  matrix. If TRUE then observations in the spectrum and elements of the
  isotope pattern are binned. Otherwise, elements of each isotope
  pattern are associated with the nearest observation in the spectrum up
  to a distance of `epsilon`. (boolean)

- epsilon:

  Tolerance of binning algorithm. (numeric)

- min_mass_charge:

  Lower bound of the mass charge window for the analysis. Any
  observations in the spectrum or isotopes of candidates with
  mass/charge ratio below this value are ignored. (numeric)

- max_mass_charge:

  Upper bound of the mass charge window for the analysis. Any
  observations in the spectrum or isotopes of candidates with
  mass/charge ratio above this value are ignored. (numeric)

- prescreen:

  Minimum count for prescreening candidate. See Prescreening for further
  information. (numeric)

- prescreen_prior:

  Prescreening threshold for prior inclusion probability. Any candidate
  with prior inclusion probability greater than or equal to this value
  will be retained in the analysis regardless of the counts in the
  associated intervals. (numeric)

- prescreen_weight:

  If TRUE then counts are weighted be the relative isotope abundance
  during the prescreening phase. See Prescreening for further
  information. (boolean)

- rounding:

  Method for rounding non-integer counts. Options are "nearest",
  "floor", "ceiling", or "bernoulli". See Rounding for further
  information. (character)

- n.chains:

  Number of parallel chains in MCMC sampler. (numeric)

- n.adapt:

  Number of iterations in MCMC sampler's adapting phase. (numeric)

- n.burnin:

  Number of iterations in MCMC sampler's burn-in phase. (numeric)

- n.sampling:

  Number of iterations in the MCMC sampler's sampling phase. (numeric)

- n.thin:

  Thinning parameter for MCMC sampler. (integer)

- prior_par:

  List of lists specifying the parameters of the prior distributions.
  List.

- model:

  Specifies the form of the model for the abundances of candidates
  present in the sample. Options are "simple" or "hierarchical". See
  Models for further information. (character)

- run_model:

  If TRUE then run MCMC sampler. If FALSE then prepare input without
  runing sampler. Boolean.

- mixtures:

  If TRUE then summarize information about the sampled mixtures.
  (boolean)

- reload_results:

  If TRUE then MCMC sampler is not run and output is retrieved from the
  file specified by `xlsx_out`. Otherwise any existing file at this
  location will be overwritten. (boolean)

- xlsx_out:

  Path including name of output file. (character)

- verbose:

  If TRUE then print progress messages. (boolean)

## Value

A list of objects produced by processing the data and fitting the model.
Objects marked (ro) are only included if `run` is `TRUE`.

- `data`

- `convergence` (ro)

- `burnin` (ro)

- `samples` (ro)

- `coefficients` (ro)

- `intercept` (ro)

- `fitted` (ro)

- `parameters` (ro)

## Input Data

The inputs `spectrum`, `candidates`, `isoinfo`, and `adducts` may either
be specified as strings or data frames. If a string is provided, then it
is interpreted as the path to a file containing the relevant information
and is loaded with
[`read_any()`](https://sjbonner.github.io/SSLAMR/reference/read_any.md).
Otherwise, it is assumed to be data frame that already contains the
requisite information in the correct format.

The following list specifies the format for these arguments (whether
contained in an external file or passed as an existing data frame):

- `spectrum` two columns: `Mass/Charge` and `Count`.

- `candidates` three columns: `Name`, `Formula`, and `Charge`.

- `isoinfo ` four columns `element`, `mass`, `abundance`, and
  `nucleons`. Abundances for the same element must sum to 1.

- `adducts` may either be a string or a data frame. If a string then it
  is interpreted as the path to a file containing the observed spectrum
  and is read with
  [`read_any()`](https://sjbonner.github.io/SSLAMR/reference/read_any.md).
  Otherwise, it is assumed to be a data frame. The data must contain
  three columns: `Name`, `Formula`, and `Action`. Action may either be
  `-` indicating that the fragment may be removed from the candidates or
  `+` indicating that the adduct may be added to the candidates.

## Isotope Data

The `isotope_data` may either be `NULL` (default), a string, or a data
frame. If `NULL` then the relative abundances for the isotope of the
candidates are computed with the `ecipex()` function. If a string then
it is interpreted as the path to a file containing the processed isotope
data and is read with
[`read_any()`](https://sjbonner.github.io/SSLAMR/reference/read_any.md).
Otherwise, it is assumed to be a data frame. The data must contain five
columns: `ID`, `Charge`, `Isotope`, `Mass`, and `Abund`.

## Models

The `model` argument determines the structure of the prior (the model)
for the abundances of candidates present in the sample. If `simple` then
abundances are assigned independent prior distributions with fixed
parameters. If `hierarchical` then the abundances are modelled as draws
from a distribution whose parameters are determined by the data and
assigned further prior distributions.

## Prescreening

Prescreening allows some candidates to be removed from the analysis
before fitting the model. Let \\y_i\\ denote the count for the i-th
entry in the observed spectrum and \\x\_{ij}\\ the relative contribution
of the j-th candidate to this value. If `prescreen_weight` is `FALSE`
then the candidate is include in modelling only if the sum of the counts
in the spectrum associated with the candidate, \$\$\sum\_{i=1}^n y_i
1(x\_{ij} \>0),\$\$ is greater than `prescreen`. If `prescreen_weight`
is `TRUE` then counts are weighted by the relative abundance of the
isotopes, scaled so that the most common isotope has a weight of 1. The
candidate is include in modelling only if the weighted sum,
\$\$\sum\_{i=1}^n y_i \frac{x\_{ij}}{\max\_{j}x\_{ij}}\$\$ is greater
than `prescreen`. Setting `prescreen=0` will include all candidates in
the model.

## Rounding

Intensities in the observed spectrum (`Count`) are assumed to be integer
counts. If non-integer values are detected then these may be rounded in
different ways. Options are:

- `nearest` round to the nearest integer in the usual way.

- `floor` round down to the next highest integer.

- `ceiling` round up to the next lowest integer.

- `bernoulli` if \\d\\ is the decimal part of the count then round up
  with probability \\d\\ and down with probability \\1-p\\.

The advantage of `bernoulli` is that the expected value of the sum of
the rounded values is equal to the sum of the original values.

## Examples

``` r
#TBD
```
