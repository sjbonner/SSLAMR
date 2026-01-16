# Mass Spectrum Simulation

Mass Spectrum Simulation

## Usage

``` r
sslamr_simulate(
  candidates,
  isoinfo = NULL,
  mc_cal_error,
  mc_sd,
  mc_min_diff,
  intercept,
  mc_min,
  mc_max
)
```

## Arguments

- candidates:

  List of candidate compounds included in the mixture.

- isoinfo:

  Custom element isotope information that either replaces or adds to the
  definitions in `ecipex()`.

- mc_cal_error:

  Mass/charge calibration error. This error is assumed to be systematic
  and is added directly to the mass/charge values for all peaks.

- mc_sd:

  Standard deviation of the mass/charge observation error.

- mc_min_diff:

  Minimum discernable mass/charge difference. Peaks separated by less
  than this value are combined.

- intercept:

  Mean count of the noise per amu.

- mc_min:

  Lower bound of the mass/charge window.

- mc_max:

  Upper bound of the mass/charge window.

## Value

A tibble with the simulated mass spectrum containing columns for the
mass/charge ratio and intensity.

## Examples

``` r
# TBD
```
