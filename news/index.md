# Changelog

## SSLAMR 0.1.15 (January 26, 2026)

- Fixed bug with labelling of groups within mixtures. Group names were
  based on the similarity of the candidates isotope patterns regardless
  of whether or not a specific candidate was present on a given
  iteration of the MCMC sampler. Groups are still identified in the same
  way, but the group names on each iteration are formed by concatenating
  the names of only those candidates that are present in the mixture.

## SSLAMR 0.1.14 (January 16, 2026)

- Happy New Year!
- Re-implement grouping based on isotope pattern for the greedy
  algorithm. Some candidates were being included in groups even if their
  parents weren’t included. This should not be allowed.

## SSLAMR 0.1.12 (May 8, 2025)

- Incorporated grouping of the candidates into the greedy algorithm.
  This avoids overpenalizing the greedy algorithm when it fits one
  candidate out of a group and discards the rest.

## SSLAMR 0.1.11 (May 6, 2025)

- Added grouping of candidates with similar isotope patterns when
  computing posterior summary statistics for mixtures.

## SSLAMR 0.1.10

- Restricted number of mixtures saved to Excel output file.

## SSLAMR 0.1.9

- Adding grouping of candidates with similar isotope patterns when
  computing posterior summary statistics for presence and abundance.

## SSLAMR 0.1.8

- Added documentation for and exported simulate() function.

## SSLAMR 0.1.7

- Added convergence diagnostics to Excel output.
- Revised convergence diagnostics to include only the final half of the
  burn-in.
- Added argument to thing Markov chains.
- Changed default prior probability of inclusion from .5 to .1.
- Further minor bug fixes.

## SSLAMR 0.1.6

- Updated prescreen to account for prior inclusion probabilities.
  Candidates are now retained if their prior inclusion probability is
  greater than prescreen_prior regardless of the peaks in the associated
  intervals. This ensures that candidates with high prior probabilities
  are retained. The value of prescreen_prior is 1 by default.

## SSLAMR 0.1.5

- Fixed bug in generating initial values when inclusion parameters have
  separate priors. Initial scenario now correctly keep any candidates
  with a prior inclusion of 1.

## SSLAMR 0.1.4

- Added functionality to allow separate priors on gamma (inclusion)
  parameters.
- Constrained model so that adducts can only appear if the parent is
  also included.

## SSLAMR 0.1.3

- Added function to summarize mixtures and abundance of components
  conditional on mixture.

## SSLAMR 0.1.2

- Updated `greedy` algorithm for comparison. Primary change is to allow
  prescreening in previous models that feed into the algorithm.
- Added SSLAMR version to output file.

## SSLAMR 0.1.1

- Added additional timing steps and saved timing information in output.
- Retained chemical formula in candidate information.
- Scaled weights when `prescreen_weight == TRUE` so that weights are
  proportional to abundance and weight of the most abundant isotope is
  1.
