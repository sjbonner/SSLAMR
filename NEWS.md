# SSLAMR 0.1.7

* Added convergence diagnostics to Excel output.
* Revised convergence diagnostics to include only the final half of the burn-in. 
* Added argument to thing Markov chains.
* Changed default prior probability of inclusion from .5 to .1.
* Further minor bug fixes.

# SSLAMR 0.1.6

* Updated prescreen to account for prior inclusion probabilities. Candidates are now retained if their prior inclusion probability is greater than prescreen_prior regardless of the peaks in the associated intervals. This ensures that candidates with high prior probabilities are retained. The value of prescreen_prior is 1 by default. 

# SSLAMR 0.1.5

* Fixed bug in generating initial values when inclusion parameters have separate priors. Initial scenario now correctly keep any candidates with a prior inclusion of 1.

# SSLAMR 0.1.4

* Added functionality to allow separate priors on gamma (inclusion) parameters.
* Constrained model so that adducts can only appear if the parent is also included.

# SSLAMR 0.1.3

* Added function to summarize mixtures and abundance of components conditional on mixture. 

# SSLAMR 0.1.2

* Updated `greedy` algorithm for comparison. Primary change is to allow prescreening in previous models that feed into the algorithm.
* Added SSLAMR version to output file.

# SSLAMR 0.1.1

* Added additional timing steps and saved timing information in output.
* Retained chemical formula in candidate information.
* Scaled weights when `prescreen_weight == TRUE` so that weights are proportional to abundance and weight of the most abundant isotope is 1.
