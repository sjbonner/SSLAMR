# SSLAMR 0.1.3

* Added function to summarize mixtures and abundance of components conditional on mixture. 

# SSLAMR 0.1.2

* Updated `greedy` algorithm for comparison. Primary change is to allow prescreening in previous models that feed into the algorithm.
* Added SSLAMR version to output file.

# SSLAMR 0.1.1

* Added additional timing steps and saved timing information in output.
* Retained chemical formula in candidate information.
* Scaled weights when `prescreen_weight == TRUE` so that weights are proportional to abundance and weight of the most abundant isotope is 1.
