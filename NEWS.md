# gasmodel (development version)
* Added the Burr distribution.
* Added the Fisk distribution.
* Added the Lomax distribution.
* Increased the running speed of the estimation.
* Corrected the log-likelihood function for the zero-inflated geometric distribution.

# gasmodel 0.5.1
* Added block bootstrap methods.

# gasmodel 0.5.0
* Added support for summarize() and plot() generic functions.
* Added computation of p-values in the gas_bootstrap() function.
* Added functionality for parallelization in the gas_bootstrap() function.

# gasmodel 0.4.0
* Added the asymmetric Laplace distribution.
* Added the Birnbaum–Saunders distribution.
* Changed the default value of par_static. Now, the first group of parameters (usually just the first parameter for univariate distributions) is dynamic while the rest are static.
* Renamed the "spec" argument to "regress". One of its value, "reg_err", renamed to "sep".
* Corrected the Fisher information for the generalized gamma and zero-inflated negative binomial distributions.
* Updated the ice_hockey_championships dataset.

# gasmodel 0.3.0
* Added diagonal and full Fisher information matrix scalings.
* Added the Beta distribution.
* Added the Dirichlet distribution.
* Added the Laplace distribution.
* Added the von Mises distribution.
* Added german_car_market_cap dataset.
* Corrected NA and Inf handling for the Plackett-Luce distribution.
* Corrected the Fisher information for the generalized gamma distribution.
* Replaced the tidyverse package with dplyr, tidyr, and ggplot2 packages in Suggests.

# gasmodel 0.2.0
* Added expected values, variances, and residuals to the return value of gas().
* Added the zero-inflated Skellam distribution.
* Modified the computation of the expected value, the variance, and the Fisher information of the Plackett-Luce distribution.
* Modified the computation of the Fisher information of the zero-inflated negative binomial distribution.
* Updated the sp500_daily dataset.

# gasmodel 0.1.0
* Initial release on CRAN.
