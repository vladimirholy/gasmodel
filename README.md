
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gasmodel

<!-- badges: start -->
<!-- badges: end -->

# Overview

A package for estimation, forecasting, and simulation of generalized
autoregressive score (GAS) models of Creal et al. (2013) and Harvey
(2013), also known as dynamic conditional score (DCS) models or
score-driven (SD) models.

Model specification allows for various conditional distributions,
different parametrizations, exogenous variables, higher score and
autoregressive orders, custom and unconditional initial values of
time-varying parameters, fixed and bounded values of coefficients, and
NA values. Model estimation is performed by the maximum likelihood
method and the Hessian matrix.

The package offers the following functions for working with GAS models:

-   `gas()` estimates GAS models.
-   `gas_simulate()` simulates GAS models.
-   `gas_forecast()` forecasts GAS models.
-   `gas_filter()` obtains filtered time-varying parameters of GAS
    models.
-   `gas_bootstrap()` bootstraps coefficients of GAS models.

Probability distributions are handled by the following functions:

-   `distr()` provides table of supported distributions.
-   `distr_density()` computes the density.
-   `distr_mean()` computes the mean.
-   `distr_var()` computes the variance.
-   `distr_score()` computes the score.
-   `distr_fisher()` computes the Fisher information.
-   `distr_random()` generates random observations.

In addition, the package contains the following datasets used in
examples:

-   `bookshop_sales` ???
-   `ice_hockey_championships` ???

## Installation

You can install the development version of gasmodel from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vladimirholy/gasmodel")
```

## Supported Distributions

Currently, there are 17 distributions available.

``` r
library(gasmodel)
print(distr(), right = FALSE, row.names = FALSE)
#>  distr_title                     param_title   distr     param type        dim   orthog default
#>  Bernoulli                       Standard      bernoulli std   binary      uni    TRUE   TRUE  
#>  Categorical                     Worth         cat       worth categorical multi FALSE   TRUE  
#>  Double Poisson                  Standard      dpois     std   count       uni    TRUE   TRUE  
#>  Exponential                     Rate          exp       rate  duration    uni    TRUE  FALSE  
#>  Exponential                     Scale         exp       scale duration    uni    TRUE   TRUE  
#>  Gamma                           Rate          gamma     rate  duration    uni   FALSE  FALSE  
#>  Gamma                           Scale         gamma     scale duration    uni   FALSE   TRUE  
#>  Generalized Gamma               Rate          gengamma  rate  duration    uni   FALSE  FALSE  
#>  Generalized Gamma               Scale         gengamma  scale duration    uni   FALSE   TRUE  
#>  Geometric                       Mean          geom      mean  count       uni    TRUE   TRUE  
#>  Geometric                       Probabilistic geom      prob  count       uni    TRUE  FALSE  
#>  Multivariate Normal             Standard      mnorm     std   real        multi FALSE   TRUE  
#>  Negative Binomial               NB2           negbin    nb2   count       uni    TRUE   TRUE  
#>  Negative Binomial               Probabilistic negbin    prob  count       uni   FALSE  FALSE  
#>  Normal                          Standard      norm      std   real        uni    TRUE   TRUE  
#>  Plackett-Luce                   Standard      pluce     std   ranking     multi FALSE   TRUE  
#>  Poisson                         Standard      pois      std   count       uni    TRUE   TRUE  
#>  Student‘s t                     Standard      t         std   real        uni   FALSE   TRUE  
#>  Weibull                         Rate          weibull   rate  duration    uni   FALSE  FALSE  
#>  Weibull                         Scale         weibull   scale duration    uni   FALSE   TRUE  
#>  Zero-Inflated Geometric         Mean          zigeom    mean  count       uni   FALSE   TRUE  
#>  Zero-Inflated Negative Binomial NB2           zinegbin  nb2   count       uni   FALSE   TRUE  
#>  Zero-Inflated Poisson           Standard      zipois    std   count       uni   FALSE   TRUE
```

## Generalized Autoregressive Score Models

The generalized autoregressive score (GAS) models of Creal et al. (2013)
and Harvey (2013), also known as dynamic conditional score (DCS) models
or score-driven (SD) models, have established themselves as a useful
modern framework for time series modeling.

The GAS models are observation-driven models allowing for any underlying
probability distribution
![p(y_t\|f_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%28y_t%7Cf_t%29 "p(y_t|f_t)")
with any time-varying parameters
![f_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_t "f_t")
for time series
![y_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_t "y_t").
They capture the dynamics of time-varying parameters using the
autoregressive term and the lagged score, i.e. the gradient of the
log-likelihood function. Exogenous variables can also be included.
Specifically, time-varying parameters
![f\_{t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_%7Bt%7D "f_{t}")
follow the recursion

![f\_{t} = \\omega + \\sum\_{i=1}^M \\beta_i x\_{ti} + \\sum\_{j=1}^P \\alpha_j S(f\_{t - j}) \\nabla(y\_{t - j}, f\_{t - j}) + \\sum\_{k=1}^Q \\varphi_k f\_{t-k},](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_%7Bt%7D%20%3D%20%5Comega%20%2B%20%5Csum_%7Bi%3D1%7D%5EM%20%5Cbeta_i%20x_%7Bti%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5EP%20%5Calpha_j%20S%28f_%7Bt%20-%20j%7D%29%20%5Cnabla%28y_%7Bt%20-%20j%7D%2C%20f_%7Bt%20-%20j%7D%29%20%2B%20%5Csum_%7Bk%3D1%7D%5EQ%20%5Cvarphi_k%20f_%7Bt-k%7D%2C "f_{t} = \omega + \sum_{i=1}^M \beta_i x_{ti} + \sum_{j=1}^P \alpha_j S(f_{t - j}) \nabla(y_{t - j}, f_{t - j}) + \sum_{k=1}^Q \varphi_k f_{t-k},")

where
![\\omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Comega "\omega")
is a vector of constants,
![\\beta_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_i "\beta_i")
are regression parameters,
![\\alpha_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_j "\alpha_j")
are score parameters,
![\\varphi_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvarphi_k "\varphi_k")
are autoregressive parameters,
![x\_{ti}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_%7Bti%7D "x_{ti}")
are exogenous variables,
![S(f_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%28f_t%29 "S(f_t)")
is a scaling function for the score, and
![\\nabla(y_t, f_t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnabla%28y_t%2C%20f_t%29 "\nabla(y_t, f_t)")
is the score given by

![\\nabla(y_t, f_t) = \\frac{\\partial \\ln p(y_t \| f_t)}{\\partial f_t}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnabla%28y_t%2C%20f_t%29%20%3D%20%5Cfrac%7B%5Cpartial%20%5Cln%20p%28y_t%20%7C%20f_t%29%7D%7B%5Cpartial%20f_t%7D. "\nabla(y_t, f_t) = \frac{\partial \ln p(y_t | f_t)}{\partial f_t}.")

The GAS models can be straightforwardly estimated by the maximum
likelihood method. For the asymptotic theory regarding the GAS models
and maximum likelihood estimation, see Blasques et al. (2014), Blasques
et al. (2018), and Blasques et al. (2022).

The use of the score for updating time-varying parameters is optimal in
an information theoretic sense. For an investigation of the optimality
properties of GAS models, see Blasques et al. (2015) and Blasques et
al. (2021).

Generally, the GAS models perform quite well when compared to
alternatives, including parameter-driven models. For a comparison of the
GAS models to alternative models, see Koopman et al. (2016) and Blazsek
and Licht (2020).

The GAS class includes many well-known econometric models, such as the
generalized autoregressive conditional heteroskedasticity (GARCH) model
of Bollerslev (1986), the autoregressive conditional duration (ACD)
model of Engle and Russel (1998), and the Poisson count model of Davis
et al. (2003). More recently, a variety of novel score-driven models has
been proposed, such as the Beta-t-(E)GARCH model of Harvey and
Chakravarty (2008), the discrete price changes model of Koopman et
al. (2018), the directional model of Harvey (2019), the bivariate
Poisson model of Koopman and Lit (2019), and the ranking model of Holý
and Zouhar (2021). For an overview of various GAS models, see Harvey
(2022).

The extensive GAS literature is listed on
[www.gasmodel.com](http://www.gasmodel.com).

## References

Blasques, F., Gorgi, P., Koopman, S. J., and Wintenberger, O. (2018).
Feasible Invertibility Conditions and Maximum Likelihood Estimation for
Observation-Driven Models. *Electronic Journal of Statistics*,
**12**(1), 1019–1052. doi:
[10.1214/18-ejs1416](https://doi.org/10.1214/18-ejs1416).

Blasques, F., Koopman, S. J., and Lucas, A. (2014). Stationarity and
Ergodicity of Univariate Generalized Autoregressive Score Processes.
*Electronic Journal of Statistics*, **8**(1), 1088–1112. doi:
[10.1214/14-ejs924](https://doi.org/10.1214/14-ejs924).

Blasques, F., Koopman, S. J., and Lucas, A. (2015).
Information-Theoretic Optimality of Observation-Driven Time Series
Models for Continuous Responses. *Biometrika*, **102**(2), 325–343. doi:
[10.1093/biomet/asu076](https://doi.org/10.1093/biomet/asu076).

Blasques, F., Lucas, A., and van Vlodrop, A. C. (2021). Finite Sample
Optimality of Score-Driven Volatility Models: Some Monte Carlo Evidence.
*Econometrics and Statistics*, **19**, 47–57. doi:
[10.1016/j.ecosta.2020.03.010](https://doi.org/10.1016/j.ecosta.2020.03.010).

Blasques, F., van Brummelen, J., Koopman, S. J., and Lucas, A. (2022).
Maximum Likelihood Estimation for Score-Driven Models. *Journal of
Econometrics*, **227**(2), 325–346. doi:
[10.1016/j.jeconom.2021.06.003](https://doi.org/10.1016/j.jeconom.2021.06.003).

Blazsek, S. and Licht, A. (2020). Dynamic Conditional Score Models: A
Review of Their Applications. *Applied Economics*, **52**(11),
1181–1199. doi:
\[10.1080/00036846.2019.1659498}(<https://doi.org/10.1080/00036846.2019.1659498>).

Bollerslev, T. (1986). Generalized Autoregressive Conditional
Heteroskedasticity. *Journal of Econometrics*, **31**(3), 307–327. doi:
[10.1016/0304-4076(86)90063-1](https://doi.org/10.1016/0304-4076(86)90063-1).

Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized
Autoregressive Score Models with Applications. *Journal of Applied
Econometrics*, **28**(5), 777–795. doi:
[10.1002/jae.1279](https://doi.org/10.1002/jae.1279).

Davis, R. A., Dunsmuir, W. T. M., and Street, S. B. (2003).
Observation-Driven Models for Poisson Counts. *Biometrika*, **90**(4),
777–790. doi:
[10.1093/biomet/90.4.777](https://doi.org/10.1093/biomet/90.4.777).

Engle, R. F. and Russell, J. R. (1998). Autoregressive Conditional
Duration: A New Model for Irregularly Spaced Transaction Data.
*Econometrica*, **66**(5), 1127–1162. doi:
[10.2307/2999632](https://doi.org/10.2307/2999632).

Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy Tails:
With Applications to Financial and Economic Time Series*. Cambridge
University Press. doi:
[10.1017/cbo9781139540933](https://doi.org/10.1017/cbo9781139540933).

Harvey, A. C. (2022). Score-Driven Time Series Models. *Annual Review of
Statistics and Its Application*, **9**(1), 321–342. doi:
[10.1146/annurev-statistics-040120-021023](https://doi.org/10.1146/annurev-statistics-040120-021023).

Harvey, A. C. and Chakravarty, T. (2008). Beta-t-(E)GARCH. *Cambridge
Working Papers in Economics*, CWPE 0840. doi:
[10.17863/cam.5286](https://doi.org/10.17863/cam.5286).

Harvey, A., Hurn, S., and Thiele, S. (2019). Modeling Directional
(Circular) Time Series. *Cambridge Working Papers in Economics*, CWPE
1971. doi: [10.17863/cam.43915](https://doi.org/10.17863/cam.43915).

Holý, V. and Zouhar, J. (2021). Modelling Time-Varying Rankings with
Autoregressive and Score-Driven Dynamics. arXiv:
[2101.04040](https://arxiv.org/abs/2101.04040).

Koopman, S. J. and Lit, R. (2019). Forecasting Football Match Results in
National League Competitions Using Score-Driven Time Series Models.
*International Journal of Forecasting*, **35**(2), 797–809. doi:
[10.1016/j.ijforecast.2018.10.011](https://doi.org/10.1016/j.ijforecast.2018.10.011).

Koopman, S. J., Lit, R., Lucas, A., and Opschoor, A. (2018). Dynamic
Discrete Copula Models for High-Frequency Stock Price Changes. *Journal
of Applied Econometrics*, **33**(7), 966–985. doi:
[10.1002/jae.2645](https://doi.org/10.1002/jae.2645).

Koopman, S. J., Lucas, A., and Scharth, M. (2016). Predicting
Time-Varying Parameters with Parameter-Driven and Observation-Driven
Models. *Review of Economics and Statistics*, **98**(1), 97–110. doi:
[10.1162/rest_a\_00533](https://doi.org/10.1162/rest_a_00533).
