
<!-- README.md is generated from README.Rmd. Please edit that file -->
pdmppoly
========

The goal of package `pdmppoly` is to simulate polynomial piecewise deterministic markov processes (polynomial PDMPs) within `R` and to provide methods to calculate or approximate their moments.

In additon, it contains all models that are used in the (not yet published) doctoral thesis of Charlotte Jana.

Installation
------------

You can install pdmppoly from github with:

``` r
# install.packages("devtools")
devtools::install_github("CharlotteJana/pdmppoly")
```

Polynomial PDMPs
----------------

A *polynomial PDMP* is a piecewise deterministic Markov process (PDMP), whose dynamics and rates are all polynomials, depending on the different variables and parameters of the PDMP.

A polynomial PDMP in `pdmppoly` is represented by an object of class `polyPdmpModel`. This class is very similar to class `pdmpModel` of package [pdmpsim](https://github.com/CharlotteJana/pdmpsim) and every method provided by [pdmpsim](https://github.com/CharlotteJana/pdmpsim) can be applied to `polyPdmpModel`-objects. The polynomials appear in the new slots `dynfunc` and `ratefunc`. They are defined using package [spray](https://github.com/RobinHankin/spray.git) and internally represented as sparse coefficient matrices.

This is a simple example modelling gene expression with positive feedback:

``` r
genePolyBF <- new("polyPdmpModel",
                  descr = "Gene regulation with positive feedback",
                  parms = list(b = 0.5, a0 = 1, a1 = 3, k10 = 1, k01 = 0.5), 
                  init = c(f = 1, d = 1), 
                  discStates = list(d = 0:1),
                  dynpolys = quote(list(
                    list(overall = -b*lone(1,2),
                         specific = list(a0, a1))
                  )),
                  ratepolys = quote(list(  
                    list(k01*lone(1,2), k10)
                  )),
                  jumpfunc = function(t, x, parms, jtype) {
                    c(x[1], 1 - x[2])
                  }, 
                  times = c(from = 0, to = 25, by = 0.1), 
                  solver = "lsodar")
```

Moment calculation
------------------

Function `momApp` calculates the moments of the distribution of of a polynomial PDMP, up to a given order. It solves a system of differential equaitons with package `deSove` to obtain the results. This system is usually not closed and has to be altered by replacing moments of higher orders with fixed values. The moments to replace can be raw or central moments, depending on the argument `centralize`. The following closure methods are available:

-   `closure = "zero"`: Replace the moments with 0
-   `closure = "normal"`: Replace with moments of a normal distribution
-   `closure = "lognormal"`: Replace with moments of a lognormal distribution
-   `closure = "gamma"`: Replace with moments of a gamma distribution

``` r
mom <- momApp(genePolyBF, maxOrder = 8,
              closure = c("zero", "normal", "zero"), 
              centralize = c(FALSE, TRUE, TRUE))
mom
#> Model: 
#> Description: Gene regulation with positive feedback
#> Parameters: b = 0.5, a0 = 1, a1 = 3, k10 = 1, k01 = 0.5
#> Initial Values: f = 1, d = 1
#> 
#> Moment approximation for moments of order > 8 leads to 
#> 
#>              method order time         d            f
#> 1        zero (raw)     1   25 0.6209922 4.484764e+00
#> 9  normal (central)     1   25 0.6792388 4.716683e+00
#> 17   zero (central)     1   25 0.6792388 4.716683e+00
#> 2        zero (raw)     2   25 0.6209922 2.194537e+01
#> 10 normal (central)     2   25 0.6792388 2.286486e+01
#> 18   zero (central)     2   25 0.6792388 2.286486e+01
#> 3        zero (raw)     3   25 0.6209922 1.076902e+02
#> 11 normal (central)     3   25 0.6792388 1.131856e+02
#> 19   zero (central)     3   25 0.6792388 1.131856e+02
#> 4        zero (raw)     4   25 0.6209922 5.468273e+02
#> 12 normal (central)     4   25 0.6792388 5.704221e+02
#> 20   zero (central)     4   25 0.6792388 5.704221e+02
#> 5        zero (raw)     5   25 0.6209922 2.773017e+03
#> 13 normal (central)     5   25 0.6792388 2.916249e+03
#> 21   zero (central)     5   25 0.6792388 2.916249e+03
#> 6        zero (raw)     6   25 0.6209922 1.452112e+04
#> 14 normal (central)     6   25 0.6792388 1.510106e+04
#> 22   zero (central)     6   25 0.6792388 1.510106e+04
#> 7        zero (raw)     7   25 0.6209922 7.470581e+04
#> 15 normal (central)     7   25 0.6792388 7.898200e+04
#> 23   zero (central)     7   25 0.6792388 7.898200e+04
#> 8        zero (raw)     8   25 0.6209922 4.061504e+05
#> 16 normal (central)     8   25 0.6792388 4.172857e+05
#> 24   zero (central)     8   25 0.6792388 4.172857e+05
```

``` r
plot(mom, plotorder = 1:2)
```

![](man/figures/README-unnamed-chunk-5-1.png)

### Compare with empirical moments

The results of `momApp` are approximated values of the moments of the distribution. To test their accuracy, it is best to simulate the model multiple times, calculate the empirical moments and compare them with the approximated moments.

``` r
simulations <- multSim(genePolyF, seeds = 1:30)
mom <- addSimulation(mom, simulations)
plot(mom, plotorder = 1:4, vars = "f")
```

Test if the distribution is unimodal
------------------------------------

The results of function `momApp` can also be used to test if the distribution of the PDMP is unimodal. Function `is.unimodal` from package `momcalc` performs this test for one dimensional distributions. It is only meaningful for distributions with a known compact support. Package `pdmppoly` provides function `modalityTest` to apply the test directly on a data.frame generated by `momApp`. It performs the test for every time value and every variable of the PDMP seperately.

``` r
testresults <- modalityTest(mom, lower = 0, upper = 35, vars = "f")
head(testresults)
#>   time     method variable     modality
#> 1  0.0 zero (raw)        f         <NA>
#> 2  0.1 zero (raw)        f not unimodal
#> 3  0.2 zero (raw)        f not unimodal
#> 4  0.3 zero (raw)        f not unimodal
#> 5  0.4 zero (raw)        f not unimodal
#> 6  0.5 zero (raw)        f not unimodal
```

``` r
plotModalities(mom, modalities = testresults)
```

![](man/figures/README-unnamed-chunk-8-1.png)

Be aware that both `momApp` and `modalityTest` do not depend on time-consuming simulations!

Analyse a polynomial PDMP
-------------------------

Package `pdmppoly` provides a function `analyseModel` that analyses a polynomial PDMP as good as possible, meaning that it performs most of the functions available on packages `pdmpsim` and `pdmppoly`. It

-   simulates the model multiple times,
-   plots the simulated data with all available plot methods,
-   calculates statistics such as mean, sd, median, etc. ,
-   performs moment approximation with all available closure methods,
-   tests if the distribution is unimodal and plots the result.

License
-------

GPL 2 or higher
