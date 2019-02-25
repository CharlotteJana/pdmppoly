
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pdmppoly

The goal of pdmppoly is to simulate polynomial piecewise deterministic
markov processes (PDMPs) within R and to provide methods to calculate or
approximate their moments.

In additon, it contains all methods that are not part of package pdmpsim
but are used in the (not yet published) doctoral thesis of Charlotte
Jana.

## Installation

You can install pdmppoly from github with:

``` r
# install.packages("devtools")
devtools::install_github("CharlotteJana/pdmppoly")
```

## Example

This is a simple example modelling gene expression with positive
feedback:

``` r
genePolyF <- new("polyPdmpModel",
     descr = "Gene regulation with positive feedback (polynomial version)",
     parms = list(b = 0.2, a = 7, k10 = 0.04, k01 = 0.02), 
     init = c(f = 1, d = 1), 
     discStates = list(d = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-b,a)))
     )),
     ratepolys = quote(list(  
       list(k01*lone(1,2), k10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     }, 
     times = c(from = 0, to = 100, by = 0.1), 
     solver = "lsodar")
```

Calculate the first 4 moments:

``` r
mcalc <- momApp(genePolyF, maxOrder = 4)$moments
```

Simulate multiple times and calculate the empirical moments:

``` r
simulations <- multSim(genePolyF, seeds = 1:30)
msim.list <- lapply(1:4, function(i) moments(simulations, order = i))
msim <- do.call(rbind, msim.list)
```
