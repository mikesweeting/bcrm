<!-- README.md is generated from README.Rmd. Please edit that file -->
bcrm
====

Implements a wide variety of one and two-parameter Bayesian CRM designs. The program can run interactively, allowing the user to enter outcomes after each cohort has been recruited, or via simulation to assess operating characteristics.

Installation
------------

You can install from CRAN with:

``` r
install.packages("bcrm")
```

Or try the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("mikesweeting/bcrm")
```

Example
-------

``` r
library(bcrm)

## Dose-escalation cancer trial example as described in Neuenschwander et al 2008.
## Pre-defined doses
dose<-c(1,2.5,5,10,15,20,25,30,40,50,75,100,150,200,250)
## Pre-specified probabilities of toxicity
## [dose levels 11-15 not specified in the paper, and are for illustration only]
p.tox0<-c(0.010,0.015,0.020,0.025,0.030,0.040,0.050,0.100,0.170,0.300,0.400,0.500,0.650
  ,0.800,0.900)
## Data from the first 5 cohorts of 18 patients
data<-data.frame(patient=1:18,dose=rep(c(1:4,7),c(3,4,5,4,2)),tox=rep(0:1,c(16,2)))
## Target toxicity level
target.tox<-0.30

## Simulate 10 replicate trials of size 36 (cohort size 3) using this design 
## with constraint (i.e. no dose-skipping) and starting at lowest dose
## True probabilities of toxicity are set to pre-specified probabilities (p.tox0) 
Power.LN.bcrm.sim<-bcrm(stop=list(nmax=36),p.tox0=p.tox0,dose=dose,ff="power"
  ,prior.alpha=list(3,0,1.34^2),target.tox=target.tox,constrain=TRUE
  ,sdose.calculate="median",pointest="mean",start=1,simulate=TRUE,nsims=10,truep=p.tox0)
#> 1 
#> 2 
#> 3 
#> 4 
#> 5 
#> 6 
#> 7 
#> 8 
#> 9 
#> 10
print(Power.LN.bcrm.sim)
#> Operating characteristics based on  10  simulations: 
#>  
#>               
#> Sample size 36
#> 
#>                             Doses
#>                              No dose      1    2.5      5     10     15
#>   Experimentation proportion      NA 0.0833 0.0833 0.0833 0.0833 0.0833
#>   Recommendation proportion        0 0.0000 0.0000 0.0000 0.0000 0.0000
#>                             Doses
#>                                  20     25     30    40   50     75
#>   Experimentation proportion 0.0833 0.0833 0.0917 0.125 0.15 0.0417
#>   Recommendation proportion  0.0000 0.0000 0.0000 0.200 0.50 0.2000
#>                             Doses
#>                                  100 150 200 250
#>   Experimentation proportion 0.00833   0   0   0
#>   Recommendation proportion  0.10000   0   0   0
#> 
#>                             Probability of DLT
#>                              [0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8] (0.8,1]
#>   Experimentation proportion     0.8     0.192   0.00833         0       0
#>   Recommendation proportion      0.2     0.700   0.10000         0       0
```
