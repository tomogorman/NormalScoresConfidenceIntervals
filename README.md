
This repository contains the computer output of simulations that estimate 
the coverage probabilities and the average width of confidence intervals 
that are based on the inversion of the generalized normal scores (GNS) test. 

The confidence interval method is described in

O'Gorman, T. W. (2020) A Method to Reduce the Width of Confidence Intervals  
by Using a Normal Scores Transformation. (Submitted for publication.)

The GNS test is described in

O'Gorman, T. W. (2020) A Generalized Normal Scores Test of Significance
for any Coefficient in a Linear Model. Communications In Statistics--Theory
and Methods (https: (https://doi.org/10.1080/03610926.2021.1987471)

In each of these simulations 100000 data sets were generated. For each of these 
data sets the traditional confidence interval and the Generalized Normal 
Scores confidence were computed. Based on these 100000 data sets the 
empirical coverage probability and the average width were computed.

The files are:

|    File name        |               Model                                     |
|---------------------|:-------------------------------------------------------:|
|nsci.ic2.bal.txt     | For two equal samples from two populations              |
|nsci.ic2.nor.txt     | Simple Linear Regression with one Normal covariate      |
|nsci.ic2.ln.txt      | Simple Linear Regression with one Lognormal covariate   |
|nsci.ic3.nor.r0.txt  | Regression with 2 normal covariates correlated r = .0   |
|nsci.ic3.ln.r0.txt   | Regression with 2 lognormal covariates correlated r = .0|
|nsci.ic3.nor.r4.txt  | Regression with 2 normal covariates correlated r = .4   | 
|nsci.ic3.ln.r4.txt   | Regression with 2 lognormal covariates correlated r = .4| 
|nsci.ic3.nor.r8.txt  | Regression with 2 normal covariates correlated r = .8   | 
|nsci.ic3.ln.r8.txt   | Regression with 2 lognormal covariates correlated r = .8| 
|nsci.ic5.nor.r4.txt  | Regression with 4 normal covariates correlated r = .4   | 
|nsci.ic5.ln.r4.txt   | Regression with 4 lognormal covariates correlated r = .4|

The second group of file give the coverage probabilities and the average
width of the confidence intervals when the covariates are fixed.
                                                                              
|    File name        |               Model                                     |
|---------------------|:-------------------------------------------------------:|
|nsci.ic2.x21.txt     | Simple Regression, 1 fixed covariate, equal spacing     |
|nsci.ic2.x22.txt     | Simple Regression, 1 fixed covariate, unequal spacing   |
|nsci.ic3.x23.txt     | Mult. Reg., 2 uncorrelated covariates, equal spacing    |
|nsci.ic3.x24.txt     | Mult. Reg., 2 correlated covariates, equal spacing      |
|nsci.ic3.x25.txt     | Mult. Reg., 2 uncorrelated covariates, unequal spacing  |
|nsci.ic3.x26.txt     | Mult. Reg., 2 correlated covariates, unequal spacing    |





The following notation is used in these files:

"n" is the number of observations.

"Distribution" is the distribution of the errors that were generated.

"OLS" the traditional confidence interval based on the t distribution. 

"CP" is the empircal coverage probability.

"NS" is the Generalized Normal Scores (GNS) confidence interval. 

"PctDec" is the percent decrease in width by using the GNS confidence interval.


