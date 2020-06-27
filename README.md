
This repository contains the computer output of simulations that estimate 
the coverage probabilities and the average width of confidence intervals 
that are based on the inversion of the generalized normal scores (GNS) test. 

The confidence interval method is described in:

O'Gorman, T. W. (2020) A Method to Reduce the Width of Confidence Intervals while Maintaining 
Their Coverage Probability. (Submitted for publication.)

The GNS test is described in:

O'Gorman, T. W. (2020) A Generalized Normal Scores Test of Significance
for any Coefficient in a Linear Model. (Submitted for publication.)

In each of these simulations 100000 data sets were generated. For each of these 
data sets the traditional confidence interval and the Generalized Normal 
Scores confidence were computed. Based on these 100000 data sets the 
empirical coverage probability and the average width were computed.

The files are:

        File name                   Model
______________________________________________________________________________

nsci.ic2.bal.txt      For two equal samples from two populations

nsci.ic2.nor.txt      Simple Linear Regression with one Normal covariate

nsci.ic2.ln.txt       Simple Linear Regression with one Lognormal covariate

nsci.ic3.nor.r0.txt   Regression with 2 normal covariates correlated r = .0

nsci.ic3.ln.r0.txt    Regression with 2 lognormal covariates correlated r = .0

nsci.ic3.nor.r4.txt   Regression with 2 normal covariates correlated r = .4
 
nsci.ic3.ln.r4.txt    Regression with 2 lognormal covariates correlated r = .4
 
nsci.ic3.nor.r8.txt   Regression with 2 normal covariates correlated r = .8

nsci.ic3.ln.r8.txt    Regression with 2 lognormal covariates correlated r = .8

nsci.ic5.nor.r4.txt   Regression with 4 normal covariates correlated r = .4

nsci.ic5.ln.r4.txt    Regression with 4 lognormal covariates correlated r = .4 




The following notation is used in these files:

"n" is the number of observations.

"Distribution" is the distribution of the errors that were generated.

"OLS" the traditional confidence interval based on the t distribution. 

"CP" is the empircal coverage probability.

"NS" is the Generalized Normal Scores (GNS) confidence interval. 

"PctDec" is the percent decrease in width by using the GNS confidence interval.



