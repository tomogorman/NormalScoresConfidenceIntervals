#
#                                 Program:  nsci.r
#                                 Revision Date: July 9, 2020
#
#  This nscilimits function computes 95% confidence intervals limits
#  for a coefficient in a linear model having fixed effects.
#  In this function, the residuals from the reduced model
#  are replaced by their normal scores, which are then standardized,
#  and the scores are added to the predicted values from the 
#  reduced model to find a new vector of y-values.
#  These y-values are then used a linear model to find the p-value for the
#  test of significance for the last parameter in the linear model.  
#  The confidence interval limits are found by inverting the normal 
#  scores test.
#
#  To understand the normal scores confidence interval method refer to:
#
#      O'Gorman, T. W. (2020) A method to reduce the width of 
#      confidence intervals while maintaining their coverage 
#      probability. (Under review.)  
#
#  To understand the normal scores test of significance refer to:
#               
#      O'Gorman, T. W. (2020) A generalized normal scores test that
#      increases the power of a test of signifiance for a 
#      coefficient in a linear model. (Under review.)                           
#
#  The confidence limits are found by calling the function
#
#       nscilimits <- function(dataset, complete, details=1)
#
#  The function arguments are:
#
#    1)  The variable "dataset" is a file that is in the form of an R data 
#        frame that includes all of the variables that will be used in the
#        analysis. It will be used in a read.table() function to 
#        create an R data frame.
#    2)  "complete" is a character string that specifies the 
#        complete (full) model, including the confidence interval variable.
#        This model uses the same syntax as the lm() function in R.
#        Note: The confidence interval will be computed for the last
#              variable in the "complete" character string.
#    3)  if details = 0 the limits will be returned by the function, but
#        no output will be printed.
#        if details = 1 (default) the basic results will be printed
#        and the limits will be returned.
#        if details = 2 more information about the search will be printed
#        and the limits will be returned.
#


#  Notes:
#
#    1) The data frame cannot contain missing values for any variables
#       used in the complete model.
#    2) This function calls the invert, ns.1tailed, and trad.1tailed 
#       functions, which are included.
#    3) This function is written in base R.  No packages are required. 
#
#  Examples:
#
#    If blood pressure data is on a dataset(bpdata) that has
#    blood pressure (bp), age (age), and a treatment indicator (group),
#    we could find the normal scores 95% confidence interval for the
#    treatment effect by using this code:
#
#      source("nsci.r")
#      bplimits <- nscilimits(dataset=bpdata, complete=c("bp~group") )
#
#    After this code is executed the bplimits will contain
#    the lower and upper limits.
#
#    We could expand this example if we needed to include age as a
#    covariate.
#
#      source("nsci.r")
#      bplimits <- nscilimits(dataset=bpdata, complete=c("bp~age+group"))
#
#  Note that the group variable was specified as the last variable in the
#  model because we wanted the confidence interval for the group effect.
#
#  These R functions were carefully checked and I believe
#  that the functions are correct.  However, the author is not
#  responsible for any errors that may still exist in the code.
#
#  Please report any issues concerning this code to T. W. O'Gorman via 
#  email at twogorman@gmail.com
#
# *******************************************************
#
nscilimits <- function(dataset, complete, details = 1)         {

if(details == 2) {
  cat("\n"," Data set:   ", dataset, "\n\n")
                  }
  ci.df  <- read.table(dataset) 
if(details == 2) {
  print(ci.df)
                   }
attach(ci.df)

vars <- strsplit(complete,"~")[[1]]
depvar <- vars[[1]]
vars <- strsplit(complete,"+", fixed = TRUE)[[1]]
nvars <- length(vars)

if(nvars == 1) {
  vars    <- strsplit(complete,"~",fixed=TRUE)[[1]]
  civar   <- vars[[2]]
  redvec  <- c(depvar,"1")
  reduced <- paste(redvec,collapse="~")
               }

if(nvars > 1){
  civar     <- vars[[nvars]]
  nreduced  <- nvars - 1
  redvec    <- vars[1:nreduced]
  reduced   <- paste(redvec,collapse="+")
               }
depvar  <- gsub(" ","",depvar)
reduced <- gsub(" ","",reduced)
civar   <- gsub(" ","",civar)
n   <- length(ci.df[,depvar])
limits <- double(2)

compmodel <- lm(as.formula(complete), data = ci.df)
if(details == 2) {
  cat("\n", "  Summary of the linear model using the raw data ", "\n")
  print(summary(compmodel))
                  }
beta  <- summary(compmodel)$coefficients[civar,1]
se    <- summary(compmodel)$coefficients[civar,2]
tstat <- summary(compmodel)$coefficients[civar,3]
ndf   <- df.residual(compmodel)

tradlow95 <- beta + qt(0.025, ndf)*se
tradhi95  <- beta + qt(0.975, ndf)*se

estlow <- beta + qt(0.01, ndf)*se
esthi  <- beta + qt(0.04, ndf)*se
pstar  <- 0.975

#            Compute lower limit

lowerlmt <- invert(ci.df, complete, reduced, civar, depvar, pstar, estlow, esthi, n, details)

estlow <- beta + qt(0.96, ndf)*se
esthi  <- beta + qt(0.99, ndf)*se
pstar  <- 0.025

#            Compute upper limit

upperlmt <- invert(ci.df, complete, reduced, civar, depvar, pstar, estlow, esthi, n, details)

limits <- c(lowerlmt, upperlmt)
if(details >= 1) {
cat("\n","*************************************************************","\n","\n")
  cat("  Input:","\n","\n")
  cat("    Function arguments for the nscilimits function:","\n")
  cat("    Data set:", dataset,"\n")
  cat("    Complete model: ",complete,"\n")
  cat("    Reduced model : ",reduced,"\n","\n")
  cat("  Output:","\n","\n")
  cat("    Normal Scores  95% Confidence Interval Limits:","\n")
  cat("    Lower limit: ", lowerlmt, " Upper limit: ", upperlmt, "\n\n")
  cat("    Traditional (Textbook) 95% Confidence Interval Limits:","\n")
  cat("    Trad. Lower limit: ", tradlow95, " Trad. Upper limit: ", tradhi95, "\n")
 
  cat("\n","*************************************************************","\n","\n")
                 } 
return(limits)                                                            }

# *************************************************************

invert <- function(ci.df, complete, reduced, civar, depvar, pstar, estlow, esthi, n, details) {

ciadj.df <- ci.df

ns.df <- ciadj.df

repeat {
  ciadj.df[,depvar] <- ci.df[ ,depvar] - estlow*ci.df[ ,civar]

plow <- ns.1tailed(ciadj.df, ns.df, complete, reduced, civar, depvar, n, details)

ciadj.df[,depvar] <- ci.df[ ,depvar] - esthi*ci.df[ ,civar]
phi  <- ns.1tailed(ciadj.df, ns.df, complete, reduced, civar, depvar, n, details)
if(plow >= pstar & phi <= pstar) break

shift <-(esthi - estlow)/2

if(phi < pstar) {
  estlow <- estlow - shift
  esthi  <- esthi  - shift
                 }
if(plow > pstar) {
  estlow <- estlow + shift
  esthi  <- esthi  + shift
                  }
         }
#
#    Interval found
#

if(details == 2){ 
cat(" Now found initial limits ", estlow, esthi, " begin 12 bisections ","\n")}

  nupdates <- 1

  repeat                        {
 estnext <- (estlow + esthi)/2
ciadj.df[,depvar] <- ci.df[ ,depvar] - estnext*ci.df[ ,civar]
pnext <- ns.1tailed(ciadj.df, ns.df, complete, reduced, civar, depvar, n, details)

if(nupdates == 12)  break

if(pnext < pstar) {
   esthi  <- estnext
   phi    <- pnext
  } else {
   estlow <- estnext
   plow   <- pnext
         }

nupdates <- nupdates + 1
                                }
return(estnext)
                                                             }

# *********************************************************

ns.1tailed <- function(ciadj.df, ns.df, complete, reduced, civar, depvar, n, details) {

red <- lm(as.formula(reduced), data=ciadj.df)
resid <- residuals(red)
predicted <- predict(red)
#                               compute traditional quantiles
resid25 <- quantile(resid,0.25,type=6)
resid50 <- quantile(resid,0.50,type=6)
resid75 <- quantile(resid,0.75,type=6)
sigmaresid <- (resid75-resid25)/1.349
rankresid <- rank(resid)
lowerprob <- rankresid/(n+1)
ns <- qnorm(lowerprob)

ns25 <- quantile(ns, 0.25, type = 6)
ns50 <- quantile(ns, 0.50, type = 6)
ns75 <- quantile(ns, 0.75, type = 6)
sigmans <- (ns75 - ns25)/1.349

standardns <- (ns- ns50)*(sigmaresid/sigmans) + resid50
ns.df[[depvar]] <- predicted + standardns
pns <- trad.1tailed(ns.df,complete, details)

return(pns)                                                  

# ********************************************************************

trad.1tailed <- function(test.df, complete, details)  {
comp.model <- lm(as.formula(complete), data = test.df)
lmsummary <-summary(comp.model)
nvars <-length(coef(lmsummary)[,3])
p <- pt(coef(lmsummary)[nvars,3],df.residual(comp.model),lower=TRUE)
return(p)
                                                  }
# *********************************************************************

