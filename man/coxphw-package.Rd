\name{coxphw-package}

\alias{coxphw-package}

\docType{package}

\title{ Weighted Estimation in Cox Regression }

\description{This package implements weighted estimation in Cox regression as proposed by
             Schemper, Wakounig and Heinze (Statistics in Medicine, 2009, \doi{10.1002/sim.3623}). Weighted Cox regression
             provides unbiased average hazard ratio estimates also in case of non-proportional hazards.
             The package provides options to estimate time-dependent effects conveniently by
             including interactions of covariates with arbitrary functions of time, with or without
             making use of the weighting option. For more details we refer to Dunkler, Ploner, Schemper and Heinze
             (Journal of Statistical Software, 2018, \doi{10.18637/jss.v084.i02}).  }
             % Additionally estimation of nonlinear effects using fractional polynomials similar to the
             % MFP algorithm (Royston, Sauerbrei, 2008) is provided. This feature can also be used to
             % estimate the interaction of a covariate with a nonlinear function of time.

\details{
\tabular{ll}{
Package: \tab \code{coxphw}\cr
Type: \tab Package\cr
Version: \tab 4.0.1\cr
Date: \tab 2018-04-16\cr
License: \tab GPL-2\cr
}

Main functions included in the \pkg{coxphw} package are
\tabular{ll}{
\code{\link{coxphw}}    \tab weigthed estimation of Cox regression: either (recommended) estimation of \cr
                        \tab average hazard ratios (Schemper et al., 2009), estimation of average regression \cr
                        \tab effects (Xu and O'Quigley, 2000), or proportional hazards regression. \cr
                        \tab \cr
\code{\link{plot}}      \tab plots the weights used in a weighted Cox regression analysis against time. \cr
\code{\link{concord}}   \tab obtains generalized concordance probabilities with confidence intervalls. \cr
%\code{\link{fp.power}}  \tab provides fractional polynomials as accessible function. \cr
\code{\link{predict}}   \tab obtains the effect estimates (of e.g. a nonlinear or a time-dependent effect)\cr
                        \tab at specified values of a continuous covariable. With \code{\link{plot.coxphw.predict}}\cr
                        \tab these relative or log relative hazard versus values of the continuous covariable \cr
                        \tab can be plotted. \cr
\code{\link{wald}}      \tab obtain Wald chi-squared test statistics and p-values for one or more regression \cr
                        \tab coefficients given their variance-covariance matrix. \cr
}

Data sets included in the \pkg{coxphw} package are
\tabular{ll}{
\code{\link{biofeedback}} \tab biofeedback treatment data \cr
\code{\link{gastric}}     \tab gastric cancer data \cr
}
}

\note{ The \code{SAS} macro \code{WCM} with similar functionality can be obtained at \cr
       \url{http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/wcmcoxphw/}.

       Important version changes: \cr
       Up to Version 2.13 \pkg{coxphw} used a slightly different syntax (arguments: \code{AHR}, \code{AHR.norobust},
       \code{ARE}, \code{PH}, \code{normalize}, \code{censcorr}, \code{prentice}, \code{breslow},
       \code{taroneware}).
       From Version 3.0.0 on the old syntax is disabled.
       From Version 4.0.0 fractional polynomials are disabled and \code{plotshape} is replaced with \code{predict}
       and \code{plot.coxphw.predict}. }

\author{Georg Heinze, Meinhard Ploner, Daniela Dunkler \cr
        {Maintainer: \email{georg.heinze@meduniwien.ac.at}} }

\references{Dunkler D, Ploner M, Schemper M, Heinze G. (2018) Weighted Cox Regression Using the {R} Package {coxphw}. \emph{JSS} \bold{84}, 1--26, \doi{10.18637/jss.v084.i02}.

Dunkler D, Schemper M, Heinze G. (2010) Gene Selection in Microarray Survival Studies Under
Possibly Non-Proportional Hazards. \emph{Bioinformatics} \bold{26}:784-90.

Lin D and Wei L (1989). The Robust Inference for the Cox Proportional Hazards Model.
\emph{J AM STAT ASSOC} \bold{84}, 1074-1078.

Lin D (1991). Goodness-of-Fit Analysis for the Cox Regression Model Based on a Class of Parameter
Estimators. \emph{J AM STAT ASSOC} \bold{86}, 725-728.

Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous Covariates:
Parsimonious Parametric Modelling. \emph{J R STAT SOC C-APPL} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-Building. A Pragmatic Approach to Regression
Analysis Based on Fractional Polynomials for Modelling Continuous Variables.} Wiley, Chichester, UK.

Sasieni P (1993). Maximum Weighted Partial Likelihood Estimators for the Cox Model. \emph{J AM STAT ASSOC}
\bold{88}, 144-152.

Schemper M (1992).  Cox Analysis of Survival Data with Non-Proportional Hazard Functions. \emph{J R STAT SOC D}
\bold{41}, 455-465.

Schemper M, Wakounig S and Heinze G (2009). The Estimation of Average Hazard Ratios by Weighted Cox
Regression. \emph{Stat Med} \bold{28}, 2473-2489, \doi{10.1002/sim.3623}.

Xu R and O'Quigley J (2000). Estimating Average Regression Effect Under Non-Proportional Hazards.
\emph{Biostatistics} \bold{1}, 423-439. }

\seealso{ \code{\link{coxphw}}, \code{\link{concord}}, \code{\link{plot.coxphw}}, \code{\link{predict.coxphw}}, \code{\link{plot.coxphw.predict}}, \code{\link{wald}} }
%\code{\link{fp.power}},

\examples{
## for examples see coxphw
}

\keyword{ package }
\keyword{ regression }
\keyword{ survival }
