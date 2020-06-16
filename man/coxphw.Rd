\name{coxphw}

\alias{coxphw}
%%\alias{coxphw.fit}

\title{ Weighted Estimation in Cox Regression }

\description{ Weighted Cox regression as proposed by Schemper et al. (2009) \doi{10.1002/sim.3623}
provides unbiased estimates of average hazard ratios also in case of non-proportional hazards. Time-dependent
effects can be conveniently estimated by including interactions of covariates with arbitrary
functions of time, with or without making use of the weighting option. }
%Additionally estimation of nonlinear effects using fractional polynomials similar to the MFP algorithm
%is provided. This feature can also be used to estimate the interaction of a covariate with
%a nonlinear function of time.

\usage{
coxphw(formula, data, template = c("AHR", "ARE", "PH"), subset, na.action,
       robust = TRUE, jack = FALSE, betafix = NULL, alpha = 0.05,
       trunc.weights = 1, control, caseweights, x = TRUE, y = TRUE,
       verbose = FALSE, sorted = FALSE, id = NULL, clusterid = NULL, \dots)
}

\arguments{
  \item{formula}{ a formula object with the response on the left of the  operator and the model terms
                  on the right. The response must be a survival object as returned by
                  \code{\link[survival:Surv]{Surv}}. \code{formula} may include \code{offset}-terms or
                  functions of time (see example).}
  \item{data}{ a data frame in which to interpret the variables named in \code{formula}. }
  \item{template}{ choose among three pre-defined templates: \code{"AHR"} requests estimation of
                   average hazard ratios (Schemper et al., 2009), \code{"ARE"} requests estimation of
                   average regression effects (Xu and O'Quigley, 2000) and \code{"PH"} requests
                   Cox proportional hazards regression.
                   Recommended and default template is \code{"AHR"}. }
  \item{subset}{expression indicating which subset of the rows of data should be used in the fit.
                All observations are included by default.}
  \item{na.action}{missing-data filtering. Defaults to \code{options$na.action}. Applied after subsetting data,
                   but applied to the all variables in the data set (not only those listed in the formula).}
  \item{robust}{ if set to TRUE, the robust covariance estimate (Lin-Wei) is used; otherwise the Lin-Sasieni
                 covariance estimate is applied. Default is TRUE. }
  \item{jack}{ if set to TRUE, the variance is based on a complete jackknife. Each individual (as
               identified by \code{id}) is left out in turn. The resulting matrix of DFBETA residuals
               D is then used to compute the variance matrix: V = D'D. Default is FALSE. }
  \item{betafix}{ can be used to restrict the estimation of one or more regression coefficients to pre-defined
                  values. A vector with one element for each model term as given in \code{formula} is expected
                  (with an identical order as in \code{formula}). If estimation of a model term is requested,
                  then the corresponding element in \code{betafix} has to be set to \code{NA}, otherwise it
                  should be set to the fixed parameter value. The default value is \code{betafix = NULL},
                  yielding unrestricted estimation of all regression coefficients.}
  \item{alpha}{ the significance level (1-\eqn{\alpha}), 0.05 as default. }
  \item{trunc.weights}{ specifies a quantile at which the (combined normalized) weights are to be truncated.
                        It can be used to increase the precision of the estimates, particularly if
                        \code{template} \code{= "AHR"} or \code{"ARE"} is used. Default is 1 (no truncation).
                        Recommended value is 0.95 for mild truncation. }
  \item{control}{ Object of class \code{\link{coxphw.control}} specifying iteration limit and other
                  control options. Default is \code{coxphw.control(...)}.  }
  \item{caseweights}{ vector of case weights, equivalent to \code{weights} in \code{\link[survival:coxph]{coxph}}.
                      If \code{caseweights} is a vector of integers, then the estimated coefficients
                      are equivalent to estimating the model from data with the individual cases replicated as
                      many times as indicated by \code{caseweights}. These weights should not be confused with the
                      weights in weighted Cox regression which account for the non-proportional hazards. }
  \item{x}{ requests copying explanatory variables into the output object. Default is TRUE. }
  \item{y}{ requests copying survival information into the output object. Default is TRUE. }
  \item{verbose}{ requests echoing of intermediate results. Default is FALSE. }
  % (particularly, for estimation with  \code{fp()}-terms).
  \item{sorted}{ if set to TRUE, the data set will not be sorted prior to passing it to FORTRAN.
                 This may speed up computations. Default is FALSE. }
  \item{id}{ a vector of subject identification integer numbers starting from 1 used only if the data are
             in the counting process format. These IDs are used to compute the robust covariance matrix.
             If \code{id = NA} (the default) the program assumes that each line of the data set refers to
             a distinct subject. }
  \item{clusterid}{a vector of cluster identification integer numbers starting from 1. These IDs are used
                   to compute the robust covariance matrix. If \code{clusterid = NA} (the default) the program
                   assumes that no cluster exist. }
  \item{\dots}{ additional arguments. }
}
%\item{alpha.fp}{ a vector of length three which specifies p-values to include a fractional
%                   polynomial \code{fp()}-term:
%                   \code{alpha.fp[1]} defines the threshold p-value for keeping at least the linear
%                   term in the model, \code{alpha.fp[2]} defines the threshold p-value for the nonlinear
%                   term, and \code{alpha.fp[3]} defines the threshold p-value for the complexity of
%                   the nonlinear term (either a fractional polynomial of degree 1 or 2). Default is
%                   \code{c(0.20, 0.05, 0.05)}.}
%  \item{fp.max}{ highest power of \code{fp()}-terms applied to all \code{fp()}-terms. Select 1 or 2,
%                 with 2 as default. }

\details{ If Cox's proportional hazards regression is used in the presence of non-proportional
hazards, i.e., with underlying time-dependent hazard ratios of prognostic factors, the average relative
risk for such a factor is under- or overestimated and testing power for the corresponding regression
parameter is reduced. In such a situation weighted estimation provides a parsimonious
alternative to more elaborate modelling of time-dependent effects. Weighted estimation in Cox regression
extends the tests by Breslow and Prentice to a multi-covariate situation as does the Cox model to Mantel's
logrank test. Weighted Cox regression can also be seen as a robust alternative to the standard Cox
estimator, reducing the influence of outlying survival times on parameter estimates.

%%Weighted estimation assigns weights to risk sets, according to the survival function estimates
%%(Prentice weights), the number of subjects at risk (Breslow weights), or according to their square
%%roots (Tarone-Ware weights). These weights are applied to the summands of the score function. The final
%%estimate is the vector of parameter values which equates the score function to 0. Since there is one
%%score function corresponding to each parameter of the model, weights may be applied to some but not
%%necessarily to all parameters of a model.

Three pre-defined templates can be requested: \cr 1) \code{"AHR"}, i.e., estimation of average hazard ratios
(Schemper et al., 2009) using Prentice weights with censoring correction and robust variance estimation; \cr
2) \code{"ARE"}, i.e., estimation of average regression effects (Xu and O'Quigley, 2000) using censoring
correction and robust variance estimation; or \cr 3) \code{"PH"}, i.e., Cox proportional hazards regression
using robust variance estimation.

Breslow's tie-handling method is used by the program, other methods to handle ties are currently not
available.

A fit of \code{coxphw} with \code{template = "PH"} will yield identical estimates as a fit of
\code{\link[survival:coxph]{coxph}} using Breslow's tie handling method and robust variance estimation
(using \code{\link[survival:cluster]{cluster}}).

If \code{robust = FALSE}, the program estimates the covariance matrix using the Lin (1991)
and Sasieni (1993) sandwich estimate \eqn{A^{-1}BA^{-1}} with \eqn{-A} and \eqn{-B} denoting the sum
of contributions to the second derivative of the log likelihood, weighted by \eqn{w(t_j)} and \eqn{w(t_j)^2},
respectively. This estimate is independent from the scaling of the weights and reduces to the inverse
of the information matrix in case of no weighting. However, it is theoretically valid only in case of
proportional hazards. Therefore, since application of weighted Cox regression usually implies a violated
proportional hazards assumption, the robust Lin-Wei covariance estimate is used by default (\code{robust = }
\code{TRUE}).

If some regression coefficients are held constant using \code{betafix}, no standard errors
are given for these coefficients as they are not estimated in the model. The global Wald test
only relates to those variables for which regression coefficients were estimated.

An \code{\link[stats:offset]{offset}} term can be included in the \code{formula} of \code{coxphw}.
In this way a variable can be specified which is included in the model but its parameter estimate is
fixed at 1.
}

%Estimation of nonlinear effects of continuous variables using \emph{fractional polynomials} similar to
%the MFP algoritm (Royston and Altman, 1994; Royston and Sauerbrei, 2008) is available. (It is similar,
%but not identical to \code{\link[mfp:mfp]{mfp}}, e.g., tests are based on robust Wald test statistics
%instead of likelihood ratio statistics.) Estimation of fractional polynomials can also be applied to a
%survival time in a time-by-covariate interaction. A pretransformation
%(by shifting to positive values and dividing by scaling factor [a power of 10] such that the standard
%deviation is approximately equal to one) may be applied to continuous variables for which fractional
%polynomials are requested. For examples see the documentation of \code{\link{plotshape}}.

\value{ A list with the following components:
 \item{coefficients}{the parameter estimates.}
 \item{var}{the estimated covariance matrix.}
 \item{df}{the degrees of freedom.}
 \item{ci.lower}{the lower confidence limits of exp(beta).}
 \item{ci.upper}{the upper confidence limits of exp(beta).}
 \item{prob}{the p-values.}
 \item{linear.predictors}{the linear predictors.}
 \item{n}{the number of observations.}
 \item{dfbeta.resid}{ matrix of DFBETA residuals. }
 \item{iter}{the number of iterations needed to converge.}
 \item{method.ties}{the ties handling method.}
 \item{PTcoefs}{ matrix with scale and shift used for pretransformation of \code{fp()}-terms. }
% \item{fp.ind}{ if \code{fp()}-terms were used, a matrix with all transformations of fractional
%               polynomials used in the fitting process. }
% \item{ind}{ if \code{fp()}-terms were used, an indicator which variables from all transformations of
%             fractional polynomials used in the fitting process were seleted in the final model. }
 \item{cov.j}{the covariance matrix computed by the jackknife method (only computed if \code{jack = TRUE}).}
 \item{cov.lw}{the covariance matrix computed by the Lin-Wei method (robust covariance)}
 \item{cov.ls}{the covariance matrix computed by the Lin-Sasieni method.}
 \item{cov.method}{the method used to compute the (displayed) covariance matrix and the standard errors.
       This method is either "jack" if \code{jack = TRUE}, or "Lin-Wei" if \code{jack = FALSE}. }
 \item{w.matrix}{a matrix with four columns according to the number of uncensored failure
                 times. The first column contains the failure times, the remaining columns (labeled
                 \code{w.raw}, \code{w.obskm}, and \code{w}) contain the raw weights, the weights
                 according to the inverse of the Kaplan-Meier estimates with reverse status indicator
                 and the normalized product of both.}
 \item{caseweights}{ if \code{x = TRUE} the case weights. }
 \item{Wald}{ Wald-test statistics. }
 \item{means}{the means of the covariates.}
 \item{offset.values}{offset values. }
 \item{dataline}{the first dataline of the input data set (required for \code{plotfp}).}
 \item{x}{ if \code{x = TRUE} the explanatory variables. }
 \item{y}{the response.}
 \item{alpha}{the significance level = 1 - confidence level.}
 \item{template}{ the requested template. }
 \item{formula}{the model formula.}
 \item{betafix}{the \code{betafix} vector.}
 \item{call}{the function call.}
 %%\item{cards}{ . }
 %%\item{parms}{ . }
 %%\item{exit.code}{ . }
 %%\item{ioarray}{ . }
 %%\item{offset}{ . }
 %%\item{method}{the estimation method (usually weighted estimation).}
 %%\item{method.ci}{the confidence interval estimation method (profile likelihood or Wald).}
 %% \item{loglik}{the null and maximimized (penalized) log likelihood. - DAS GIBT ES NICHT !! }
}

\note{ The \code{SAS} macro \code{WCM} with similar functionality is offered for download at
       \url{http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/wcmcoxphw/} .

       Up to Version 2.13 \pkg{coxphw} used a slightly different syntax (arguments: \code{AHR},
       \code{AHR.norobust}, \code{ARE}, \code{PH}, \code{normalize}, \code{censcorr}, \code{prentice},
       \code{breslow}, \code{taroneware}). From Version 3.0.0 on the old syntax is disabled.
       From Version 4.0.0 estimation of fractional polynomials is disabled. }

\references{Dunkler D, Ploner M, Schemper M, Heinze G. (2018) Weighted Cox Regression Using the {R} Package {coxphw}. \emph{JSS} \bold{84}, 1--26, \doi{10.18637/jss.v084.i02}.

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

Schemper M (1992).  Cox Analysis of Survival Data with Non-Proportional Hazard Functions. \emph{J  R STAT SOC D}
\bold{41}, 455-465.

Schemper M, Wakounig S and Heinze G (2009). The Estimation of Average Hazard Ratios by Weighted Cox
Regression. \emph{STAT MED} \bold{28}, 2473-2489. \doi{10.1002/sim.3623}

Xu R and O'Quigley J (2000). Estimating Average Regression Effect Under Non-Proportional Hazards.
\emph{Biostatistics} \bold{1}, 423-439. }

\author{Georg Heinze, Meinhard Ploner, Daniela Dunkler}

\seealso{ \code{\link{concord}}, \code{\link{plot.coxphw}}, \code{\link{predict.coxphw}}, \code{\link{plot.coxphw.predict}}, \code{\link[survival:coxph]{coxph} }}

\examples{
data("gastric")

# weighted estimation of average hazard ratio
fit1 <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = "AHR")
summary(fit1)
fit1$cov.lw     # robust covariance
fit1$cov.ls     # Lin-Sasieni covariance


# unweighted estimation, include interaction with years
# ('radiation' must be included in formula!)
gastric$years <- gastric$time / 365.25
fit2 <- coxphw(Surv(years, status) ~ radiation + years : radiation, data = gastric,
               template = "PH")
summary(fit2)


# unweighted estimation with a function of time
data("gastric")
gastric$yrs <- gastric$time / 365.25

fun <- function(t) { (t > 1) * 1 }
fit3 <- coxphw(Surv(yrs, status) ~ radiation + fun(yrs):radiation, data = gastric,
               template = "PH")

# for more examples see vignette or predict.coxphw
}
%# does not work any longer
%# select best fp(1) for interaction with time using fp-functionality
%# (don't include 'radiation' in formula if using fp(time)!)
%%#fit3 <- coxphw(Surv(time, status) ~ fp(time) : radiation, data = gastric, template = "PH",
%#               alpha.fp = c(1, 0.05, 0.05), fp.max = 1)
%#summary(fit3)

\keyword{ survival }
\keyword{ regression }
