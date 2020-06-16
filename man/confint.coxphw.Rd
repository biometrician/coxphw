\name{confint.coxphw}

\alias{confint.coxphw}

\title{ Confidence Intervals for Model Parameters }

\description{ Computes confidence intervals for one or more parameters in a model fitted by 
              \code{coxphw}. Objects of this class have methods for the functions \code{summary}, 
              \code{print}, \code{coef}, \code{vcov}, \code{plot}, and \code{confint}.}

\usage{ \method{confint}{coxphw}(object, parm, level = 0.95, ...) }

\arguments{
  \item{object}{ a fitted model object of class \code{coxphw}. }
  \item{parm}{ a specification of which parameters are to be given confidence intervals, either a 
               vector of numbers or a vector of names. If missing, all parameters are considered. }
  \item{level}{ the confidence level required. }
  \item{\dots}{ additional argument(s) for methods. }
}

%%\details{  ~~ If necessary, more details than the description above ~~ }

\value{ A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. 
        These will be labeled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 
        97.5\%). }

%%\references{~put references to the literature/web site here ~ }

\author{ Daniela Dunkler }

%%\note{  ~~further notes~~ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{coxphw}} }

%%\examples{ }

\keyword{ utilities }
