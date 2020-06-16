\name{summary.coxphw}

\alias{summary.coxphw}

\title{ Summary Method for Objects of Class \code{coxphw} }

\description{ This class of objects is returned by the \code{coxphw} function. 
  Objects of this class have methods for the functions \code{summary}, \code{print}, \code{coef}, 
  \code{vcov}, \code{plot}, and \code{confint}. }

\usage{ \method{summary}{coxphw}(object, print = TRUE, ...) }

\arguments{
  \item{object}{ object of class \code{coxphw}. }
  \item{print}{ print summary. Default is TRUE. }
  \item{\dots}{ further arguments. }
}

\details{ If some regression coefficients are held fixed by \code{betafix}, no standard errors and related quantities
are given for these coefficients as they are not estimated in the model. The global Wald test 
only relates to those variables for which regression coefficients were estimated.}

%%\value{ ~Describe the value returned }

%%\references{~put references to the literature/web site here ~ }

\author{ Georg Heinze, Daniela Dunkler }

%%\note{  ~~further notes~~ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{coxphw}} }

%%\examples{ }

\keyword{ utilities }
