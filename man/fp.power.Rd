\name{fp.power}

\alias{fp.power}

\title{ Provides Fractional Polynomials as Accessible Function }

\description{ Provides fractional polynomials as accessible function. }

\usage{ fp.power(z, a, b = NULL)  }

\arguments{
  \item{z}{a scalar or vector of positive numerical values. }
  \item{a}{first power. }
  \item{b}{optional second power. }
}

\details{ The function returns fp(\code{a}) of \code{z} (and optionally \code{fp(b)} of \code{z}). }                           

\value{ A matrix with one or two columns (if a second power \code{b} was specified), and number of 
        rows equal to the length of \code{z}. The columns are sorted by descending power. }
 
\author{Georg Heinze}

\seealso{ \code{\link{coxphw}} }

\references{ Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous Covariates: Parsimonious Parametric Modelling. \emph{J R STAT SOC C-APPL} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-Building. A Pragmatic Approach to Regression 
Analysis Based on Fractional Polynomials for Modelling Continuous Variables.} Wiley, Chichester, UK. }

\examples{
fp.power(z = c(1, 4, 6), a = 1)
fp.power(z = c(1, 4, 6), a = 0.5)
fp.power(z = c(1, 4, 6), a = 0.5, b = 0.5)
fp.power(z = c(1, 4, 6), a = 0, b = 2)
}

\keyword{survival}
\keyword{math}
