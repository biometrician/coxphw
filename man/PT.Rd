\name{PT}

\alias{PT}

\title{ Pretransformation function }

\description{ Provides automatic pretransformation of variables (to well-scaled and nonzero values). }

\usage{ PT(z) }

\arguments{
  \item{z}{a vector of numerical values. }
}

\details{ The function transforms a variable by shifting to positive values, and dividing by scaling 
          factor (a power of 10) such that the standard deviation is approximately equal to 1. }

\value{ (\code{z} + shift) / scale }

\author{Georg Heinze}

\seealso{ \code{\link{coxphw}} }

\keyword{math}

\examples{
PT(z = c(-6, -1, 4, 6))
}
