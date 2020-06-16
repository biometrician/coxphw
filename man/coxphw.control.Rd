\name{coxphw.control}

\alias{coxphw.control}

\title{Ancillary arguments for controling coxphw fits}

\description{
This is used to set various numeric parameters controlling a Cox model fit using \code{coxphw}.
Typically it would only be used in a call to \code{coxphw}.
}
\usage{
coxphw.control(iter.max = 200, maxhs = 5, xconv = 1e-4, gconv = 1e-4, maxstep = 1,
               round.times.to = 0.00001, add.constant = 0, pc = TRUE, pc.time = TRUE,
               normalize = TRUE) %, fp.iter = 10)
}
\arguments{
  \item{iter.max}{ maximum number of iterations to attempt for convergence. Default is 200. }
  \item{maxhs}{ maximum number of step-halvings per iterations. Default is 5. The increments of the
                parameter vector in one Newton-Rhaphson iteration step are halved, unless the new
                pseudo-likelihood is greater than the old one, maximally doing \code{maxhs} halvings. }
  \item{xconv}{ specifies the maximum allowed change in standardized parameter estimates to declare
                convergence. Default is 0.0001. }
  \item{gconv}{ specifies the maximum allowed change in score function to declare convergence.
                Default is 0.0001. }
  \item{maxstep}{ specifies the maximum change of (standardized) parameter values allowed in one
                  iteration. Default is 1. }
  \item{round.times.to}{ rounds survival times to the nearest number that is a multiple of
                         \code{round.times.to}. This may improve numerical stability if very small
                         survival times are used (mostly in simulations). Default is 0.00001. }
  \item{add.constant}{ this number will be added to all times. It can be used if some survival times
                       are exactly 0. Default is 0. }
  \item{pc}{ if set to TRUE, it will orthogonalize the model matrix to speed up convergence. Default is TRUE. }
             %during \code{fp()}-mode
  \item{pc.time}{ if set to TRUE, it will orthogonalize time-dependent covariates (i.e., interactions of
                  covariates with functions of time) to speed up convergence. Default is TRUE. }
  \item{normalize}{ if set to TRUE, weights are normalized such that their sum is equal to the number of
                    events. This may speed up or enable convergence, but has no consequences on the
                    estimated regression coefficients. }
%  \item{fp.iter}{ number of iterations for the big fractional polynomial-loop, default is 10. This
%                  may be important if several \code{fp()}-terms are included in the model as inclusion
%                  or exclusion of one term based on the threshold p-value may increase or decrease the
%                  significance of others. }
}

\value{A list containing the values of each of the above constants
}
\author{ Daniela Dunkler }
\seealso{\code{\link{coxphw}}
  }
\keyword{survival}
