\name{plot.coxphw}

\alias{plot.coxphw}

\title{ Plot Weights of Weighted Estimation in Cox Regression }

\description{ This function plots the weights used in a weighted Cox regression analysis against time. }

\usage{ \method{plot}{coxphw}(x, rank = FALSE, log = FALSE, legendxy = NULL,...) }

\arguments{
  \item{x}{ a \code{coxphw} object. }
  \item{rank}{ if set to TRUE, plots the weights against ranked time. Default is FALSE. }
  \item{log}{ if set to TRUE, shows logarithm of weights. Default is FALSE. }
  \item{legendxy}{ an optional vector of length 2 of the x and y co-ordinates to be used to position the legend. }
  \item{...}{ additional arguments for plotting }
}

\details{ The function plots the survival weights, i.e., the left-continuous survivor function
estimates, the censoring weights, i.e., estimates of the follow-up distribution obtained by
Kaplan-Meier estimation with reversed meaning of the status indicator and the combined
normalized weights, i.e. the product of the survival and the censoring weights, rescaled to a mean of 1.}

\value{ No output value. }

\author{Georg Heinze, Daniela Dunkler}

\note{ In \pkg{coxphw} version 4.0.0 the new \code{plot} function replaces the old \code{plotw} function. }

\seealso{ \code{\link{coxphw}} }

\keyword{survival}

\examples{
data("gastric")

# weighted estimation of average hazard ratio
fit1 <- coxphw(Surv(time, status) ~ radiation, data = gastric,  template = "AHR")
plot(fit1)

# estimation of average regression effect by inverse probability of censoring weights;
# truncate weights at 95th percentile
fit2 <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = "ARE",
               trunc.weights = 0.95)
plot(fit2)
}
