\name{plot.coxphw.predict}

\alias{plot.coxphw.predict}

\title{ Plot the Relative or Log Relative Hazard Versus Values of a Continuous Covariable.}

\description{ This function visualizes a nonlinear or a time-dependent effect of a \code{predict.coxphw}
              object. }
              %(possibly estimated by the fractional polynomial functionality of \code{coxphw})

\usage{
\method{plot}{coxphw.predict}(x, addci = TRUE,  xlab = NULL, ylab = NULL, \dots)
}

\arguments{
  \item{x}{ an output object of \code{coxphw}.} %, which may include \code{fp()}-terms. }
  \item{addci}{ confidence intervalls are obtained. Default is TRUE. }
  \item{xlab}{ label for x-axis of plot, uses variable specified in \code{x} in \code{predict} as default. }
  \item{ylab}{ label for y-axis of plot, uses as appropriate either "relative hazard" or "log relative hazard" as default. }
  \item{\dots}{ further parameters, to be used for plots (e.g., scaling of axes). }
}

\details{ This function can be used to depict the estimated nonlinear or time-dependent
          effect of an object of class \code{predict.coxphw}. It supports simple nonlinear
          effects as well as interaction effects of continuous variables with binary
          covariates (see examples section in \code{\link{predict.coxphw}}. ).
%          \code{plot.x} must include only positive numbers greater than 0.
}

\value{ No output value. }

%\value{ A list with the follwing components:
% \item{xbeta}{ a matrix with estimates of (log) relative hazard. }
% \item{cilower}{ lower confidence limits. }
% \item{ciupper}{ upper confidence limits. }
% \item{p}{ p-value, only if \code{plot = FALSE} and \code{pval = TRUE}. }
% \item{alpha}{ the significance level = 1 - confidence level. }
% \item{plot.x}{ the data values of the continuous variable. }
% \item{exp}{ an indicator if log relative hazard or relative hazard was obtained. }
%}

\note{ In \pkg{coxphw} version 4.0.0 the old \code{plotshape} function is replaced with
       \code{predict.coxphw} and \code{plot.coxphw.predict}. }

\author{Georg Heinze, Meinhard Ploner, Daniela Dunkler}

\references{ Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous
Covariates: Parsimonious Parametric Modelling. \emph{Applied StatisticsJ R STAT SOC C-APPL} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-building. A pragmatic approach to regression
analysis based on fractional polynomials for modelling continuous variables.} Wiley, Chichester, UK. }

\seealso{ \code{\link{coxphw}}, \code{\link{predict.coxphw}} }

\examples{
set.seed(30091)
n <- 300
x <- 1:n
true.func <- function(x) 3 * (x / 100)^{2} - log(x / 100) - 3 * x / 100
x <- round(rnorm(n = x) * 10 + 40, digits = 0)
time <- rexp(n = n, rate = 1) / exp(true.func(x))
event <- rep(x = 1, times = n)
futime <- runif(n = n, min = 0, max = 309000)
event <- (time < futime) + 0
time[event == 0] <- futime[event == 0]
my.data <- data.frame(x, time, event)

fitahr <- coxphw(Surv(time, event) ~ x, data = my.data, template = "AHR")

# estimated function
plotx <- quantile(x, probs = 0.05):quantile(x, probs = 0.95)
plot(predict(fitahr, type = "shape", newx = plotx, refx = median(x), x = "x",
             verbose = FALSE))

# true function
lines(x = plotx, true.func(plotx) - true.func(median(plotx)), lty = 2)

legend("topright", legend=c("AHR estimates", "true"), bty = "n", lty = 1:2, inset = 0.05)

# for more examples see predict.coxphw
}

\keyword{survival}
