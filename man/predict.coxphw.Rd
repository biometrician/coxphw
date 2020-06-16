\name{predict.coxphw}
\alias{predict.coxphw}

\title{Predictions for a weigthed Cox model}

\description{ This function obtains the effect estimates (e.g. of a nonlinear or a
              time-dependent effect) at specified values of a continuous
              covariable for a model fitted by \code{\link{coxphw}}. It prints the
              relative or log relative hazard. Additonally, the linear predictor lp
              or the risk score exp(lp) can be obtained. }

\usage{ \method{predict}{coxphw}(object, type = c("shape", "slice.time", "slice.z", "slice.x", "lp", "risk"),
        x = NULL, newx = NA, refx = NA, z = NULL, at = NULL, exp = FALSE,
        se.fit = FALSE, pval = FALSE, digits = 4, verbose = FALSE, \dots)
}

\arguments{
  \item{object}{ an output object of \code{coxphw}. }
  \item{type}{the type of predicted value. Choices are:
    \code{"lp"} for the linear predictors, \cr
     \code{"risk"} for the risk scores \code{exp(lp)} \cr
     \code{"shape"} for visualizing a nonlinear effect of \code{x}, \cr
    \code{"slice.x"}  for slicing an interaction of type \code{fun(x)*z} at values of \code{x},
    \code{"slice.z"} for slicing an interaction of type \code{fun(x)*z} at a value of \code{z}, \cr
    \code{"slice.time"} for slicing a time-by-covariate interaction of type \code{z + fun(time):z} at values of \code{time} \cr
    See details. }
  \item{x}{ name of the continuous or time variable (use "") for \code{type = "shape"},
            \code{"slice.time"}, \code{"slice.x"}, or \code{"slice.z"}. }
  \item{newx}{ the data values for \code{x} for which the effect estimates should be
               obtained (e.g., \code{30:70}) for \code{type = "shape"},
               \code{"slice.time"}, \code{"slice.x"}, or \code{"slice.z"}. }
  \item{refx}{ the reference value for variable \code{x} for \code{type = "shape"} or \code{"slice.z"}.
               The log relative hazard at this value will be 0. (e.g., \code{refx= 50}).}
  \item{z}{ variable which is in interaction with \code{x} (use "") for \code{"slice.time"},
            \code{"slice.x"}, or \code{"slice.z"}. }
  \item{at}{ if \code{type = "slice.z"} at which level ("slice") of \code{z} should the effect estimates of the
              \code{x} be obtained. }
  \item{exp}{ if set to TRUE (default), the log relative hazard is given, otherwise the relative hazard
              is requested for \code{type = "shape"}, \code{"slice.time"}, \code{"slice.x"},
              or \code{"slice.z"}. }
  \item{se.fit}{if set to TRUE, pointwise standard errors are produced for the predictions for
                \code{type = "shape"}, \code{"slice.time"}, \code{"slice.x"}, or \code{"slice.z"}. }
  \item{pval}{ if set to TRUE add Wald-test p-values to effect estimates at values of
               \code{newx} for \code{type = "shape"}, \code{"slice.time"}, \code{"slice.x"},
               or \code{"slice.z"}. Default is set to FALSE. }
  \item{digits}{ number of printed digits. Default is 4. }
  \item{verbose}{ if set to TRUE (default), results are printed.}
  \item{\dots}{ further parameters. }
}

\details{ This function can be used to depict the estimated nonlinear or time-dependent
          effect of an object of class \code{coxphw}. It supports simple nonlinear
          effects as well as interaction effects of a continuous variable with a binary
          covariate or with time (see examples section).

          If the effect estimates of a simple nonlinear effect of \code{x} without
          interaction is requested with \code{type = "shape"}, then \code{x} (the usually
          continuous covariate), \code{refx} (the reference value of \code{x}) and \code{newx}
          (for these values of \code{x} the effect estimates are obtained) must be given.

          If the effect estimates of an interaction of \code{z} with \code{x} are requested
          with \code{type = "slice.x"}, then \code{x} (the usually continuous variable),
          \code{z} (the categorical variable) and \code{newx} (for these values of \code{x}
          the effect estimates are obtained) must be given.

          If the effect estimates of an interaction of \code{z} with \code{x} for one level of \code{z}
          are requested with \code{type = "slice.z"}), then \code{x} (the usually continuous variable),
          \code{z} (the categorical variable), \code{at} (at which level of \code{z}),
          \code{refx} (the reference value of \code{x}), and \code{newx} (for these values of \code{x}
          the effect estimates are obtained) must be given.

          If the effect estimates of an interaction of \code{z} with \code{time} are requested
          with \code{type = "slice.time"}, then \code{x} (the \code{time}), \code{z} (the
          categorical variable) and \code{newx} (for these values of \code{x} the effect
          estimates are obtained) must be given.

          Note that if the model formula contains time-by-covariate interactions, then the linear predictor
          and the risk score are obtained for the failure or censoring time of each subject.
}
%          \code{plot.x} must include only positive numbers greater than 0.

\value{ If \code{type = "shape"}, \code{"slice.time"}, \code{"slice.x"}, or \code{"slice.z"}
        a list with the following components:
 \item{estimates}{ a matrix with estimates of (log) relative hazard and corresponding confidence limits. }
 \item{se}{ pointwise standard errors, only if \code{se.fit = TRUE}. }
 \item{p}{ p-value, only if \code{pval = TRUE}. }
 \item{alpha}{ the significance level = 1 - confidence level. }
 \item{exp}{ an indicator if log relative hazard or relative hazard was obtained. }
 \item{x}{ name of \code{x}. }

If  \code{type = "lp"} or \code{"risk"}, a vector.
}

\note{ In \pkg{coxphw} version 4.0.0 the old \code{plotshape} function is replaced with
       \code{predict.coxphw} and \code{plot.coxphw.predict}. }

\author{Georg Heinze, Meinhard Ploner, Daniela Dunkler}

\references{ Royston P and Altman D (1994). Regression Using Fractional Polynomials of Continuous
Covariates: Parsimonious Parametric Modelling. \emph{Applied Statistics} \bold{43}, 429-467.

Royston P and Sauerbrei W (2008). \emph{Multivariable Model-building. A pragmatic approach to regression
analysis based on fractional polynomials for modelling continuous variables.} Wiley, Chichester, UK. }

\seealso{\code{\link{coxphw}}, \code{\link{plot.coxphw.predict}} }


\examples{
### Example for type = "slice.time"
data("gastric")
gastric$yrs <- gastric$time / 365.25

# check proportional hazards
fitcox <- coxph(Surv(yrs, status) ~ radiation + cluster(id), data = gastric, x = TRUE,
                method = "breslow")
fitcox.ph <- cox.zph(fit = fitcox, transform = "identity")


## compare and visualize linear and log-linear time-dependent effects of radiation
fit1 <- coxphw(Surv(yrs, status) ~ yrs * radiation, data = gastric, template = "PH")
summary(fit1)

predict(fit1, type = "slice.time", x = "yrs", z = "radiation", newx = c(0.5, 1, 2),
        verbose = TRUE, exp = TRUE, pval = TRUE)


fit2 <- coxphw(Surv(yrs, status) ~ log(yrs) * radiation, data = gastric, template = "PH")
summary(fit2)

predict(fit2, type = "slice.time", x = "yrs", z = "radiation", newx = c(0.5, 1, 2),
        verbose = TRUE, exp = TRUE, pval = TRUE)


plotx <- seq(from = quantile(gastric$yrs, probs = 0.05),
             to = quantile(gastric$yrs, probs = 0.95), length = 100)
y1 <- predict(fit1, type = "slice.time", x = "yrs", z = "radiation", newx = plotx)
y2 <- predict(fit2, type = "slice.time", x = "yrs", z = "radiation", newx = plotx)

plot(x = fitcox.ph, se = FALSE, xlim = c(0, 3), las = 1, lty = 3)
abline(a = 0, b = 0, lty = 3)
lines(x = plotx, y = y1$estimates[, "coef"], col = "red", lty = 1, lwd = 2)
lines(x = plotx, y = y2$estimates[, "coef"], col = "blue", lty = 2, lwd = 2)
legend(x = 1.7, y = 1.6, title = "time-dependent effect", title.col = "black",
       legend = c("LOWESS", "linear", "log-linear"), col = c("black", "red", "blue"),
       lty = c(3, 1:2), bty = "n", lwd = 2, text.col = c("black", "red", "blue"))



### Example for type = "shape"
set.seed(512364)
n <- 200
x <- 1:n
true.func <- function(x) 2.5 * log(x) - 2
x <- round(runif(x) * 60 + 10, digits = 0)
time <- round(100000 * rexp(n= n, rate = 1) / exp(true.func(x)), digits = 1)
event <- rep(x = 1, times = n)
my.data <- data.frame(x,time,event)

fit <- coxphw(Surv(time, event) ~ log(x) + x, data = my.data, template = "AHR")

predict(fit, type = "shape", newx = c(30, 50), refx = 40, x = "x", verbose = TRUE)

plotx <- seq(from = quantile(x, probs = 0.05),
             to = quantile(x, probs = 0.95), length = 100)
plot(predict(fit, type = "shape", newx = plotx, refx = 40, x = "x"))


### Example for type = "slice.x" and "slice.z"
set.seed(75315)
n <- 200
trt <- rbinom(n = n, size = 1, prob = 0.5)
x <- 1:n
true.func <- function(x) 2.5 * log(x) - 2
x <- round(runif(n = x) * 60 + 10, digits = 0)
time <- 100 * rexp(n = n, rate = 1) / exp(true.func(x) /
                                      4 * trt - (true.func(x) / 4)^2 * (trt==0))
event <- rep(x = 1, times = n)
my.data <- data.frame(x, trt, time, event)

fun<-function(x) x^(-2)
fit <- coxphw(Surv(time, event) ~ x * trt + fun(x) * trt , data = my.data,
              template = "AHR", verbose = FALSE)

# plots the interaction of trt with x (the effect of trt dependent on the values of x)
plotx <- quantile(x, probs = 0.05):quantile(x, probs = 0.95)
plot(predict(fit, type = "slice.x", x = "x", z = "trt",
             newx = plotx, verbose = FALSE), main = "interaction of trt with x")

# plot the effect of x in subjects with trt = 0
y0 <- predict(fit, type = "slice.z", x = "x", z = "trt", at = 0, newx = plotx,
              refx = median(x), verbose = FALSE)
plot(y0, main = "effect of x in subjects with trt = 0")


# plot the effect of x in subjects with trt = 1
y1 <- predict(fit, type = "slice.z", x = "x", z = "trt", at = 1, newx = plotx,
              refx = median(x), verbose = FALSE)
plot(y1, main = "effect of x in subjects with trt = 1")


# Example for type = "slice.time"
set.seed(23917)
time <- 100 * rexp(n = n, rate = 1) / exp((true.func(x) / 10)^2 / 2000 * trt + trt)
event <- rep(x = 1, times = n)
my.data <- data.frame(x, trt, time, event)
plot.x <- seq(from = 1, to = 100, by = 1)

fun  <- function(t) { PT(t)^-2 * log(PT(t)) }
fun2 <- function(t) { PT(t)^-2 }
fitahr <- coxphw(Surv(time, event) ~ fun(time) * trt + fun2(time) * trt + x,
                 data = my.data, template = "AHR")
yahr <- predict(fitahr, type = "slice.time", x = "time", z = "trt", newx = plot.x)

fitph <- coxphw(Surv(time, event) ~ fun(time) * trt + fun2(time) * trt + x,
                data = my.data, template = "PH")
yph <- predict(fitph, type = "slice.time", x = "time", z = "trt", newx = plot.x)

plot(yahr, addci = FALSE)
lines(yph$estimates$time, yph$estimates$coef, lty = 2)
legend("bottomright", legend = c("AHR", "PH"), bty = "n", lty = 1:2,
       inset = 0.05)
}

\keyword{survival}
