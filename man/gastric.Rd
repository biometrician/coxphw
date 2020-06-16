\name{gastric}

\alias{gastric}

\docType{data}

\title{ Gastric Cancer Data }
\description{ A data set of survival times of patients with locally advanced, nonresectable gastric
              carcinoma. The patients were either treated with chemotherapy plus radiation or
              chemotherapy alone. }

\usage{data(gastric)}

\format{
  A data frame with 90 observations on the following 4 variables:
  \describe{
    \item{\code{id}}{ unique patient id. }
    \item{\code{radiation}}{ treatment of either 0 = chemotherapy alone or 1 = chemotherapy plus radiation. }
    \item{\code{time}}{ survival time in days. }
    \item{\code{status}}{ 0 = censored or 1 = death. }
  }
}

%%\details{  ~~ If necessary, more details than the __description__ above ~~ }

\source{ Stablein DM, Carter J, Novak JW. (1981) Analysis of Survival Data with Nonproportional Hazard
Functions. \emph{Controlled Clinical Trials} \bold{2}:149-159. OR

\url{http://www.mayo.edu/research/documents/gastrichtml/DOC-10027680} }

\references{  Gastrointestinal Tumor Study Group. (1982) A Comparison of Combination
Chemotherapy and Combined Modality Therapy for Locally Advanced Gastric Carcinoma. \emph{Cancer}
\bold{49}:1771-7. }

\examples{
data("gastric")
plot(survfit(Surv(time, status) ~ radiation, data = gastric), lty = 1:2, las = 1,
     xscale = 365.25, xlab = "time (years)", ylab = "survival distribution function")

coxphw(Surv(time, status) ~ radiation, data = gastric, template = "AHR")
}

\keyword{datasets}
