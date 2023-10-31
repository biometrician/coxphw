\name{concord}

\alias{concord}

\title{ Compute Generalized Concordance Probabilities for Objects of Class \code{coxphw} or  
        \code{coxph} }

\description{ Compute generalized concordance probabilities with accompanying 
              confidence intervals for objects of class \code{coxphw} or \code{coxph}. }

\usage{ concord(fit, digits = 4) }
        
\arguments{
  \item{fit}{ an object of class \code{coxphw}. }
  \item{digits}{ integer indicating the number of decimal places to be used. Default is 4. }
}
    
\details{ The generalized concordance probability is defined as \eqn{P(T_i < T_j | x_i = x_j + 1)}
          with \eqn{T_i} and \eqn{T_j} as survival times of randomly chosen subjects with covariate 
          values \eqn{x_i} and \eqn{x_j}, respectively. Assuming that \eqn{x_i} and \eqn{x_j} are
          1 and 0, respectively, this definition includes a two-group comparison.
          
          If proportional hazards can be assumed, the generalized concordance probability can also
          be derived from Cox proportional hazards regression (\code{coxphw} with \code{template = "PH"} 
          or \code{coxph}) or weighted Cox regression as suggested by Xu and O'Quigley (2000)
          (\code{coxphw} with \code{template = "ARE"}).

          If in a fit to \code{coxphw} the \code{betafix} argument was used, then for the 
          fixed parameters only the point estimates are given. }    
          
\value{ A matrix with estimates of the generalized concordance probability with 
        accompanying confidence intervalls for each explanatory variable in the model. }

\author{Daniela Dunkler}

\references{ Dunkler D, Schemper M, Heinze G. (2010) Gene Selection in Microarray Survival Studies Under 
Possibly Non-Proportional Hazards. \emph{Bioinformatics} \bold{26}:784-90. 

Xu R and O'Quigley J (2000). Estimating Average Regression Effect Under Non-Proportional Hazards.
\emph{Biostatistics} \bold{1}, 423-439. }

\seealso{ \code{\link{coxphw}} }

\examples{
data("gastric")
fit <- coxphw(Surv(time, status) ~ radiation, data = gastric, template = "AHR")
concord(fit)
}

\keyword{survival}
