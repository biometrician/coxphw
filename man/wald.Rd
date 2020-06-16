\name{wald}

\alias{wald}

\title{ Wald-Test for Model Coefficients }

\description{ Obtain Wald chi-squared test statistics and p-values for one or more regression
              coefficients given their variance-covariance matrix. }

\usage{ wald(coeff, cov, index = NULL, h0 = NULL) }

\arguments{
  \item{coeff}{ a vector with regression coefficients. }
  \item{cov}{ a variance-covariance matrix. }
  \item{index}{ an integer specifiying which parameters should be jointly tested. Default is to test
                all parameters given in \code{coeff} and \code{cov}.}
  \item{h0}{ a vector with the same length as \code{coeff} stating the null hypothesis for the test.
             Default is 0 for all coefficients. }
}

\details{ The test is based on the assumption that the coefficients asymptotically follow
          a (multivariate) normal distribution with mean \code{coeff} and a variance-covariance
          matrix \code{cov}. }

\value{ A vector with the follwing components:
 \item{chi2}{ the Wald-test statistic. }
 \item{df}{ degress of freedom. }
 \item{p}{ p-value. }
}

\author{Daniela Dunkler}

%\references{  }

\seealso{ \code{\link{coxphw}} }

%\examples{ }

\keyword{survival}
