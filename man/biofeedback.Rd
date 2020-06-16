\name{biofeedback}

\alias{biofeedback}

\docType{data}

\title{  Biofeedback Treatment Data }

\description{ In this study the effect of biofeedback treatment on time until treatment \code{success} 
was evaluated in patients suffering from aspiration after head and neck 
surgery. The outcome of interest was the time from start of treatment until the patient achieved full
swallowing rehabilitation (\code{thdur}). Patients were randomized into two groups (\code{bfb}): 
one group of patients received videoendoscopic biofeedback treatment; the other group received the 
conservative treatment including thermal stimulation with ice and exercises for the lips, tongue, 
laryngeal closure and elevation. Treatment was started as soon as the healing process after surgery 
was finished (\code{thealing}).
}

\usage{data(biofeedback)}

\format{
  A data frame with 33 observations on the following 6 variables:
  \describe{
    \item{\code{id}}{ the patient id. }
    \item{\code{success}}{ of treatment within the first 100 days; either 0 = no success or 1 = success. }
    \item{\code{thdur}}{ the duration of therapy in days. }
    \item{\code{bfb}}{ indicates the treatment group; either 0 = conservative or 1 = biofeedback. }
    \item{\code{theal}}{ time from surgery to start of therapy in days. }
    \item{\code{log2heal}}{ log2-transformed time from surgery to start of therapy. }
  }
}

%%\details{ If necessary, more details than the __description__ above ~~ }

\source{  Data were supplied by Dr. D.-M. Denk, who gave permission to freely distribute 
  the data and use them for non-commercial purposes. }

\references{ Denk, D.-M. & Kaider, A. (1997). Videoendoscopic Biofeedback: A Simple Method to Improve 
the Efficacy of Swallowing Rehabilitation of Patients After Head and Neck Surgery. \emph{ORL J OTO-RHINO-LARY}
\bold{59}, 100-105. }

\examples{
data("biofeedback")

plot(survfit(Surv(thdur, success) ~ bfb, data = biofeedback), lty = 1:2, las = 1,
     xlab = "time (days)", ylab = "propability of success")
     
coxphw(Surv(thdur, success) ~ bfb, data = biofeedback,  template = "AHR")     
}

\keyword{datasets}
