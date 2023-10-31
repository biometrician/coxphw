\name{NEWS}

\title{NEWS file for the \code{coxphw} package}

% \section{UPCOMING CHANGES VERSION 4.0.3 }{
% \itemize{
%  \item { predict: add newdata, na.action=na.pas }
%  \item { offset: with numbers  }
%  \item { clean rebuild using roxygen, testthat, etc.  }
 %}
%}

\section{CHANGES \code{coxphw} VERSION 4.0.3 (2023-10-31, compiled under R 4.3.1)}{
   \itemize{
      \item { changes in coxohw.f90, removed Visual Fortran comments }
      \item { delete licence file,  spell checking, add link to git, etc. }
 }
}


\section{CHANGES \code{coxphw} VERSION 4.0.2 (2020-06-16, compiled under R 4.0.1)}{
   \itemize{
      \item { add testthat }
      \item { add information on the sample size in \code{print} and \code{summary} }
      \item { only observations with missing values in variables mentioned in the formula are deleted. }
      \item { add option \code{clusterid} in \code{coxphw} }
 }
}

 \section{CHANGES \code{coxphw} VERSION 4.0.1 (2018-04-16, compiled under R 3.4.4)}{
 \itemize{
  \item { citation and JSS example code is updated due to JSS publication 2018 }
 }
}

\section{CHANGES in \code{coxphw} VERSION 4.0.0 (2015-12-30, compiled under R 3.2.2) }{
 \itemize{
 \item { plotshape function replaced by \code{predict.coxphw} and \code{plot.coxphw.predict}, \code{print.coxhpw.predict} added. }
 \item { \code{coxphw.control} with arguments iter.max, maxhs, xconv, gconv, maxstep, round.times.to, add.constant, pc, pc.time, normalize, and fp.iter is added. }
 \item { \code{coxphw} renamed argument \code{print} into \code{verbose}. \code{maxit} is now renamed as \code{iter.max} as in the survival package. Usage of * in the formula of \code{coxphw} is now consistent to other \pkg{R} packages like \pkg{lm}. Variable names with dots and underlines also work now. offset argument was moved to formula. }
 \item { \code{plotw} function is now renamed to \code{plot.coxphw}. }
 \item { \code{linear.predictors} in \code{coxphw} were corrected if \code{coxphw.control(pc = TRUE)}. }
 \item { Issue with factors in \code{coxphw} is solved. }
 \item { Issue in \code{coxph} in \pkg{survival} (simultaneous use of \code{tt} and \code{cluster}) is solved. }
 }
}

\section{CHANGES in coxphw VERSION 3.0.0 (2014-03-15, compiled under R 3.0.3)}{
 \itemize{
  \item { All fp-transformation functions are declared as internal functions. }
  \item { Add concordance probability to summary. }
  \item { Add functions: coef, concord, plotshape, wald. }
  \item { Include data: gastric and biofeedback. }
  \item { Appearance of plotw plots is slightly changed. }
  \item { Input is simplified: instead of the options AHR, AHR.norobust, ARE, and PH, templates are used now. }
 }
}

