% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C2PL_iw2.R
\name{C2PL_iw2}
\alias{C2PL_iw2}
\title{IW-GVEM Algorithm for Confirmatory M2PL Analysis}
\usage{
C2PL_iw2(
  data,
  model = matrix(1, ncol(data)),
  criterion = "BIC",
  iter = 200,
  eps = 0.001,
  c = 1,
  S = 10,
  M = 10,
  lr = 0.1,
  SE.level = NULL
)
}
\arguments{
\item{data}{An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})}

\item{model}{A \eqn{J\times K} binary matrix of loading indicators  (all items load on the only dimension by default)}

\item{criterion}{Information criterion for model selection, one of \code{'GIC'} (recommended), \code{'BIC'}, or \code{'AIC'}}

\item{iter}{Maximum number of iterations}

\item{eps}{Termination criterion on numerical accuracy}

\item{c}{Constant for computing GIC}

\item{S}{Sample size for approximating the expected lower bound}

\item{M}{Sample size for approximating a tighter lower bound}

\item{lr}{Learning rate for the Adam optimizer}

\item{SE.level}{Accuracy level of Gaussian quadrature for \code{mvQuad} to compute standard errors (SEs are not computed if \code{SE.level} is \code{NULL})}
}
\value{
An object of class \code{vemirt_DIF}, which is a list containing the following elements:
  \item{N}{Number of respondents}
  \item{niter0}{Number(s) of iterations for initialization}
  \item{fit}{The only element of \code{all}}
  \item{best}{Equal to \code{1}}
  \item{all}{A list of model which has one element:}
  \item{ ...$niter}{Number(s) of iterations}
  \item{ ...$SIGMA}{Person-level posterior covariance matrices}
  \item{ ...$MU}{Person-level posterior mean vectors}
  \item{ ...$Sigma}{Population covariance matrix}
  \item{ ...$Mu}{Population mean vector}
  \item{ ...$a}{Slopes}
  \item{ ...$b}{Intercepts}
  \item{ ...$SE.a}{Standard errors of \code{a}}
  \item{ ...$SE.b}{Standard errors of \code{b}}
  \item{ ...$ll}{Estimated lower bound of log-likelihood}
  \item{ ...$l0}{Number of nonzero elements in \code{model}}
  \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
  \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
  \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
}
\description{
IW-GVEM Algorithm for Confirmatory M2PL Analysis
}
\examples{
\dontrun{
with(C2PL_data, C2PL_iw2(data, model, SE = TRUE))}
}
\seealso{
\code{\link{C2PL_gvem}}, \code{\link{C2PL_iw}}, \code{\link{D2PL_gvem}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
}
\author{
Weicong Lyu <weiconglyu@um.edu.mo>
}
