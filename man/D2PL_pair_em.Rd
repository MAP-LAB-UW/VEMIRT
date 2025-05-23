% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/D2PL_pair_em.R
\name{D2PL_pair_em}
\alias{D2PL_pair_em}
\title{EM Algorithm with ADMM for DIF Detection Using Group Pairwise Truncated \eqn{L_1} Penalty in 2PL Models}
\usage{
D2PL_pair_em(
  data,
  group = rep(1, nrow(data)),
  Lambda0 = if (length(unique(group)) == 1) 0 else seq(0.5, 1.5, by = 0.1),
  Tau = if (length(unique(group)) == 1) 0 else c(Inf, seq(0.05, 0.3, by = 0.05)),
  rho0 = 0.5,
  level = 10,
  criterion = "BIC",
  iter = 200,
  eps = 0.001,
  c = 1,
  verbose = TRUE
)
}
\arguments{
\item{data}{An \eqn{N\times J} binary matrix of item responses (missing responses should be coded as \code{NA})}

\item{group}{An \eqn{N} dimensional vector of group indicators from \code{1} to \code{G} (all respondents are in the same group by default)}

\item{Lambda0}{A vector of \code{lambda0} values for truncated \eqn{L_1} penalty (\code{lambda} equals \code{sqrt(N) / G * lambda0})}

\item{Tau}{A vector of \code{tau} values for truncated \eqn{L_1} penalty (becomes \eqn{L_1} penalty when \code{tau} equals \code{Inf})}

\item{rho0}{A value of \code{rho} for augmented Lagrangian in ADMM (\code{tau} equals \code{sqrt(N) / G * tau0})}

\item{level}{Accuracy level of Gaussian quadrature for \code{mvQuad}}

\item{criterion}{Information criterion for model selection, one of \code{'BIC'} (recommended), \code{'AIC'}, or \code{'GIC'}}

\item{iter}{Maximum number of iterations}

\item{eps}{Termination criterion on numerical accuracy}

\item{c}{Constant for computing GIC}

\item{verbose}{Whether to show the progress}
}
\value{
An object of class \code{vemirt_DIF}, which is a list containing the following elements:
  \item{N}{Number of respondents}
  \item{niter0}{Number(s) of iterations for initialization}
  \item{fit}{The best (with lowest information criterion) model, which is an element of \code{all}}
  \item{best}{The index of \code{fit} in \code{all}}
  \item{all}{A list of models which has the same length as \code{Lambda0}:}
  \item{ ...$lambda0}{Corresponding element in \code{Lambda0}}
  \item{ ...$lambda}{\code{sqrt(N) / G * lambda0}}
  \item{ ...$tau}{Corresponding element in \code{Tau}}
  \item{ ...$rho0}{Same as \code{rho0} in input}
  \item{ ...$rho}{\code{sqrt(N) / G * rho0}}
  \item{ ...$niter}{Number(s) of iterations}
  \item{ ...$Sigma}{Group-level covariance matrices}
  \item{ ...$Mu}{Group-level mean vectors}
  \item{ ...$a}{Slopes}
  \item{ ...$b}{Intercepts}
  \item{ ...$d.a}{Group pairwise differences of slopes}
  \item{ ...$d.b}{Group pairwise differences of intercepts}
  \item{ ...$u.a}{Lagrangian multipliers of corresponding elements in \code{d.a}}
  \item{ ...$u.b}{Lagrangian multipliers of corresponding elements in \code{d.b}}
  \item{ ...$ll}{Log-likelihood}
  \item{ ...$l0}{Number of nonzero D2PL parameters in \code{gamma} and \code{beta}}
  \item{ ...$AIC}{Akaike Information Criterion: \code{-2*ll+l0*2}}
  \item{ ...$BIC}{Bayesian Information Criterion: \code{-2*ll+l0*log(N)}}
  \item{ ...$GIC}{Generalized Information Criterion: \code{-2*ll+c*l0*log(N)*log(log(N))}}
}
\description{
EM Algorithm with ADMM for DIF Detection Using Group Pairwise Truncated \eqn{L_1} Penalty in 2PL Models
}
\examples{
\dontrun{
with(D2PL_data, D2PL_pair_em(data, group, Tau = c(Inf, seq(0.01, 0.05, by = 0.01))))}
}
\seealso{
\code{\link{D2PL_em}}, \code{\link{D2PL_gvem}}, \code{\link{D2PL_lrt}}, \code{\link{coef.vemirt_DIF}}, \code{\link{print.vemirt_DIF}}, \code{\link{summary.vemirt_DIF}}
}
\author{
Weicong Lyu <weiconglyu@um.edu.mo>
}
