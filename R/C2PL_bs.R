#' Bootstrap Version of GVEM Confirmatory Analysis for M2PL
#' @description A bootstrap version of GVEM (i.e., GVEM-BS) can be implemented
#' to correct the bias on item parameters and compute standard errors under confirmatory M2PL models
#' @param gvem_result a list that includes exploratory or confirmatory GVEM
#' results for M2PL models.
#' @param boots the number of bootstrap samples; default is 5
#'
#' @return a list containing the following objects:
#'   \item{boots_a}{item discrimination parameters corrected by bootstrap sampling, a \eqn{J \times K} \code{matrix}}
#'   \item{boots_b}{item difficulty parameters corrected by bootstrap sampling, a vector of length \eqn{J}}
#'   \item{sd_a}{stardard errors of item discrimination parameters, a \eqn{J \times K} matrix}
#'   \item{sd_b}{stardard errors of item difficulty parameters, a vector of length \eqn{J}}
#'
#' @author Jiaying Xiao <jxiao6@uw.edu>
#' @seealso \code{\link{C2PL_gvem}}, \code{\link{C2PL_iw}}
#' @export
#'
#' @examples
#' \dontrun{
#' gvem_result <- with(C2PL_data, C2PL_gvem(data, model))
#' C2PL_bs(gvem_result, boots=10)}
C2PL_bs <- function(gvem_result,boots=5){
  domain <- dim(gvem_result$ra)[2]
  indic <- gvem_result$Q_mat
  item <- dim(gvem_result$ra)[1]
  N1 <- dim(gvem_result$mu_i)[2]
  mu <- rep(0,domain)
  #est_gvem <- cbind(gvem_result$ra,gvem_result$rb)
  rs2<-vector("list",boots)
  for (j in 1:boots) {
    set.seed(j*111)
    theta<-MASS::mvrnorm(N1,mu,gvem_result$rsigma)
    u2<-1/(1+exp(-(theta%*%t(gvem_result$ra)-do.call(rbind, replicate(N1, list(t(gvem_result$rb)))))))
    u2<-1*(u2>matrix(runif(N1*item),nrow=N1))
    rs2[[j]]=C2PL_gvem(u2,indic)
  }
  #check the result
  new_a<-sapply(rs2,'[',"ra")
  new_b<-sapply(rs2,'[',"rb")
  boots_a<-apply(simplify2array(new_a), 1:2, mean)
  boots_b<-Reduce("+", new_b) / length(new_b)
  sd_a <-apply(simplify2array(new_a), 1:2, sd)
  sd_b <-apply(simplify2array(new_b), 1, sd)
  #bootstrap correction: 2*gvem - boots
  boots_a<-2*gvem_result$ra - boots_a
  boots_b<-2*gvem_result$rb - boots_b
  new.vemirt_FA(list(boots_a=boots_a,boots_b=boots_b,sd_a=sd_a,sd_b=sd_b))
}

