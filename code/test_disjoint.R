#' test_disjoint
#' 
#' Perform testing of independence based on two sets of normal variables
#' 
#' This command performs testing of independence based on two sets of normal variables based on Anderson (1962).
#' Any linear transformations of the two sets of variables are considered.
#' 
#' 
#' @usage test_disjoint(X, Y, N_sim=100, sig = 0.05)
#' @param X the first set of variables to be tested.
#' @param Y the second set of variables to be tested
#' @param N_sim number of simulation performed to approximate the Wilk's Lambda distribution
#' @param seed random seed
#' @param sig significance level for testing independence
#' @return 
#' \item{p value}{P value associated with the test} \item{lambda}{Wilk's Lambda test statistics}
#' @author Selena Wang
#' @examples
#' 
#' @export test_disjoint
#' 
test_disjoint<-function(X, Y, N_sim=100000, sig = 0.05){
  
  R <- det(cov(cbind(X,Y)))
  R_network <- det(cov(X))
  R_responese <- det(cov(Y))
  
  t <- R/(R_network*R_responese)
  
  K<- ncol(X)
  D<- ncol(Y)
  N<-nrow(X)
  
  beta_values<- matrix(NA, nrow = N_sim, ncol = D)
  
  for (d in 1:D){
    beta_values[,d]=rbeta(N_sim, shape1 = (N-K-D+d)/2, shape2 = K/2)
    
  }
  
  lambda_dis <- apply(beta_values, 1, prod)
  
  # quantile(lambda_dis, probs = .05)
  # 
  # sum(lambda_dis<t)/N_sim
  # hist(lambda_dis)
  # ecdf(lambda_dis)(quantile(lambda_dis, probs = sig))
  # 
  list("p" = ecdf(lambda_dis)(t), "lambda" = t)
  
}

tmp<-sqrt(solve(t(X) %*% X)) %*% (t(X) %*% Y) %*% solve(t(Y) %*% Y) %*% (t(Y) %*% X) %*% sqrt(solve(t(X) %*% X))
sqrt(eigen(tmp))
