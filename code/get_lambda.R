
get_lambda<-function(X, Y, lambda_dis){
  
  R <- det(cov(cbind(X,Y)))
  R_network <- det(cov(X))
  R_responese <- det(cov(Y))
  
  t <- R/(R_network*R_responese)
 
  list("p" = ecdf(lambda_dis)(t), "lambda" = t)
  
}
