#' plot goodness of fit for latent variable models
#' 
#' Goodness of fit statistics 
#' 
#' 
#' @usage plot_gof(model)
#' @param GOF a matrix of GOF statistics, true values in the first row
#' @return a plot of goodness of fit statistics 
#' @author Selena Wang
#' @examples
#' 
#' 
#' 
#' @export plot_gof
plot_gof <- function(GOF){
    for(k in 1:ncol(GOF))
    {
      hist(GOF[-1,k],xlim=range(GOF[,k]),main="",prob=TRUE,
           xlab=colnames(GOF)[k],ylab="",yaxt="n")
      abline(v=GOF[1,k],col="red")
    }
  
}

