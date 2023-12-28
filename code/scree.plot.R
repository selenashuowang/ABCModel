#' @title Identify the number of dimensions in the model using scree plots
#'
#' @description identify the number of dimensions in the model using scree plots
#'
#' @param M Matrix: containing the posterior estimates that are used for identifying dimensions
#' @param D number: maximum number of dimensions
#' @param plot logical: plot results?
#' @param eigen logical: plot eigen values? default is plotting the proportion of variation
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{results} eigen values and proportions of variation
#'  }
#'
#' @export
#' @author Selena Wang
#'
#' @examples





scree.plot<-function(M, D, eigen = FALSE, plot=FALSE ){
  
  f.norm=sqrt(sum(M^2))^2

  
  prop <-matrix(NA, nrow = D, ncol=3)
  colnames(prop) = c("Eigen value","Proportion of variation","Dimension")
  rownames(prop) = seq(1,D)
  prop[,"Dimension"] =seq(1,D)
  prop[,"Eigen value"] = eigen(M %*% t(M), symmetric = TRUE)$values[1:D]
  prop[,"Proportion of variation"] = eigen(M %*% t(M), symmetric = TRUE)$values[1:D]/f.norm
  
  if(plot & !eigen){
    barplot(prop[,"Proportion of variation"], main=NULL, ylab = "Proportion of variation",
            xlab="Dimension", 
            beside=TRUE, ylim = c(0,1))
  }
  
  if(plot & eigen){
    
    plot(prop[,"Dimension"], prop[,"Eigen value"], ylab = "Eigen value",
         xlab="Dimension",)
    
    lines(prop[,"Dimension"][order(prop[,"Dimension"])], prop[,"Eigen value"][order(prop[,"Dimension"])], pch=16)

  }
  
  prop
}