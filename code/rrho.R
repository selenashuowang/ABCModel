#' Metropolis update for dyadic correlation
#' 
#' Metropolis update for dyadic correlation
#' 
#' 
#' @usage rrho(Z, rho, s2 = 1,offset=0, asp=NULL)
#' @param Z n X n normal relational matrix
#' @param rho current value of rho
#' @param s2 current value of s2
#' @param offset matrix of the same dimension as Z. It is assumed that 
#' Z-offset is equal to dyadic noise, so the offset should contain any 
#' additive and multiplicative effects (such as 
#' \code{ U\%*\%t(U) +  outer(a,b,"+")  }   ) 
#' @param asp use arc sine prior (TRUE) or uniform prior (FALSE) 
#' 
#' @return a new value of rho
#' @author Selena Wang
#' @export rrho
rrho <-
function(Z,rho,s2=1,offset=0,asp=NULL)
{
  if(is.null(asp)){ asp<-TRUE } 

  E<-Z-offset
  EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] )/sqrt(s2)
  emcp<-sum(EM[,1]*EM[,2])
  emss<-sum(EM^2)

  m<- nrow(EM)
  sr<- 2*(1-cor(EM)[1,2]^2)/sqrt(m)
  
  if(sr==0){
    sr=.001
  }

  rho1<-rho+sr*qnorm( runif(1,pnorm( (-1-rho)/sr),pnorm( (1-rho)/sr)))

  lhr<-(-.5*(m*log(1-rho1^2)+(emss-2*rho1*emcp)/(1-rho1^2)))-
       (-.5*(m*log(1-rho^2 )+(emss-2*rho*emcp )/(1-rho^2 )))   +
       asp*( (-.5*log(1-rho1^2)) - (-.5*log(1-rho^2)) )

  if(log(runif(1))<lhr) { rho<-rho1 }
  min(abs(rho),.995)*sign(rho)
}
