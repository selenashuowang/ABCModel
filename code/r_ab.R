#' Conditional simulation of additive effects
#' 
#' Simulates from the joint full conditional distribution of (a,b)
#' 
#' @param Z n X n normal relational matrix
#' @param Sab row and column covariance
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param offset a matrix of the same dimension as Z. It is assumed that 
#' Z-offset follows a SRRM, so the offset should contain any multiplicative 
#' effects (such as \code{U\%*\% t(U) } )
#' @return  \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Selena Wang
#' @export r_ab_fc
r_ab<-function(Z,Sab,rho,s2=1,offset=0)
{
  Z<-Z-offset
  n<-nrow(Z)
  ###

  ### decorrelation
  Se<-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2<-mhalf(solve(Se))
  td<-iSe2[1,1] ; to<-iSe2[1,2]
  Sabs<-iSe2%*%Sab%*%iSe2
  tmp<-eigen(Sabs)
  k<-sum(zapsmall(tmp$val)>0 )

  Zs<-td*Z+to*t(Z)
  zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) 
  ###


  ### row and column reduction
  ab<-matrix(0,nrow(Z),2)
  if(k>0)
  {
    G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
    K<-matrix(c(0,1,1,0),2,2)
    A<-n*t(G)%*%G + diag(k)
    B<-t(G)%*%K%*%G
    iA0<-solve(A)
    C0<- -solve(A+ n*B)%*%B%*%iA0

    iA<-G%*%iA0%*%t(G)
    C<-G%*%C0%*%t(G)

    }
  
  ###


  ###
 
  ### simulate a, b 
  if(k>0) 
  {
    E<- Zs
    er<-rowSums(E) ; ec<-colSums(E) ; es<-sum(ec) 
    m<-t(t(crossprod(rbind(er,ec),t(iA0%*%t(G)))) + rowSums(es*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab<- w%*%t(G)%*%solve(iSe2) 
  }

list(a=ab[,1],b=ab[,2] )  
}
