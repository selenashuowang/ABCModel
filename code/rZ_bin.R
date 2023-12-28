#' Simulate Z based on a probit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix X
#' 
#' 
#' @usage rZ_bin(Z, EZ, rho, X)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param X square binary relational matrix
#' @param s2 network variance
#' @return a square matrix, the new value of Z
#' @author Selena Wang
#' @export rZ_bin
rZ_bin <-
function(Z,EZ,rho,X, s2=1)
{ 
  # simulates Z under the contraints
  # (1)  X[i,j]=1   => Z[i,j]>0
  # (2)  X[i,j]=0   => Z[i,j]<0
  
  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  X[is.na(X)]<- -1
  for(x in c((-1):1))
  { 
    lb<-c(-Inf,-Inf,0)[x+2] ; ub<-c(Inf,0,Inf)[x+2]
    
    for(tri in 1:2)
    { 
      if(tri==1){ up<-ut & X==x }
      if(tri==2){ up<-lt & X==x }
      
      ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
      zup<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
      zerr<-which(abs(zup)==Inf)
      if(length(zerr)>0){ zup[zerr]<-(Z[up])[zerr] }
      Z[up]<-zup
    }
  }
  
  ##
  c<-(sqrt(1+rho) + sqrt(1-rho))/2
  d<-(sqrt(1+rho) - sqrt(1-rho))/2
  E<-matrix(rnorm(nrow(X)^2, sd= sqrt(s2)),nrow(X),nrow(X))
  ZP<-EZ + c*E + d*t(E) 
  A<-( (X== -1) | ( sign(ZP) == sign(X-.5)) ) ; diag(A)<-TRUE
  A<-A & t(A)
  Z[A]<-ZP[A]
  #Z[!A]=0
  ##

  ## this line now redundant because of previous chunk
  diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))
  ##

  Z
}

