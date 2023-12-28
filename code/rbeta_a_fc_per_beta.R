#' Conditional simulation of intercept and regression coefficients
#' 
#' Simulates from the joint full conditional distribution of (beta,intercept)
#' in a brain connectivity model
#' 
#' @param Fl a list of V X V normal connectivity matrix
#' @param W N x Q covariate matrix
#' @param s2 dyadic variance
#' @param offset a V by V matrix, UU 
#' @param S0 prior precision matrix for regression parameters
#' @param beta0 prior mean vector for regression parameters 
#' @param ivA prior inverse variance for the intercept parameters 
#' 
#' @return \item{beta}{regression coefficients} \item{a}{subject-specific intercept}
#' @author Selena Wang
#' @export rbeta_a_fc
rbeta_a_fc_per_beta<-
function(Fl,W=NULL,s2=1,offset=offset,ivA=NULL,beta0=NULL,S0=NULL)
{
  V<-nrow(Fl[[1]]) ; N <-length(Fl) 

  if(is.null(ivA))
  { 
    ivA <- matrix(1, nrow = 1, ncol = 1)
    
  } 
  
  
  ### set priors 
  if(!is.null(W) & is.null(S0))
  { 
    # g-prior plus small ridge in case XX is singular
    S0 <- diag(1,ncol(W))
  } 
  
  if(!is.null(W) & is.null(beta0))
  { 
    beta0 <- matrix(rep(0,ncol(W)), nrow = ncol(W),ncol=1)
    
  } 
  
  one_vector <- matrix(1,nrow = V*(V-1)/2, ncol = 1)
  #G <- t(kronecker(diag(1, nrow = N), t(one_vector)))
  
  to<-as.numeric(sqrt(solve(s2)))
  a=matrix(NA, nrow = N, ncol = 1)
  
  
  iQA = solve( to^2*t(one_vector) %*% one_vector + ivA) 
  
  for (i in 1:length(Fl)){
    R=Fl[[i]]-offset
    r=matrix(R[upper.tri(R, diag = FALSE)])
    r[is.na(r)]=0
    rs<- to*r
    

      m<-rs 
      
      lA <- to*t(one_vector)%*%m
      a[i,] <- iQA%*%lA + t(chol(iQA))%*%rnorm(1) 

    

    
  }
  
  
  # 
  # r=matrix(NA, nrow=0, ncol=1)
  # for (i in 1:length(Fl)){
  #   R=Fl[[i]]-offset
  #   r=rbind(r,matrix(R[upper.tri(R, diag = FALSE)]))
  #   
  # 
  # }
  # 
  # r[is.na(r)]=0
  # rs<- to*r
  # 
  # iQA = solve( to^2*t(G) %*% G + ivA) 
  # 
  # 
  # 
  # ### updata beta check this
  # if(!is.null(W)){
  #   S <- to^2*G %*% iQA %*% t(G)
  #   
  #   Wg <- kronecker(W,one_vector)
  #   iQB <- solve(to^2*t(Wg) %*% S %*% Wg + S0)
  #   lB <- to^2*t(Wg) %*% S %*% r + S0 %*% beta0
  #   beta <- iQB%*%lB + t(chol(iQB))%*%rnorm(ncol(W))
  #   
  #   m<- rs -to* Wg%*%beta
  #   lA <- to*t(G)%*%m
  #   a <- iQA%*%lA + t(chol(iQA))%*%rnorm(N) 
  #   
  #   return(list(beta=beta,a=a ) )
  # }else{
  #   m<-rs 
  #   lA <- to*t(G)%*%m
  #   a <- iQA%*%lA + t(chol(iQA))%*%rnorm(N) 
  #   return(list(beta=NULL,a=a ) )
  # }

  return(list(beta=NULL,a=a ))
 
}
