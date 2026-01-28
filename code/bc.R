#' Attribute informed brain connectivity
#' 
#' An MCMC algorithm for fitting the joint network and item response theory (JNIRT) model 
#' to network data as well as item responses
#' 
#' This command provides posterior inference for parameters in JNIRT models 
#' assuming normal or binary data types
#' 
#' "nrm": assuming the data is normally distributed, can be network or item responses
#' 
#' "bin": assuming the data is binary, can be network or item responses
#' 
#' @usage jnirt(X, Y, family_network, family_responses, K=0, rvar = FALSE ,
#' cvar = FALSE,  dcor = FALSE, model = UU, 
#' intercept=TRUE, seed = 1, nscan =
#' 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE,
#' prior=list())
#' @param X a list of V x V brain connectivity data. 
#' @param W a matrix of N x Q covariates for the connectivity data.
#' @param K integer: dimension of the multiplicative effects (can be zero) in the connectivity data.
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param prior list: A list of hyperparameters for the prior distribution
#' @return 
#' \item{VC}{posterior samples of the variance parameters}
#' \item{COV}{posterior samples of the covariance parameters}
#' \item{BETAPM}{posterior mean of the regression coefficient parameters for the connectivity data}
#' \item{APM}{posterior mean of connectivity intercepts} \item{U}{posterior estimates of multiplicative
#' row effects u} \item{U_1}{the last iteration of the multiplicative
#' row effects u} 
#' \item{UVPM}{posterior mean of UV} 
#' \item{XPM}{mean of the generated data based on posterior samples} 
#' \item{X}{observed X} 
#' \item{UVC}{posterior samples of elements of the connectivity UU} 
#' \item{EFlPM}{posterior mean estimates of X}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' \item{input}{input values} 
#' @author Selena Wang
#' @examples
#' 
#' @export bc
#' 





bc<- function(X, W,  K=2,
               indices = NULL,
               seed = 1, nscan = 10000, burn = 500, odens = 25,
               print = TRUE, gof=TRUE, plot=TRUE, 
               prior=list())
{ 
  ## record model set up 
  input<-list(K=K, nscan=nscan, burn=burn, odens=odens, prior=prior, indices = indices)
  # set random seed
  set.seed(seed)
  
  # starting Fl values
  Fl<-X

  N<-length(X)
  V<-nrow(X[[1]])
  
  if(is.null(prior$Sutheta0)){ prior$Sutheta0<-diag(K) } 
  if(is.null(prior$etautheta)){ prior$etautheta<-(K+2) } 
  
  # starting intercept values
  a<-matrix(sapply(Fl, mean, na.rm=TRUE ))
  
  # starting beta values
  if(!is.null(W)){ beta<-matrix(rep(0,ncol(W)), nrow = ncol(W),ncol=1) }else{beta<-NULL} 
  
  s2<-mean(sapply(1:length(X), function(x) mean((Fl[[x]] - a[x,])^2, na.rm=TRUE)), na.rm=TRUE)
  
  # U
  tmp<-sapply(1:length(X), function(x) Fl[[x]] - a[x,], simplify = FALSE)
  for(i in 1:length(X)){tmp[[i]][is.na(tmp[[i]])]=0}
  E<-Reduce('+', tmp)/N
  
  U<-matrix(0,V,K) 
  if(K>0) 
  {  
    sE<-svd(E)
    U<-sE$u[,1:K,drop=FALSE]%*%diag(sqrt(sE$d[1:K]),nrow=K)
  }
  
  # Initialize output matrices
  if(!is.null(W)){  
    BETA<- matrix(0,nrow = 0, ncol = ncol(W))
  }else{
    BETA<- matrix(0,nrow = 0, ncol = 0)
  }
  
  BETAPS <- beta * 0 
  XPM<-EFlPM<-list()
  UVPS <- U %*% t(U) * 0 
  APS<- rep(0,length(X))  
  names(APS)<- names(X)
  rownames(U)<-rownames(X[[1]])
  
  # Initialize UVC exactly like ABC model
  UVC<-matrix(nrow=0,ncol=(V*(V-1)/2)) 
  
  VC<-matrix(nrow=0,ncol=1+length(c(seq(1,(K)*(K+1)/2)))) 
  colnames(VC) <- c(paste("Su",seq(1,(K)*(K+1)/2),sep=""),
                   "ve_connectivity") 
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  for (s in 1:(nscan + burn)) 
  { 
    if(!is.null(W)){
      EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + as.numeric(W[x,] %*% beta) + U %*% t(U)), simplify = FALSE)
    } else {
      EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + U %*% t(U)), simplify = FALSE)
    }
    
    # update s2
    s2<-rs2(Fl,offset = EFl, nu2=prior$nu2,s20=prior$s20)  
    
    # update beta, a 
    tmp <- rbeta_a_fc_per(Fl,W=W,s2=s2,offset=U%*%t(U),ivA=prior$ivA,beta0=prior$beta0,S0=prior$S0) 
    beta <- tmp$beta 
    a <- tmp$a
    
    ## update variances
    tmp <- rSu(U, Su0=prior$Sutheta0,etau=prior$etautheta) 
    Su=tmp$Su
    
    # update U,V
    if (K > 0)
    {
      U<-rU_con(Fl,U, Su, s2, offset=sapply(EFl, function(x) x-U%*%t(U), simplify = FALSE))
      U <- U - colMeans(U)
      
      if(s == 1){U_target<-U}
      if(s > 1){
        tmp <- Procrustes(U, U_target,
                          translate = FALSE,
                          dilate = FALSE,
                          sumsq = FALSE)
        U<-tmp$X.new
      }
    }
    
    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  
      # save results
      BETA<-rbind(BETA, as.vector(beta)) 
      VC<-rbind(VC, c(Su[upper.tri(Su, diag = T)], s2)) 
      BETAPS<-BETAPS+beta
      
      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(U)
      
      # Store UVC exactly like ABC model
      tmp<-U %*% t(U)
      tmp.1<-NULL
      tmp.1 <- c(tmp.1, tmp[upper.tri(tmp,diag = FALSE)])
      UVC<-rbind(UVC, c(tmp.1))
      
      APS <- APS + a
    }
  } # end MCMC   
  
  # posterior means 
  BETAPM <- BETAPS/nrow(VC)
  APM <- APS/nrow(VC)
  UVPM <- UVPS/nrow(VC)
  
  if(!is.null(W)){
    EFlPM<-sapply(1:length(X), function(x) (APM[x,] + as.numeric(W[x,] %*% BETAPM) + UVPM), simplify = FALSE)
  } else {
    EFlPM<-sapply(1:length(X), function(x) (APM[x,] + UVPM), simplify = FALSE)
  }
  
  names(APM)<-names(X)
  rownames(UVPM)<-colnames(UVPM)<-rownames(X[[1]])
  
  # asymmetric output 
  UDV<-eigen(UVPM)
  U_1<-UDV$vectors[,seq(1,K,length=K)]%*%diag(sqrt(UDV$values[seq(1,K,length=K)]),nrow=K)
  
  rownames(U)<-rownames(X[[1]]) 
  
  fit <- list(BETAPM=BETAPM, VC=VC, APM=APM, U=U_1, UVPM=UVPM, EFlPM=EFlPM,
              X=X, UVC=UVC, input=input)
  
  class(fit) <- "bc"
  fit
}