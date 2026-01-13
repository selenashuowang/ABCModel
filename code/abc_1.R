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
#' @usage faster version of abc model
#' @param X a list of V x V brain connectivity data. 
#' @param Y a list of V x P attribute data.
#' @param W a matrix of N x Q covariates for the connectivity data.
#' @param H a matrix of N x Q1 covariates for the attribute data.
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
#' \item{GAMMAPM}{posterior mean of the regression coefficient parameters for the attribute data}
#' \item{THETAPM}{posterior samples of the latent person variable}
#' \item{APM}{posterior mean of connectivity intercepts} \item{BPM}{posterior
#' mean of attribute intercepts} \item{U}{posterior estimates of multiplicative
#' row effects u} \item{U_1}{the last iteration of the multiplicative
#' row effects u} 
#' \item{UVPM}{posterior mean of UV} 
#' \item{Theta_1}{the last iteration of the Theta estimate} 
#' \item{X}{observed X} 
#' \item{Y}{observed Y} 
#' \item{UVC}{posterior samples of elements of the connectivity UU} 
#' \item{TAC}{posterior samples of elements of the Theta} 
#' \item{EFlPM}{posterior mean estimates of X}
#'  \item{ETlPM}{posterior mean estiamtes of Y}
#' values of four goodness-of-fit statistics}
#' \item{input}{input values} 
#' @author Selena Wang
#' @examples
#' 
#' @export abc_1
#' 




abc_1<- function(X, Y,W, H, K=2,
               indices = NULL, indices_irt = NULL,
               seed = 1, nscan = 10000, burn = 500, odens = 25,
               print = TRUE, gof=TRUE, plot=TRUE, 
               prior=list())
{ 
  ## record model set up 
  
  input<-list(K=K, nscan=nscan, burn=burn, odens=odens, prior=prior, indices = indices, indices_irt = indices_irt)
  # set random seed
  set.seed(seed)
  
  
  
  
  
  # starting Fl values
  Fl<-X
  Tl<-Y
  
  N<-length(X)
  V<-nrow(X[[1]])
  P<-ncol(Y[[1]])
  
  
  
  if(is.null(prior$Sutheta0)){ prior$Sutheta0<-diag(P+K) } 
  if(is.null(prior$etautheta)){ prior$etautheta<-(P+K+2) } 
  
  # starting intercept values
  a<-matrix(sapply(Fl, mean, na.rm=TRUE ))
  b<-matrix(sapply(Tl, mean, na.rm=TRUE))
  
  
  # starting beta values
  if(!is.null(W)){ beta<-matrix(rep(0,ncol(W)), nrow = ncol(W),ncol=1) }else{beta<-NULL} 
  if(!is.null(H)){ gamma<-matrix(rep(0,ncol(H)), nrow = ncol(H),ncol=1) } else{gamma<-NULL}
  
  
  
  
  s1<-mean(sapply(1:length(X), function(x) mean((Tl[[x]] - b[x,])^2, na.rm=TRUE)), na.rm=TRUE)
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
  
  # Theta
  
  tmp<-sapply(1:length(Y), function(x) (Tl[[x]] - b[x,]), simplify = FALSE)
  for(i in 1:length(X)){tmp[[i]][is.na(tmp[[i]])]=0}
  
  E<-Reduce('+', tmp)/N
  
  Theta<-matrix(0,V,P) 
  if(P>0) 
  {  
    sETl<-svd(E)
    Theta<-sETl$u[,1:P,drop=FALSE]%*%diag(sqrt(sETl$d[1:P]),nrow=P)
  }
  
  
  
  
  
  # output items
  
  if(!is.null(W)){  BETA<- matrix(0,nrow = 0, ncol = ncol(W))}else{
    BETA<- matrix(0,nrow = 0, ncol = 0)
  }
  if(!is.null(H)){  GAMMA<- matrix(0,nrow = 0, ncol = ncol(H))}else{
    GAMMA<- matrix(0,nrow = 0, ncol = 0)
  }
  BETAPS <- beta * 0 
  GAMMAPS <- gamma * 0 
  THETAPS <- Theta * 0
  XPM<-YPM<-EFlPM<-ETlPM<-list()
  UVPS <- U %*% t(U) * 0 
  APS<-BPS<- rep(0,length(X))  
  names(APS)<-names(BPS)<- names(X)
  rownames(U)<-rownames(Theta)<-rownames(X[[1]])
  #XPS<-YPS<-list()
  
  #GOF<-list()
  #for(i in 1:length(X)){
    #gofXY<-c(gofstats_c(X[[i]]), gofstats_a(Y[[i]]))
    #GOF[[i]]<-matrix(gofXY,1,length(gofXY))  
    #rownames(GOF[[i]])<-"obs"
    #colnames(GOF[[i]])<-names(gofXY)
    
    #YPS[[i]]<-matrix(0,nrow=V,ncol=P) ; dimnames(YPS[[i]])<-dimnames(Y[[i]]) 
    #XPS[[i]]<-matrix(0,nrow=V,ncol=V) ; dimnames(XPS[[i]])<-dimnames(X[[i]]) 
    
  #}
  
  
  
  if(is.null(indices)){
    indices<-matrix(sample(1:nrow(X[[1]]),min(round(nrow(X[[1]])/5),5)*2, replace = FALSE), nrow=2)
    
  }
  names_n<-NULL
  for(i in 1:ncol(indices)){names_n<-c(names_n,paste("UV",indices[1,i], indices[2,i],sep = ","))}
  
  if(is.null(indices_irt)){
    indices_irt<-matrix(c(sample(1:nrow(X[[1]]),min(round(nrow(X[[1]])/5),5), replace = FALSE),
                          sample(1:ncol(Y[[1]]),min(round(nrow(X[[1]])/5),5), replace = TRUE)), nrow=2, byrow = TRUE)
  }
  
  names_i<-NULL
  for(i in 1:ncol(indices_irt)){names_i<-c(names_i,paste("ThetaAlppha",indices_irt[1,i], indices_irt[2,i],sep = ","))}
  TAC<-matrix(nrow=0,ncol=ncol(indices_irt)) 
  UVC<-matrix(nrow=0,ncol=ncol(indices)) 
  
  colnames(UVC) <- names_n
  colnames(TAC) <- names_i
  
  
  VC<-matrix(nrow=0,ncol=2+length(c(seq(1,(K)*(K+1)/2), seq(1,(P)*(P+1)/2)))) 
  
  COV<-matrix(nrow = 0, ncol = P*K)
  
  colnames(VC) <- c(paste("Su",seq(1,(K)*(K+1)/2),sep=""),
                    paste("Stheta",seq(1,(P)*(P+1)/2),sep=""), "ve_connectivity","ve_attributes") 
  colnames(COV) <- paste("Sutheta",seq(1,(P)*K),sep="")
  
  
  
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  for (s in 1:(nscan + burn)) 
  { 
    
    #check this
    # update Tl
    
    if(is.null(H)){ETl<-sapply(1:length(Y), function(x) (as.numeric(b[x,]) + Theta), simplify = FALSE)}
    
    if(!is.null(H)){ETl<-sapply(1:length(Y), function(x) (as.numeric(b[x,]) + as.numeric(H[x,] %*% gamma) + Theta), simplify = FALSE)}
    


    # ETl=list()
    # for(i in 1:length(Y)){
    #   if(!is.null(H)){ETl[[i]] <- as.numeric(b[i,]) + as.numeric(H[i,] %*% gamma) + Theta}else{
    #     ETl[[i]] <- as.numeric(b[i,]) + Theta
    #   }
    #   #Tl[[i]] <- rTl_nrm(Tl[[i]], ETl[[i]],s1,Y[[i]])
    # }
    
    # update Fl
    if(!is.null(W)){EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + as.numeric(W[x,] %*% beta) + U %*% t(U)), simplify = FALSE)}
    if(is.null(W)){EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + U %*% t(U)), simplify = FALSE)}
    
 
    # EFl=list()
    # for(i in 1:length(X)){
    #   if(!is.null(W)){EFl[[i]] <- as.numeric(a[i,]) + as.numeric(W[i,] %*% beta) + U %*% t(U)}else{
    #     EFl[[i]] <- as.numeric(a[i,]) + U %*% t(U)
    #   }
    #   #Fl[[i]] <- rFl_nrm(Fl[[i]], EFl[[i]],s2,X[[i]])
    #   
    # }
    
    
    # update s2/s1
    s2<-rs2(Fl,offset = EFl, nu2=prior$nu2,s20=prior$s20)  
    s1<-rs1(Tl,offset = ETl, nu1=prior$nu1,s10=prior$s10)   
    
    

    
    tmp<-rbeta_a_fc_per(Fl,W=W,s2=s2,offset=U%*%t(U),ivA=prior$ivA,beta0=prior$beta0,S0=prior$S0) 
    beta <- tmp$beta 
    
    a<-tmp$a
    
    # update gamma,  b
    tmp1=rgamma_b_fc_per(Tl,H=H,s1=s1,offset=Theta,ivB=prior$ivB,gamma0=prior$gamma0,S1=prior$S1)
    gamma <- tmp1$gamma
    b <- tmp1$b
    
    
    ## update variances
    tmp <- rSu(cbind(U,Theta), Su0=prior$Sutheta0,etau=prior$etautheta) 
    Su = matrix(tmp$Su[1:K,1:K], nrow=K, ncol=K)
    Stheta = matrix(tmp$Su[(K+1):(K+P),(K+1):(K+P) ], nrow=P, ncol=P)
    Sutheta =matrix(tmp$Su[(K+1):(K+P),1:K ], nrow = P, ncol = K)
    S=tmp$Su
    
    
    
    
    
    
    # update U,V
    if (K > 0)
    {
      U<-rU(Fl,U,Theta, Stheta, Sutheta, Su, s2, offset=sapply(EFl, function(x) x-U%*%t(U), simplify = FALSE))
           
      U <- U -   colMeans(U)
      
      
       if(s ==1){U_target<-U}

       if(s>1){
         tmp <- Procrustes(U, U_target,
                           translate = FALSE,
                           dilate = FALSE,
                           sumsq = FALSE)
         U<-tmp$X.new



       }
      
    }
    
    
    
    
    # update Theta
    Theta <-rTheta(Tl, Theta ,U, Stheta, Sutheta, Su, s1, offset_Y=sapply(ETl, function(x) x-Theta, simplify = FALSE))
    
    
    
    Theta <- Theta -   colMeans(Theta)
    
    
    
    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  
      # save results
      
      BETA<-rbind(BETA, as.vector(beta)) 
      VC<-rbind(VC, c( Su[upper.tri(Su, diag = T)],
                       Stheta[upper.tri(Stheta, diag = T)], s2, s1)) 
      COV<-rbind(COV, as.vector(Sutheta))
      
      BETAPS<-BETAPS+beta
      GAMMAPS<-GAMMAPS+gamma
      
      
      
      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(U)
      
      tmp<-U %*% t(U)
      tmp.1<-NULL
      for(i in 1:ncol(indices)){tmp.1 <- c(tmp.1 ,tmp[indices[1,i],indices[2,i]])}
      UVC<-rbind(UVC, c(tmp.1))
      THETAPS <- THETAPS + Theta 
      
      tmp<-Theta 
      tmp.1<-NULL
      for(i in 1:ncol(indices_irt)){tmp.1 <- c(tmp.1 ,tmp[indices_irt[1,i],indices_irt[2,i]])}
      TAC<-rbind(TAC, c(tmp.1))
      APS <- APS + a
      BPS <- BPS + b 
      # 
      # Xs<-list()
      # Ys<-list()
      # for (i in 1:length(X)){
      #   Xs[[i]]<-simX_nrm(EFl[[i]],s2)
      #   Ys[[i]]<-simY_nrm(ETl[[i]],s1)
      #   # update posterior sum
      #   XPS[[i]]<-XPS[[i]]+Xs[[i]]
      #   YPS[[i]]<-YPS[[i]]+Ys[[i]]
      #   
      #   # save posterior predictive GOF stats
      #   #if(gof){ GOF[[i]]<-rbind(GOF[[i]],c(gofstats_c(Xs[[i]]),gofstats_a(Ys[[i]]))) }
      #   
      # }
    }
    
    
  } # end MCMC   
  
  # output 
  
  
  # posterior means 
  GAMMAPM<-GAMMAPS/nrow(VC)
  BETAPM<-BETAPS/nrow(VC)
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)
  UVPM<-UVPS/nrow(VC)
  THETAPM <-THETAPS/nrow(VC)
  
  # XPM<-sapply(XPS, function(x) x/nrow(VC), simplify = FALSE)
  # YPM<-sapply(YPS, function(x) x/nrow(VC), simplify = FALSE)
  
  
  if(!is.null(W)){EFlPM<-sapply(1:length(X), function(x) (APM[x,] + as.numeric(W[x,] %*% BETAPM) + UVPM), simplify = FALSE)}
  if(is.null(W)){EFlPM<-sapply(1:length(X), function(x) (APM[x,] + UVPM), simplify = FALSE)}
  
  if(!is.null(H)){ETlPM<-sapply(1:length(X), function(x) (BPM[x,] + as.numeric(H[x,] %*% GAMMAPM) + THETAPM), simplify = FALSE)}
  if(is.null(H)){ETlPM<-sapply(1:length(X), function(x) (BPM[x,] + THETAPM), simplify = FALSE)}
  
  
  # for(i in 1:length(X)){
  #   XPM[[i]]<-XPS[[i]]/nrow(VC) 
  #   YPM[[i]]<-YPS[[i]]/nrow(VC) 
  #   
  #   if(!is.null(W)){ EFlPM[[i]]<-APM[i,] + as.numeric(W[i,] %*% BETAPM) + UVPM }else(
  #     EFlPM[[i]]<-APM[i,] + UVPM
  #   )
  #   if(!is.null(H)){ ETlPM[[i]]<-BPM[i,] + as.numeric(H[i,] %*% GAMMAPM) + THETAPM }else{
  #     ETlPM[[i]]<-BPM[i,] + THETAPM
  #   }
  #   
  # }
  # 
  
  names(APM)<-names(BPM)<-names(X)
  rownames(UVPM)<-colnames(UVPM)<-rownames(THETAPM)<-rownames(X[[1]])
  
  
  # asymmetric output 
  UDV<-eigen(UVPM)
  U_1<-UDV$vectors[,seq(1,K,length=K)]%*%diag(sqrt(UDV$values[seq(1,K,length=K)]),nrow=K)
  
  rownames(U)<-rownames(X[[1]]) 
  
  fit <- list(BETAPM=BETAPM,GAMMAPM=GAMMAPM, VC=VC,COV=COV, APM=APM,BPM=BPM,U_1=U,U=U_1,UVPM=UVPM,THETAPM = THETAPM, EFlPM=EFlPM,ETlPM=ETlPM,
              X=X,Y=Y, UVC=UVC, TAC=TAC, Theta_1=Theta, input=input, indices=indices, indices_irt=indices_irt)
  
  class(fit) <- "abc"
  fit
}



