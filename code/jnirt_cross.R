#' JNIRT
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
#' @param X an N x N square network. See family below for
#' various data types.
#' @param Y an M x M item response matrix. 
#' @param rvar logical: whether to fit row random effects
#' @param cvar logical: whether to  fit column random effects  
#' @param dcor logical: whether to fit a dyadic correlation 
#' @param K integer: dimension of the multiplicative effects (can be zero) in the social network
#' @param D integer: dimension of the item responses
#' @param family_network character: "nrm" or "bin" for the social network
#' @param family_responses character: "nrm" or "bin" for the item responses 
#' @param model character: "bilinear" or "UV" 
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param prior list: A list of hyperparameters for the prior distribution
#' @return \item{INT}{posterior samples of the intercept}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{COV}{posterior samples of the covariance parameters}
#' \item{BETA}{posterior samples of the item intercept parameters}
#' \item{ALPHA}{posterior samples of the item discrimination parameters}
#' \item{THETA}{posterior samples of the latent person variable}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} \item{V}{posterior mean of multiplicative column effects v}
#' \item{UVPM}{posterior mean of UV} 
#' \item{AlphaTheta}{posterior mean of AlphaTheta} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix}
#'  \item{EH}{estimate of expectation of H
#' matrix}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Selena Wang
#' @examples
#' 
#' @export jnirt
jnirt_cross<- function(X, Y,family_network,family_responses, K=2,D=1, model = "bilinear",
               rvar = FALSE , cvar = FALSE, dcor = FALSE, indices = NULL, indices_irt = NULL,
               seed = 1, nscan = 10000, burn = 500, odens = 25,
               print = TRUE, gof=TRUE, 
               prior=list())
{ 
  ## record model set up 
  
  input<-list(K=K, D=D, model=model, rvar=rvar, cvar=cvar,
              dcor=dcor, nscan=nscan, burn=burn, odens=odens, prior=prior, indices = indices, indices_irt = indices_irt)
  # set random seed
  set.seed(seed)
  X[is.na(X)]=0
  Y[is.na(Y)]=0
  
  # set diag to NA
  diag(X) <- NA 
  
  # force binary if binary family specified
  if(is.element(family_network,c("bin"))) { X<-1*(X>0) } 
  if(is.element(family_responses,c("bin"))) { Y<-1*(Y>0) } 
  
  # g-prior setting for normal data 
  if(family_network=="nrm" & is.null(prior$g))
  { 
    prior$g<-sum(!is.na(X))*var(c(X),na.rm=TRUE)
  }
  
  if(family_network=="bin" & model == "bilinear")
  {  
    xdist<-table(X)
    xmode<-as.numeric(names(xdist)[ xdist==max(xdist) ])[1] 
    ## eg, in a sparse binary network, xmode will be zero 
    XB<-1*(X!=xmode) 
    xbar<-mean(XB,na.rm=TRUE) ; mu<-qnorm(xbar)
    E<- (XB - xbar)/dnorm(qnorm(xbar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$Sutheta0)){ prior$Sutheta0<-diag(D+K) } 
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$etautheta)){ prior$etautheta<-(D+K+2) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(X)) }
  }
  
  if(family_network=="bin" & model == "UV")
  {  
    xdist<-table(X)
    xmode<-as.numeric(names(xdist)[ xdist==max(xdist) ])[1] 
    ## eg, in a sparse binary network, xmode will be zero 
    XB<-1*(X!=xmode) 
    xbar<-mean(XB,na.rm=TRUE) ; mu<-qnorm(xbar)
    E<- (XB - xbar)/dnorm(qnorm(xbar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$Sutheta0)){ prior$Sutheta0<-diag(D+K*2) } 
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$etautheta)){ prior$etautheta<-(D+K*2+2) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(X)) }
  }
  

  
  # starting Z values
  if(family_network=="nrm") { Z<-X }
  if(family_network=="bin")
  { 
    Z<-matrix(zscores(X),nrow(X),nrow(X)) 
    z01<- .5* ( max(Z[X==0],na.rm=TRUE) + min(Z[X==1],na.rm=TRUE) ) 
    Z<-Z - z01
  } 
  
  if(family_responses=="nrm") { H<-Y }
  if(family_responses=="bin")
  { 
    H<-matrix(zscores(Y),nrow(Y),ncol(Y)) 
    h01<- .5* ( max(H[Y==0],na.rm=TRUE) + min(H[Y==1],na.rm=TRUE) ) 
    H<-H - h01
  } 
  
  
  
  
  # starting values for missing entries 
  mu<-mean(Z,na.rm=TRUE) 
  a<-rowMeans(Z,na.rm=TRUE) ; b<-colMeans(Z,na.rm=TRUE)  
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  ZA<-mu + outer(a,b,"+") 
  Z[is.na(Z)]<-ZA[is.na(Z)] 
  
  mui<-mean(H,na.rm=TRUE) 
  ai<-rowMeans(H,na.rm=TRUE) ; bi<-colMeans(H,na.rm=TRUE)  
  ai[is.na(ai)]<-0 ; bi[is.na(bi)]<-0 
  HA<-mui + outer(ai,bi,"+") 
  H[is.na(H)]<-HA[is.na(H)] 
  

  #### other starting values
  
  # intercept
  intercept<-mean(Z)

  # a,b,Sab
  E<-Z-intercept 
  a<-rowMeans(E,na.rm=TRUE)*rvar ; b<-colMeans(E,na.rm=TRUE)*cvar
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
  
 
  
  # beta
  beta<-colMeans(H,na.rm=TRUE)
  EH <- H- beta
  

  
  # s2, rho, s1
  E<-E-outer(a,b,"+")  
  s1 <-1
  s2<-1
  if(family_network=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  if(family_responses=="nrm"){s1<-mean(EH^2)}
  
  if(rho==1){
    rho=.5
  }
  
  # U,V 
  U<-V<-matrix(0,nrow(X),K) 
  if(K>0) 
  {  
    sE<-svd(E)
    U<-sE$u[,1:K,drop=FALSE]%*%diag(sqrt(sE$d[1:K]),nrow=K)
    V<-sE$v[,1:K,drop=FALSE]%*%diag(sqrt(sE$d[1:K]),nrow=K)
    
    if (model == "bilinear"){
      V <-U
    }
  }
  
  # Alpha Theta
  
  Alpha<-matrix(0,ncol(Y),D) 
  Theta<-matrix(0,nrow(Y),D) 
  if(D>0) 
  {  
    sEH<-svd(EH)
    Theta<-sEH$u[,1:D,drop=FALSE]%*%diag(sqrt(sEH$d[1:D]),nrow=D)
    Alpha<-sEH$v[,1:D,drop=FALSE]%*%diag(sqrt(sEH$d[1:D]),nrow=D)
  }
  

  ## prior for item parameters
  

  if(is.null(prior$mxi)){ prior$mxi<-c(rep(1,D),0) }
  if(is.null(prior$Sigmaxi)){ prior$Sigmaxi<-diag(D+1) } 

  
  
  
  
  # output items
  INT <- matrix(nrow = 0, ncol = 1)
  BETA<- matrix(nrow = 0, ncol = ncol(Y))
  ALPHATHETAPS <- Theta %*% t(Alpha) * 0
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(X))  
  names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(X)
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y) 
  XPS<-matrix(0,nrow(X),ncol(X)) ; dimnames(XPS)<-dimnames(X) 
  gofXY<-c(gofstats_network(X), gofstats_responses(Y))
  GOF<-matrix(gofXY,1,length(gofXY))  
  rownames(GOF)<-"obs"
  colnames(GOF)<-names(gofXY)
  
  if(is.null(indices)){
    indices<-matrix(sample(1:nrow(X),10, replace = FALSE), nrow=2)
    
  }
  names_n<-NULL
  for(i in 1:ncol(indices)){names_n<-c(names_n,paste("UV",indices[1,i], indices[2,i],sep = ","))}
  
  if(is.null(indices_irt)){
    indices_irt<-matrix(c(sample(1:nrow(X),5, replace = FALSE),
                          sample(1:ncol(Y),5, replace = FALSE)), nrow=2, byrow = TRUE)
  }

  names_i<-NULL
  for(i in 1:ncol(indices_irt)){names_i<-c(names_i,paste("ThetaAlppha",indices_irt[1,i], indices_irt[2,i],sep = ","))}
  TAC<-matrix(nrow=0,ncol=ncol(indices_irt)) 
  UVC<-matrix(nrow=0,ncol=ncol(indices)) 
  
  colnames(UVC) <- names_n
  colnames(TAC) <- names_i
  
  if(model == "bilinear")
  {
    VC<-matrix(nrow=0,ncol=6+length(c(seq(1,(K)*(K+1)/2), seq(1,(D)*(D+1)/2)))) 
    
    COV<-matrix(nrow = 0, ncol = D*K)
    
    colnames(VC) <- c("va", "cab", "vb",
                      paste("Su",seq(1,(K)*(K+1)/2),sep=""),
                      paste("Stheta",seq(1,(D)*(D+1)/2),sep=""), "rho", "ve_network","ve_responses") 
    colnames(COV) <- paste("Sutheta",seq(1,(D)*K),sep="")
    
  }
  
  # names of parameters, symmetric case
  if(model == "UV")
  {
    VC<-matrix(nrow=0,ncol=6+length(c(seq(1,(2*K)*(2*K+1)/2), seq(1,(D)*(D+1)/2)))) 
    COV<-matrix(nrow = 0, ncol = D*K*2)
    
    colnames(VC) <- c("va", "cab", "vb",
                      paste("Su",seq(1,(2*K)*(2*K+1)/2),sep=""),
                      paste("Stheta",seq(1,(D)*(D+1)/2),sep=""), "rho", "ve_network","ve_responses") 
    colnames(COV) <- paste("Sutheta",seq(1,(D)*K*2),sep="")
  }    
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  for (s in 1:(nscan + burn)) 
  { 
    # update Z 

    EZ<-intercept + outer(a, b, "+") + U %*% t(V)
    if(family_network=="nrm"){ Z<-rZ_nrm(Z,EZ,rho,s2,X) } 
    if(family_network=="bin"){ Z<-rZ_bin(Z,EZ,rho,X) }
    

    # update H
    EH = matrix(rep(1,nrow(X))) %*% t(beta)+ Theta%*%t(Alpha)

    if(family_responses=="nrm"){ H<-rH_nrm(H, EH,s1,Y) } 
    if(family_responses=="bin"){ H<-rH_bin(H, EH,Y) }
    

    
    # update s2
    if(family_network=="nrm"){ s2<-rs2(Z-EZ,rho,nu0=prior$nu0,s20=prior$s20)  } 
    if(family_responses=="nrm"){s1<-rs1(H-EH,nu1=prior$nu1,s10=prior$s10)}   

    # update rho
    if(dcor){ rho <- rrho(Z-EZ,rho,s2,asp = prior$asp)} 
    

    # update beta, a b
    tmp <- rbeta_ab(Z-U%*%t(V),Sab,rho,
                    W=design_array(intercept=TRUE,n=nrow(Z)),s2,iV0=prior$iV0,m0=prior$m0,
                       g=prior$g) 
    intercept <- tmp$beta 
    a <- tmp$a * rvar
    b <- tmp$b * cvar 
    
    ## update variances
    
    if (model == "bilinear"){
        tmp <- rSu(cbind(U,Theta), Su0=prior$Sutheta0,etau=prior$etautheta) 
        Su = matrix(tmp$Su[1:K,1:K], nrow=K, ncol=K)
        Stheta = matrix(tmp$Su[(K+1):(K+D),(K+1):(K+D) ], nrow=D, ncol=D)
        Sutheta =matrix(tmp$Su[(K+1):(K+D),1:K ], nrow = D, ncol = K)
        S=tmp$Su

    }
    
    if (model == "UV"){
        tmp <- rSu(cbind(U,V,Theta), Su0=prior$Sutheta0,etau=prior$etautheta) 
        Su = matrix(tmp$Su[1:(K*2),1:(K*2)], nrow=K*2, ncol=K*2)
        Stheta = matrix(tmp$Su[(K*2+1):(K*2+D),(K*2+1):(K*2+D) ], nrow=D, ncol=D)
        Sutheta =matrix(tmp$Su[(K*2+1):(K*2+D),1:(K*2) ], nrow = D, ncol = K*2)
        S=tmp$Su

    }
    


    
    # update U,V
    if (K > 0)
    {
      if(model == "bilinear"){ 
        U<-rU(Z,U, t(Theta), Stheta, Sutheta, Su, s2,rho, offset=(outer(a,b,"+") + intercept) )
        U <- U -   colMeans(U)
        # if(s ==1){
        #   U_target <- U
        #   
        # }else{
        #   U <- rU_transformation(S=U,S_target=U_target) 
        # }
        V<-U
        
      }
      if(model == "UV")
      {
        
        ## need to check

        UV<-rUV(Z,Theta, U,V,S,rho,s2,offset=(intercept+outer(a,b,"+")))
        
        #UV<-rUV_network(Z,U,V,Suv=Su,rho,s2,offset=(intercept+outer(a,b,"+")))
        U<-UV$U ; V<-UV$V


      }

    }
    
    
    # update Theta
    if(model == "bilinear"){
      Theta <-rTheta(H, beta, t(Alpha), t(Theta), U, Stheta, Sutheta, Su, s1)

    }
    if(model == "UV"){
      Theta <-rTheta(H, beta, t(Alpha), t(Theta), cbind(U,V), Stheta,Sutheta, Su, s1)
    }
    
    Theta <- Theta -   colMeans(Theta)
    Theta <- t(mhalf(solve(cov(Theta))) %*% t(Theta))

    
    # update item parameters
    tmp.2 <-rXi(H, beta, Alpha, Theta, prior$mxi, prior$Sigmaxi, s1)
    beta <-tmp.2$beta
    Alpha <-tmp.2$Alpha
    
    
    # update Sab - full SRM
    if(rvar & cvar)
    {
      Sab<-rSab(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
    }
    
    # update Sab - rvar only
    if (rvar & !cvar ) 
    {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(X))/2, (1 + sum(a^2))/2)
    }
    
    # update Sab - cvar only
    if (!rvar & cvar ) 
    {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(X))/2, (1 + sum(b^2))/2)
    }
    
   

    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  
      # save results

      INT<-rbind(INT, intercept)
      BETA<-rbind(BETA, as.vector(beta)) 
      VC<-rbind(VC, c(Sab[upper.tri(Sab, diag = T)], Su[upper.tri(Su, diag = T)],
                        Stheta[upper.tri(Stheta, diag = T)], rho,s2, s1)) 
      COV<-rbind(COV, as.vector(Sutheta))
      
      

      
      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      
      tmp<-U %*% t(V)
      tmp.1<-NULL
      for(i in 1:5){tmp.1 <- c(tmp.1 ,tmp[indices[1,i],indices[2,i]])}
      UVC<-rbind(UVC, c(tmp.1))
      ALPHATHETAPS <- ALPHATHETAPS + Theta %*% t(Alpha)
      
      tmp<-Theta %*% t(Alpha)
      tmp.1<-NULL
      for(i in 1:5){tmp.1 <- c(tmp.1 ,tmp[indices_irt[1,i],indices_irt[2,i]])}
      TAC<-rbind(TAC, c(tmp.1))
      APS <- APS + a
      BPS <- BPS + b 
      
      # simulate from posterior predictive 
      EZ<-intercept + outer(a, b, "+") + U %*% t(V) 
      ones=matrix(1,nrow(X),1)
      EH<-ones %*% matrix(beta,nrow = 1, ncol=ncol(Y)) + Theta %*% t(Alpha)
      # 
      # 
      # need to add
      if(family_network=="bin") { Xs<-simX_bin(EZ,rho) }
      if(family_network=="nrm") { Xs<-simX_nrm(EZ,rho,s2) }
      
      if(family_responses=="bin") { Ys<-simY_bin(EH) }
      if(family_responses=="nrm") { Ys<-simY_nrm(EH,s1) }

      # update posterior sum
      XPS<-XPS+Xs
      YPS<-YPS+Ys

      # save posterior predictive GOF stats
      if(gof){ Xs[is.na(X)]<-NA ; Ys[is.na(Y)]<-NA; GOF<-rbind(GOF,c(gofstats_network(Xs),gofstats_responses(Ys))) }


      #print MC progress
      if (print) 
      {
        cat(s,round(mean(INT),2),":",round(apply(VC,2,mean),2),"\n")
        if (have_coda & nrow(VC) > 3 & length(INT)>0) 
        {
          cat(round(coda::effectiveSize(INT)), "\n")
        }
      }
    }
    
    
  } # end MCMC   
  
  # output 
  

  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)
  UVPM<-UVPS/nrow(VC)
  ALPHATHETAPM <-ALPHATHETAPS/nrow(VC)
  XPM<-XPS/nrow(VC) 
  YPM<-YPS/nrow(VC) 
  EZ<-mean(INT) + outer(APM,BPM,"+")+UVPM  
  ones=matrix(1,nrow(X),1)
  EH<-ones %*% apply(BETA,2,mean) + ALPHATHETAPM 
  
  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(X)[[1]]
  
  rownames(ALPHATHETAPM)<-dimnames(X)[[1]]
  colnames(ALPHATHETAPM)<-dimnames(Y)[[2]]
  dimnames(XPM)<-dimnames(EZ)<-dimnames(X)
  dimnames(YPM)<-dimnames(EH)<-dimnames(Y)
  colnames(BETA)<-paste("easiness",seq(1,ncol(Y)),sep = "")
  
  # asymmetric output 
  if(model == "UV") 
  {
    UDV<-svd(UVPM)
    U_1<-UDV$u[,seq(1,K,length=K)]%*%diag(sqrt(UDV$d[seq(1,K,length=K)]),nrow=K)
    V_1<-UDV$v[,seq(1,K,length=K)]%*%diag(sqrt(UDV$d[seq(1,K,length=K)]),nrow=K)
    
    ALPHATHETA<-svd(ALPHATHETAPM)
    Theta_1<-ALPHATHETA$u[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d[seq(1,D,length=D)]),nrow=D)
    Alpha_1<-ALPHATHETA$v[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d[seq(1,D,length=D)]),nrow=D)
    
    rownames(Theta) <- rownames(U)<-rownames(V)<-rownames(X) 
    rownames(Alpha)<-colnames(Y)
    fit <- list(BETA=BETA,INT = INT, VC=VC,COV=COV, APM=APM,BPM=BPM,U=U,V=V,U_1=U_1,V_1=V_1,UVPM=UVPM,ALPHATHETAPM = ALPHATHETAPM, EZ=EZ,EH=EH,
                XPM=XPM,YPM=YPM,  Theta_1=Theta_1, Alpha_1=Alpha_1,
                Theta=Theta, Alpha=Alpha, GOF=GOF,family_network=family_network,
                family_responses=family_responses,X=X,Y=Y, UVC=UVC,TAC=TAC, input=input)
  }
  if(model == "bilinear") 
  {
    UDV<-eigen(UVPM)
    U_1<-UDV$vectors[,seq(1,K,length=K)]%*%diag(sqrt(UDV$values[seq(1,K,length=K)]),nrow=K)
    
    ALPHATHETA<-svd(ALPHATHETAPM)
    Theta_1<-ALPHATHETA$u[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d)[seq(1,D,length=D)],nrow=D)
    Alpha_1<-ALPHATHETA$v[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d)[seq(1,D,length=D)],nrow=D)
    
    rownames(Theta) <- rownames(U)<-rownames(X) 
    rownames(Alpha)<-colnames(Y)
    fit <- list(BETA=BETA,INT = INT, VC=VC,COV=COV, APM=APM,BPM=BPM,U_1=U_1,U=U,UVPM=UVPM,ALPHATHETAPM = ALPHATHETAPM, EZ=EZ,EH=EH,
                XPM=XPM,YPM=YPM, Theta_1=Theta_1, Alpha_1=Alpha_1,Theta=Theta, Alpha=Alpha,GOF=GOF,family_network=family_network,
                family_responses=family_responses,X=X,Y=Y, UVC=UVC, TAC=TAC, input=input)
  }
  
  class(fit) <- "jnirt"
  fit
}



