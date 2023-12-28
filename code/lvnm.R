#' lvnm
#' 
#' An MCMC algorithm for fitting the latent variable network model 
#' 
#' This command provides posterior inference for parameters in a latent variable network model
#' assuming normal or binary data types
#' 
#' "nrm": assuming the data is normally distributed
#' 
#' "bin": assuming the data is binary
#' 
#' @usage lvnm(X, Y, family_network, family_responses, K=0, rvar = FALSE ,
#' cvar = FALSE,  dcor = FALSE, model = UU, 
#' intercept=TRUE, seed = 1, nscan =
#' 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE,
#' prior=list())
#' @param X an N x N square network. See family below for
#' various data types.
#' @param rvar logical: whether to fit row random effects
#' @param cvar logical: whether to  fit column random effects  
#' @param dcor logical: whether to fit a dyadic correlation 
#' @param K integer: dimension of the multiplicative effects (can be zero) in the social network
#' @param family_network character: "nrm" or "bin" for the social network
#' @param model character: "bilinear" or "UV" 
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param gof logical: calculate goodness of fit statistics?
#' @param odens output density for the Markov chain
#' @param print logical: print results while running?
#' @param prior list: A list of hyperparameters for the prior distribution
#' @return \item{INT}{posterior samples of the intercept}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} \item{V}{posterior mean of multiplicative column effects v}
#' \item{UVPM}{posterior mean of UV} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Selena Wang
#' @examples
#' 
#' @export lvnm
lvnm<-function (X, family_network, K=2, model = "bilinear",
                rvar = FALSE , cvar = FALSE, dcor = FALSE, independent=FALSE,
                seed = 1, nscan = 10000, burn = 500, odens = 25,
                print = TRUE, gof=TRUE, indices = NULL,
                prior=list())
{ 
  input<-list(K=K, model=model, rvar=rvar, cvar=cvar,independent=independent,
              dcor=dcor, nscan=nscan, burn=burn, odens=odens, prior=prior, indices = indices)
  # set random seed
  set.seed(seed)
  
  
  # set diag to NA
  diag(X) <- NA 
  
  # force binary if binary family specified
  if(is.element(family_network,c("bin"))) { X<-1*(X>0) } 
  
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
    if(is.null(prior$Su0)){ prior$Su0<-diag(K) *vscale } 
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$etau1)){ prior$etau1<-round((K+2) *vdfmlt) }  
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
    if(is.null(prior$Su0)){ prior$Su0<-diag(K*2) *vscale } 
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$etau1)){ prior$etau1<-round((K*2+2) *vdfmlt) }  
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
  
  
  
  # starting values for missing entries 
  mu<-mean(Z,na.rm=TRUE) 
  a<-rowMeans(Z,na.rm=TRUE) ; b<-colMeans(Z,na.rm=TRUE)  
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  ZA<-mu + outer(a,b,"+") 
  Z[is.na(Z)]<-ZA[is.na(Z)] 
  
  
  #### other starting values
  
  # intercept
  intercept<-mean(Z)
  
  # a,b,Sab
  E<-Z-intercept 
  a<-rowMeans(E,na.rm=TRUE)*rvar ; b<-colMeans(E,na.rm=TRUE)*cvar
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
  
  
  # s2, rho, s1
  E<-E-outer(a,b,"+")  
  s2<-1
  if(family_network=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
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
  

  
  # output items
  INT <- matrix(nrow = 0, ncol = 1)
  UVPS <- U %*% t(V) * 0 
  APS<-BPS<- rep(0,nrow(X))  
  names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(X)
  XPS<-matrix(0,nrow(X),ncol(X)) ; dimnames(XPS)<-dimnames(X) 
  gofX<-c(gofstats_network(X))
  GOF<-matrix(gofX,1,length(gofX))  
  rownames(GOF)<-"obs"
  colnames(GOF)<-names(gofX)
  
  if(is.null(indices)){
    indices<-matrix(sample(1:nrow(X),10, replace = FALSE), nrow=2)
    
  }
  names<-NULL
  for(i in 1:ncol(indices)){names<-c(names,paste("VectorProduct",indices[1,i], indices[2,i],sep = ","))}
  
  UVC<-matrix(nrow=0,ncol=ncol(indices)) 
  colnames(UVC) <- names
  
  
  if(model == "bilinear")
  {
    VC<-matrix(nrow=0,ncol=5+(K)*(K+1)/2) 
    
    colnames(VC) <- c("va", "cab", "vb",
                      paste("Su",seq(1,(K)*(K+1)/2),sep=""), "rho", "ve_network") 
  }
  
  if(model == "UV")
  {
    VC<-matrix(nrow=0,ncol=5+(2*K)*(2*K+1)/2) 
    
    colnames(VC) <- c("va", "cab", "vb",
                      paste("Su",seq(1,(2*K)*(2*K+1)/2),sep=""), "rho", "ve_network") 
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
    
    
    # update s2
    if(family_network=="nrm"){ s2<-rs2(Z-EZ,rho,nu0=prior$nu0,s20=prior$s20)  } 
    
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
      
      
      if(independent){
        tmp <- rSu_inp(U)
        
      }else{
        tmp <- rSu(U, Su0=prior$Su0,etau=prior$etau1)
        
      }
      Su <-matrix(tmp$Su, nrow = K, ncol = K)
      
      
      
    }
    
    if (model == "UV"){
      
      if(independent){
        tmp <- rSu_inp(cbind(U,V))
        
      }else{
        tmp <- rSu(cbind(U,V), Su0=prior$Su0,etau=prior$etau1)
        
      }
      
      Su <-matrix(tmp$Su, nrow = K*2, ncol = K*2)
      
      
      
    }
    
    
    
    
    
    # update U,V
    if (K > 0)
    {
      if(model == "bilinear"){ 
        U<-rU_network(Z,U,Su, s2=s2,rho, offset=(outer(a,b,"+") + intercept))
        U <- U -   colMeans(U)
        V<-U
        
      }
      if(model == "UV")
      {
        
        UV<-rUV_network(Z,U,V,Suv=Su,rho,s2,offset=(intercept+outer(a,b,"+")))
        U<-UV$U ; V<-UV$V
        
        U <- U -   colMeans(U)
        V <- V -   colMeans(V)
   

      }
      
    }
    
    
    
    # update Sab 
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
      VC<-rbind(VC, c(Sab[upper.tri(Sab, diag = T)], Su[upper.tri(Su, diag = T)], rho, s2)) 
      
      
      # update posterior sums of random effects
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b 
      
      tmp<-U %*% t(V)
      
      tmp.1<-NULL
      
      for(i in 1:ncol(indices)){tmp.1 <- c(tmp.1 ,tmp[indices[1,i],indices[2,i]])}
      
      UVC<-rbind(UVC, c(tmp.1))
      
      # simulate from posterior predictive 
      EZ<-intercept + outer(a, b, "+") + U %*% t(V) 
      
      
      # 
      if(family_network=="bin") { Xs<-simX_bin(EZ,rho) }
      if(family_network=="nrm") { Xs<-simX_nrm(EZ,rho,s2) }
      
      # update posterior sum
      XPS<-XPS+Xs
      
      # save posterior predictive GOF stats
      if(gof){ Xs[is.na(X)]<-NA ; GOF<-rbind(GOF,c(gofstats_network(Xs))) }
      # 
      # print MC progress 
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
  XPM<-XPS/nrow(VC) 
  
  
  EZ<-mean(INT) + outer(APM,BPM,"+")+UVPM  
  
  
  names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(X)[[1]]
  
  
  dimnames(XPM)<-dimnames(EZ)<-dimnames(X)
  
  # asymmetric output 
  if(model == "UV") 
  {
    UDV<-svd(UVPM)
    U_1<-UDV$u[,seq(1,K,length=K)]%*%diag(sqrt(UDV$d)[seq(1,K,length=K)],nrow=K)
    V_1<-UDV$v[,seq(1,K,length=K)]%*%diag(sqrt(UDV$d)[seq(1,K,length=K)],nrow=K)
    
    rownames(U)<-rownames(V)<-rownames(X) 
    fit <- list(INT = INT, VC=VC, APM=APM,BPM=BPM,U=U,V=V,U_1=U_1,V_1=V_1,UVPM=UVPM,
                EZ=EZ, XPM=XPM, GOF=GOF,family_network=family_network,X=X, UVC=UVC, input=input)
  }
  if(model == "bilinear") 
  {
    UDV<-eigen(UVPM)
    U_1<-UDV$vectors[,seq(1,K,length=K)]%*%diag(sqrt(UDV$values[seq(1,K,length=K)]),nrow=K)
    
    rownames(U)<-rownames(X) 
    fit <- list(INT = INT, VC=VC, APM=APM,BPM=BPM,U_1=U_1,U=U,UVPM=UVPM,EZ=EZ,
                XPM=XPM,GOF=GOF,family_network=family_network,X=X, UVC=UVC, input=input)
  }
  
  class(fit) <- "lvnm"
  fit
}



