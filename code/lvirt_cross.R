#' LVIRT
#' 
#' An MCMC algorithm for fitting the latent variable item response theory (LVIRT) model 
#' 
#' 
#' This command provides posterior inference for parameters in LVIRT models 
#' assuming normal or binary data types
#' 
#' "nrm": assuming the data is normally distributed
#' 
#' "bin": assuming the data is binary
#' 
#' @usage lvirt(Y, family_responses, D=0,  model = twoparameter, 
#' seed = 1, nscan =
#' 10000, burn = 500, odens = 25, print = TRUE, gof=TRUE,
#' prior=list())
#' @param Y an N x M item response matrix. 
#' @param D integer: dimension of the item responses
#' @param family_responses character: "nrm" or "bin" for the item responses 
#' @param model character: "Rasch" or "twoparameter" 
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param prior list: A list of hyperparameters for the prior distribution
#' @return 
#' \item{VC}{posterior samples of the variance parameters}
#' \item{BETA}{posterior samples of the item intercept parameters}
#' \item{ALPHA}{posterior samples of the item discrimination parameters}
#' \item{THETA}{posterior samples of the latent person variable}
#' \item{AlphaTheta}{posterior mean of AlphaTheta} 
#'  \item{EH}{estimate of expectation of H
#' matrix}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' @author Selena Wang
#' @examples
#' 
#' @export jnirt
lvirt_cross<- function( Y,family_responses, D=1, 
               seed = 1, nscan = 10000, burn = 500, odens = 25,
               print = TRUE, gof=TRUE, indices_irt = NULL,
               prior=list())
{ 
  
  input<-list(D=D, nscan=nscan, burn=burn, odens=odens, prior=prior,indices_irt = indices_irt)
  # set random seed
  set.seed(seed)
  
  Y[is.na(Y)]=0
  
  # force binary if binary family specified
  if(is.element(family_responses,c("bin"))) { Y<-1*(Y>0) } 
  
  # g-prior setting for normal data 
  
  if(family_responses=="nrm") { H<-Y }
  if(family_responses=="bin")
  { 
    H<-matrix(zscores(Y),nrow(Y),ncol(Y)) 
    h01<- .5* ( max(H[Y==0],na.rm=TRUE) + min(H[Y==1],na.rm=TRUE) ) 
    H<-H - h01
  } 
  
  
  
  
  # starting values for missing entries 
  
  mui<-mean(H,na.rm=TRUE) 
  ai<-rowMeans(H,na.rm=TRUE) ; bi<-colMeans(H,na.rm=TRUE)  
  ai[is.na(ai)]<-0 ; bi[is.na(bi)]<-0 
  HA<-mui + outer(ai,bi,"+") 
  H[is.na(H)]<-HA[is.na(H)] 
  

  #### other starting values
  
  
  
  # beta
  beta<-colMeans(H,na.rm=TRUE)
  EH <- H- beta
  

  
  # s2, rho, s1
  s1 <-1
  if(family_responses=="nrm"){s1<-mean(EH^2)}
  

  
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
  BETA<- matrix(nrow = 0, ncol = ncol(Y))
  ALPHATHETAPS <- Theta %*% t(Alpha) * 0
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y) 
  gofXY<-c( gofstats_responses(Y))
  GOF<-matrix(gofXY,1,length(gofXY))  
  rownames(GOF)<-"obs"
  colnames(GOF)<-names(gofXY)
  

  VC<-matrix(nrow=0,ncol=1+length(c(seq(1,(D)*(D+1)/2)))) 

  colnames(VC) <- c(paste("Stheta",seq(1,(D)*(D+1)/2),sep=""), "ve_responses") 
  
   if(is.null(indices_irt)){
    indices_irt<-matrix(c(sample(1:nrow(X),5, replace = FALSE),
                          sample(1:ncol(Y),5, replace = FALSE)), nrow=2, byrow = TRUE)
  }
  
  names_i<-NULL
  for(i in 1:ncol(indices_irt)){names_i<-c(names_i,paste("ThetaAlppha",indices_irt[1,i], indices_irt[2,i],sep = ","))}
  TAC<-matrix(nrow=0,ncol=ncol(indices_irt)) 
  
  colnames(TAC) <- names_i
  

  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  for (s in 1:(nscan + burn)) 
  { 


    # update H
    EH = matrix(rep(1,nrow(X))) %*% t(beta)+ Theta%*%t(Alpha)

    if(family_responses=="nrm"){ H<-rH_nrm(H, EH,s1,Y) } 
    if(family_responses=="bin"){ H<-rH_bin(H, EH,Y) }
    

    
    # update s2
    if(family_responses=="nrm"){s1<-rs1(H-EH,nu1=prior$nu1,s10=prior$s10)}   


    
    ## update variances
    

    tmp <- rSu(Theta, Su0=prior$Sutheta0,etau=prior$etautheta) 
    Stheta=tmp$Su


    # update Theta

    Theta <-rTheta_IRT(H, beta, t(Alpha), t(Theta) , s1, Stheta)

    Theta <- Theta -   colMeans(Theta)
    Theta <- t(mhalf(solve(cov(Theta))) %*% t(Theta))

    
    # update item parameters
    tmp.2 <-rXi(H, beta, Alpha, Theta, prior$mxi, prior$Sigmaxi, s1)
    beta <-tmp.2$beta
    Alpha <-tmp.2$Alpha
    
    
    

    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  
      # save results

      BETA<-rbind(BETA, as.vector(beta)) 
      VC<-rbind(VC, c(Stheta[upper.tri(Stheta, diag = T)], s1)) 

      
      # update posterior sums of random effects
      ALPHATHETAPS <- ALPHATHETAPS + Theta %*% t(Alpha)

      # simulate from posterior predictive 
      ones=matrix(1,nrow(Y),1)
      EH<-ones %*% matrix(beta,nrow = 1, ncol=ncol(Y)) + Theta %*% t(Alpha)
      # 
      # 

      if(family_responses=="bin") { Ys<-simY_bin(EH) }
      if(family_responses=="nrm") { Ys<-simY_nrm(EH,s1) }

      # update posterior sum
      YPS<-YPS+Ys
      
      tmp<-Theta %*% t(Alpha)
      tmp.1<-NULL
      for(i in 1:5){tmp.1 <- c(tmp.1 ,tmp[indices_irt[1,i],indices_irt[2,i]])}
      TAC<-rbind(TAC, c(tmp.1))

      # save posterior predictive GOF stats
      if(gof){ Ys[is.na(Y)]<-NA ; GOF<-rbind(GOF,c(gofstats_responses(Ys))) }


      #print MC progress
      if (print) 
      {
        cat(s,round(apply(VC,2,mean),2),"\n")
        if (have_coda & nrow(VC) > 3 ) 
        {
          cat(round(coda::effectiveSize(VC)), "\n")
        }
      }
    }
    
    
  } # end MCMC   
  
  # output 
  

  # posterior means 

  ALPHATHETAPM <-ALPHATHETAPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  ones=matrix(1,nrow(Y),1)
  EH<-ones %*% apply(BETA,2,mean) + ALPHATHETAPM 
  

  colnames(ALPHATHETAPM)<-dimnames(Y)[[2]]
  dimnames(YPM)<-dimnames(EH)<-dimnames(Y)
  colnames(BETA)<-paste("easiness",seq(1,ncol(Y)),sep = "")
  
  # asymmetric output 

  ALPHATHETA<-svd(ALPHATHETAPM)
  Theta_1<-ALPHATHETA$u[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d[seq(1,D,length=D)]),nrow=D)
  Alpha_1<-ALPHATHETA$v[,seq(1,D,length=D)]%*%diag(sqrt(ALPHATHETA$d[seq(1,D,length=D)]),nrow=D)
    
  rownames(Theta) <- rownames(Y)
  rownames(Alpha)<-colnames(Y)
  fit <- list(BETA=BETA,VC=VC,ALPHATHETAPM = ALPHATHETAPM, EH=EH,
                YPM=YPM,  Theta_1=Theta_1, Alpha_1=Alpha_1,
                Theta=Theta, Alpha=Alpha, GOF=GOF,
                family_responses=family_responses,Y=Y, TAC=TAC, input=input)

  class(fit) <- "lvirt"
  fit
}



