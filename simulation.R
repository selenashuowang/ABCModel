library(foreign)
library(dplyr)
library(gdata)
library(psych)
library(matrixcalc)

library(MASS)
library(expm)
library(boot)
library(Matrix)

library(CovTools)
library(ggplot2)
library(plyr)
library(lvm4net)
library(pROC)

library(sbm)

sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) {
  files <- list.files(folderName, full.names=TRUE)
  
  # Grab only R files
  files <- files[ grepl("\\.[rR]$", files) ]
  
  if (!length(files) && showWarnings)
    warning("No R files in ", folderName)
  
  for (f in files) {
    if (verbose)
      cat("sourcing: ", f, "\n")
    ## TODO:  add caught whether error or not and return that
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose)
  }
  return(invisible(NULL))
}

## change it to code folder
sourceEntireFolder("/gpfs/gibbs/project/zhao_yize/sw2384/code", verbose=FALSE, showWarnings=TRUE)

#### args 1) generate data N 50 versus 200
N<-10

## args 2) two levels of V, 20 versus 70
V<-20
#P<-2
P<-1
K<-2

# W<-matrix(1,N,1)
# H<-matrix(1,N,1)

W<-NULL
H<-NULL

a_t<-matrix(c(-1,1), nrow = N, ncol = 1)
b_t<-matrix(c(1,-1), nrow = N, ncol = 1)


## args 4) Covariance of U and theta

A <- diag(1,K+P)
A[1,3]=0.9
A[3,1]=0.9



seed <- 100

set.seed(seed)

UTheta <- mvrnorm(n = V, mu=rep(0,K+P), Sigma=A)
U_t <- UTheta[,1:K]
Theta_t <- UTheta[,(K+1):(K+P)]


#U_t<-mvrnorm(n = V, mu=rep(0,K), Sigma=diag(1,K))
#Theta_t<-mvrnorm(n = V, mu=rep(0,P), Sigma=diag(2,P))

# beta_t<-matrix(1,nrow = ncol(W),ncol=1)
# gamma_t<-matrix(2,nrow = ncol(H),ncol=1)

beta_t=NULL
gamma_t=NULL

#connectivity variance, five levels 0.5,1,5
s2_t=as.numeric(0.5)
#attribute variance
s1_t=0.5

N
V
s2_t
s1_t
A
seed

X<-list()
Y<-list()
for(i in 1:N){
  
  errormatrix=matrix(0, nrow = V, ncol = V)
  errormatrix[upper.tri(errormatrix,diag = FALSE)]=rnorm(V*(V-1)/2, sd=sqrt(s2_t))
  errormatrix=t(errormatrix)+errormatrix
  diag(errormatrix)=rnorm(V, sd=sqrt(s2_t))
  
  X[[i]]<-as.numeric(a_t[i,])  + U_t%*% t(U_t) +errormatrix
  diag(X[[i]])=NA
  Y[[i]]<-as.numeric(b_t[i,])  + Theta_t +matrix(rnorm(V*P, sd=sqrt(s1_t)),V,P)
  
}


true_Uv=U_t%*% t(U_t)
obs_x=true_Uv[upper.tri(true_Uv, diag = FALSE)]


names(X)=as.character(seq(1,length(X)))
names(Y)=as.character(seq(1,length(X)))
ids=names(X)


sampled.id=sample(ids, round(length(ids)*.5), replace = FALSE)

train.id=ids[!ids %in% sampled.id]

X.array=simplify2array(X)

X.array[,,sampled.id]=NA

X_pred=lapply(seq(dim(X.array)[3]), function(x) X.array[ , , x])


## fit model
df=abc(X=X_pred, Y=Y,W=NULL, H=NULL, K=2,
            indices = NULL, indices_irt = NULL,
            seed = 1, nscan = 5000, burn = 100, odens = 10,
            print = FALSE, gof=FALSE, plot=FALSE,
            prior=list())






