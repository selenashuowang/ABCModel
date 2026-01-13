# ABC Model tutorial
## Overview


The Attributes-informed Brain Connectivity Model (ABC) model is specifically designed 



## Installation

This toolkit is implemented in R. Follow these steps for setup:


1. Clone or download the repository to your local machine.
2. Open R and navigate to the directory containing the toolkit.
3. Run the ABC model.



## Required Data


The toolkit is designed to analyze data consisting of brain connectivity and region level attribute information. Specifically, it requires two key files:

* `X`: a list of $V \times V$ brain connectivity data.
* `Y`: a list of $V \times P$ attribute data.

Additionally, simulated example data is available in the directory `data/X.RData` and `data/Y.RData` for demonstration purposes.


## Key parameters

* `W` a matrix of $N \times Q$ covariates for the connectivity data.
* `H` a matrix of $N \times Q1$ covariates for the attribute data.
* `K` number determining the latent dimension of multiplicative effects
* `nscan` number of iterations of the Markov chain (beyond burn-in)
* `burn` burn-in for the Markov chain

Note that sufficient burn-in is need to reach optimal covariance parameter estimates. See details in the method paper. 

(Add Note on latent dimension)

## Usage


The main functionality of the ABC Toolkit is encapsulated in the `abc.r` script, which performs MCMC estimation of each --------. Example data are in the `data` folder. The `abc.r` can do the following connectivity data analysis:

### Run ABC model

``` {r}
library(abc.r)

setwd("./data") #directory of example data 

load(file='X.rda')
load(file='Y.rda')


model1=abc(X, Y,W=NULL, H=NULL, K = 2,
                   seed = 1, nscan = 10, burn = 1, odens = 1,
                   prior=list())

```


### Obtain unbiased estimates of covariance parameters.

From the saved model results, we can obtain estimated mean of $V \times V$ brain latent connectivity as `model1$UVPM`, posterior samples of the covariance parameters can be obtained via `model1$COV`, etc. See details in the pacakge documentation. 

ABC Simulation Example
================

## Introduction

This is a simulation study using the ABC package to analyze brain
connectivity data and attribute information. The simulation includes data
generation, model fitting, and result visualization. The goal is to
demonstrate how latent space models can be used to analyze the
relationship between brain connectivity patterns and region level attribute measures.

## Load prerequisite packages

``` r
library(MASS)
library(psych)
library(coda)
library(magic)
```

    ## Loading required package: abind

## Source ABC code

We first define a helper function to source all R files from the abc package directory. This function will load all the necessary functions for our analysis.

``` r
#set seed and load the latentSNA code
set.seed(18)

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

sourceEntireFolder("code", verbose=FALSE, showWarnings=TRUE)  
```

## Generate Simulated Data

In this section, we generate simulated data that mimics brain
connectivity patterns and behavioral outcomes. We set up:

- A sample of 1000 subjects
- 20 brain regions
- 8 significant brain regions with strong correlations (0.9)
- A single region attribute 
- One-dimensional latent space

First, we set up basic parameters:

``` r
N<-1000 #N is sample size
ids=seq(1,N)

V<-20 # V is number of brain regions
P<-1 # P is number of attributes

K<-1 # K and D are latent space
D<-1

W<-NULL# W and H are covariates
H<-NULL
```

Next, we create the covariance structure with significant brain regions:

``` r
a_t<-matrix(0, nrow = N, ncol = 1)
S <- diag(1,V+D)
n_signa=V

id_siga=c(1,2,3,4,5,6,7,8) # sepcify significant brain regions

id=list()
id[[1]]=id

# Set correlation structure for significant regions
for (each in id_siga){
  
  S[each,id_siga[!id_siga %in% each]]=.9
}

S[(V+1),id_siga]=.9
S[id_siga,(V+1)]=.9

# Extract submatrices
Su = matrix(S[1:V,1:V], nrow=V, ncol=V)
Stheta = matrix(S[(V+1):(V+D),(V+1):(V+D) ], nrow=D, ncol=D)
Sutheta =matrix(S[(V+1):(V+D),1:V], nrow = D, ncol = V)

# Print Sutheta to show that brain regions 1-8 are significantly associated with behavior (correlation = 0.9)
print(Sutheta)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16]
    ## [1,]  0.9  0.9  0.9  0.9  0.9  0.9  0.9  0.9    0     0     0     0     0     0     0     0
    ##      [,17] [,18] [,19] [,20]
    ## [1,]     0     0     0     0

Generate latent variables using multivariate normal distribution:

``` r
UTheta <- mvrnorm(n = N, mu=rep(0,(V+D)), Sigma=S, empirical = FALSE)
U.array=array(NA, dim = c(V,K,N))
U.array[,1,]=t(UTheta[,1:(V)])
Theta_t <- data.matrix(UTheta[,(V+1):(V+D)])
rownames(Theta_t)=ids
```

Set model parameters and generate connectivity matrices:

``` r
beta_t=NULL
gamma_t=NULL

#connectivity variance
s2_t=0.1
#attribute variance
s1_t=0.1

Alpha_t=matrix(1, nrow = P, ncol = 1)
theoretical.str=Alpha_t
b_t=matrix(0, nrow = P, ncol = 1)

X<-list()
for(i in 1:N){
  
  errormatrix=matrix(0, nrow = V, ncol = V)
  errormatrix[upper.tri(errormatrix,diag = FALSE)]=rnorm(V*(V-1)/2, sd=sqrt(s2_t))
  errormatrix=t(errormatrix)+errormatrix
  diag(errormatrix)=rnorm(V, sd=sqrt(s2_t))
  
  X[[i]]<-as.numeric(a_t[i,])  + U.array[,,i]%*% t(U.array[,,i]) +errormatrix
  #diag(X[[i]])=NA
  
}
```

Finally, generate attribute information, with significant brain regions influencing values:

``` r
# Generate region-level attribute information
Y <- vector("list", N)

for (i in 1:N) {
  
  # Mean structure replicated across V brain regions
  mean_i <- matrix(rep(b_t, each = V), nrow = V, ncol = P) +
            matrix(rep(Theta_t[i, ] %*% t(Alpha_t), each = V),
                   nrow = V, ncol = P)
  
  # Add region-specific noise
  Y[[i]] <- mean_i + matrix(
    rnorm(V * P, sd = sqrt(s1_t)),
    nrow = V,
    ncol = P
  )
}

```
## Fit ABC Model and Store Results

Here we show the model fitting process for completeness. However, since
this process takes several hours to run, we have pre-computed the
results and saved them for analysis. The process includes:

1.  Data Preparation:
    - Assigning proper names to observations in X and Y
    - Creating test and training sets (10% test, 90% training)
    - Storing the full dataset for later comparison
2.  Model Specification:
    - MCMC parameters: 15,000 iterations with 15,000 burn-in
    - Random seed set to 1 for reproducibility
    - No additional covariates (W=NULL, H=NULL)

To ensure robust results, we recommend running 10 independent
simulations and selecting the best-performing one based on convergence
diagnostics and model fit metrics.

Here is the code for reference:

``` r
train_ratio <- 0.8
n <- length(X)
train_indices <- sample(seq_len(n), size = floor(train_ratio * n))
test_indices <- setdiff(seq_len(n), train_indices)

X_train <- X[train_indices]
X_test <- X[test_indices]
Y_train <- Y[train_indices]
Y_test <- Y[test_indices]
md=abc(X=X_train, Y=Y_train,W=NULL, H=NULL, K=1, nscan = 15000, burn = 1000)

res=list("model"=md,'testX' = X_test,'testY' = Y_test)
saveRDS(res,'simulationABC.rds')
```
## Explanation of Results

The abc model uses MCMC to estimate parameters, producing several
key outputs:

``` r
result<-readRDS("simulation.rds")
result$model$UVC |> dim()
```
500 190

``` r
result$model$UVPM |> dim()
```
20 20 
1.  UVC matrix (600 x 190 ): Provides all connectivity esitmates for connectivity edge  
2.  UVPM (20 x 20): Contains the scaled estimated connectivity for each edge.
3.  EFlPM (800 x 1): Provides estimated connectivity values for each participant in our training set 

While the model produces other variables, we focus on these three key outputs as they are most crucial for evaluating brain-connectivity trends.

## Determining model performance/fit

Model fit is determined via correlation between individual estimated connectivity and the connectivity values in our 20% testing set. This correlation can be used for selecting optimal latent space dimensions and determining best performing models for final analysis using k-fold cross validation. 

Model correlation can be determined via: 

```r
model1 <- result$model
x_test <- result$testX
l <- length(x_test)
est <- model1$EFlPM
est_split <- est[1:l]
vec1 <- unlist(x_test)
vec2 <- unlist(est_split)
correlations <- cor(vec1, vec2)
```
The correlation for our simulated data:
```r
correlations
```
[1]  0.4929648

This value indicates quality model fit.
## Plotting Estimated connectivity matrix

```r
library(corrplot)
corrplot(result$model$UVPM,
        method = "color",
        t.pos = "n", #Removes Axis Numbers
        cl.cex = 1.25, #Size of Legend Text
        is.corr = FALSE
         )
```
![](examplefiles/CorrplotSimulation.png)

## Detmining Credible Intervals

The `UVC` value presents each connectivity estimate for each edge ($V * (V-1)/2$). These estimates can be used to compute a 95% credible interval for each edge. This credible interval can be compared with models ran using different populations to compared populations. 

```r
hi.g1 <- t(apply(result$model$UVC, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
#Returns a $V * (V-1)/2$ x 2 matrix which details the lower and upper bounds of the 95% credible interval. 
```
This credible interval values is currently in the form of a vector, which does not tell us which edge belongs to. We can reformat this into a VxV matrix, with each index storing the lower and upper bound 95% CI. 

```r
loc_matrix <- matrix(0, nrow = V, ncol = V)
  
  # Start filling the upper triangle downwards, column by column
    index <- 1
    for (j in 2:V) { # Start from the second column
      for (i in 1:(j - 1)) { # Only fill below the diagonal
        loc_matrix[i, j] <- index
        index <- index + 1
     }
    }

CIArray <- array(NA, dim = c(V, V, 2))

for (i in 1:V) {
  for (j in 1:V) {
    n <- loc_matrix[i, j]
    if (n == 0) next
    
    CIArray[i, j, ] <- hi.g1[n, ]
  }
}
```
The lowerbound can be accessed via `CIArray[i,j,1]` and the upper bound can be accessed via `CIArray[i,j,2]`

This CIMatrix can be compared with CI's of different populations, with non-overlapping CI's being considered a statistically significant difference for that edge. 