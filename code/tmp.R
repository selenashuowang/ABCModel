zscores<-function(y,ties.method="average")
{
  qnorm( rank(y,na.last="keep",ties.method=ties.method)/(1+sum(!is.na(y))) )
}



Procrustes <-function (X, Xstar, translate = FALSE, dilate = FALSE, sumsq = FALSE) 
{
  if (any(!is.logical(translate), length(translate) != 1)) 
    stop("'translate' must be a single logical indicator", 
         call. = FALSE)
  if (any(!is.logical(dilate), length(dilate) != 1)) 
    stop("'dilate' must be a single logical indicator", call. = FALSE)
  if (any(!is.logical(sumsq), length(sumsq) != 1)) 
    stop("'sumsq' must be a single logical indicator", call. = FALSE)
  if ((N <- nrow(X)) != nrow(Xstar)) 
    stop("X and Xstar do not have the same number of rows", 
         call. = FALSE)
  if ((P <- ncol(X)) != (P2 <- ncol(Xstar))) {
    if (P < P2) {
      warning("X padded out to same number of columns as Xstar\n", 
              call. = FALSE, immediate. = TRUE)
      X <- cbind(X, matrix(0L, nrow = N, ncol = P2 - P))
    }
    else stop("X cannot have more columns than Xstar", call. = FALSE)
  }
  if (P2 == 0) 
    stop("Xstar must contain at least one column", call. = FALSE)
  if (anyNA(Xstar) || anyNA(X)) 
    stop("X and Xstar are not allowed to contain missing values", 
         call. = FALSE)
  J <- if (translate) 
    diag(N) - 1/N
  else diag(N)
  C <- if (translate) 
    crossprod(Xstar, J) %*% X
  else crossprod(Xstar, X)
  if (!all(P == 1, P2 == 1)) {
    svdX <- svd(C)
    R <- tcrossprod(svdX$v, svdX$u)
  }
  else R <- 1L
  d <- if (dilate) 
    sum(C * R)/sum(crossprod(J, X) * X)
  else 1L
  tt <- if (translate) 
    crossprod(Xstar - d * X %*% R, matrix(1L, N, 1))/N
  else 0L
  X.new <- d * X %*% R + if (translate) 
    matrix(tt, N, P2, byrow = TRUE)
  else tt
  return(c(list(X.new = X.new), list(R = R), if (translate) list(t = tt), 
           if (dilate) list(d = d), if (sumsq) list(ss = sum((X[, 
                                                                seq_len(P2), drop = FALSE] - X.new)^2L))))
}