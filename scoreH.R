scoreH <- function(y, x, z, type="art", tol=1e-08, maxiter=100, thres=10,
                      x_intercept=TRUE, z_intercept=TRUE) {
  
  h <- switch(type,
              "art" = function(tau) 1+exp(tau),
              "har" = function(tau) exp(tau))
  
  dh <- switch(type,
              "art" = function(tau) exp(tau),
              "har" = function(tau) exp(tau))
  
  if (x_intercept == TRUE) {
    X <- model.matrix(y~x)
  } else {
    X <- as.matrix(x)
  }
  
  if (z_intercept == TRUE) {
    Z <- model.matrix(y~z)
  } else {
    Z <- as.matrix(z)
  }
  n <- length(y)
  p <- ncol(X)
  r <- ncol(Z)
  b <- coef(lm.fit(X, y))
  e <- y - X %*% b
  sigma2 <- (t(e) %*% e)/n
  gamma <- switch(type,
                  "har" = coef(lm.fit(Z, log(e^2))),
                  "art" = coef(lm.fit(Z, e^2/mean(e^2)-1)))
  theta <- rbind(b, gamma)
  
  for (k in 1:maxiter) {
    tau <- Z %*% gamma
    tau[tau > thres] <- thres
    tau[tau < -thres] <- -thres
    w <- h(tau)^(-0.5)
    Xw <- X * drop(w)
    yw <- y * drop(w)
    bk <- coef(modk <- lm.fit(Xw, yw))
    ek <- resid(modk)
    sigma2k <- (t(ek) %*% ek)/n
    sigma2k[sigma2k == 0] <- 1e-12
    qk <- dh(tau)/h(tau)
    Zqk <- Z * drop(qk)
    yqk <- switch(type,
                  "har" = (y-X %*% bk)^2/h(tau) - 1,
                  "art" = (y-X %*% bk)^2/(drop(sigma2k) * h(tau)) - 1)
    yqk[is.na(yqk) | is.infinite(yqk)] <- NA_real_
    gammak <- gamma + coef(lm.fit(Zqk, yqk))
    thetak <- rbind(bk, gammak)
    
    err <- sum((thetak-theta)^2)/sum(theta^2)
    
    if (err < tol) {
      b <- bk
      het <- gammak
      ypred <- X %*% b
      res <- y - ypred
      W <- diag(drop(w))
      weights <- diag(W)
      #Xw <- X * weights
      #yw <- y * weights
      tau <- Z %*% het
      sigma2 <- drop(sigma2k)
      vcov.b <- solve(t(X) %*% W %*% X)
      sigma.hc0 <- diag(drop(res^2))
      vcov.white <- solve(t(X) %*% X) %*% t(X) %*% sigma.hc0 %*% X %*% solve(t(X) %*% X)
      ypred_w <- Xw %*% b
      VarPred_w <- diag(sigma2*(1 + Xw %*% solve(t(X) %*% W %*% X) %*% t(Xw)))
      SEPred_w <- sqrt(VarPred_w)
      VarPred <- diag(sigma2*(drop(h(tau)) + X %*% solve(t(X) %*% W %*% X) %*% t(X)))
      SEPred <- sqrt(VarPred)
      alpha <- 0.01
      tstat <- abs(qt(alpha/2, df=n-p))
      return(list(
        b = bk,
        het = gammak,
        sigma2 = switch(type,
                        "har" = exp(gammak[1]),
                        "art" = sigma2k),
        vcov.b = vcov.b,
        vcov.white = vcov.white,
        weights = weights,
        yw = yw,
        xw = Xw[,-1],
        ypred = ypred,
        ypred_w = ypred_w,
        SEPred = SEPred,
        SEPred_w = SEPred_w,
        tstat = tstat
      ))
    } else {
      gamma <- gammak
      theta <- thetak
    }
  }
  stop(paste("O algoritmo não convergiu. A última iteração teve erro de", err))
}

