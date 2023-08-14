fwdH <- function(y, x, z, type="art", bsb=NULL, init=NULL) {
  h <- switch(type,
              "art" = function(tau) 1+exp(tau),
              "har" = function(tau) exp(tau))
  
  dh <- switch(type,
               "art" = function(tau) exp(tau),
               "har" = function(tau) exp(tau))
  
  X <- model.matrix(y~x)
  Z <- model.matrix(y~z)
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Z)
  
  if(is.null(init)) {
    if (n < 40) {
      init <- p + 1
    } else {
      init <- min(c(3*p+1, floor(0.5*(n+p+1))))
    }
  }
  
  if (is.null(bsb)) {
    Ra <- TRUE
    nwhile <- 0
    while (Ra & nwhile < 100) {
      bsb <- sample(1:n, p)
      Xb <- X[bsb, ]
      Zb <- Z[bsb, ]
      Ra <- !(qr(Xb)$rank == p)
      nwhile <- nwhile + 1
    }
    if (nwhile == 100) {
      print("Incapaz de selecionar um subconjunto inicial com matriz de especificação de posto completo.")
    }
    yb <- y[bsb]
  } else {
    Xb <- X[bsb, ]
    Zb <- Z[bsb, ]
    yb <- y[bsb]
  }
  
  m0 <- unit <- length(bsb)
  
  if (init < p) {
    init <- p
  } else if (init < m0) {
    init <- m0
  } else if (init >= n) {
    init <- n - 1
  }
  
  R <- cbind(1:n, matrix(0, n, 1))
  seq100 <- 100*seq(1, ceiling(n/100), by=1)
  bgls <- cbind(init:n, matrix(NA_real_, n-init+1, p+1))
  blast <- c()
  S2 <- cbind(init:n, matrix(NA_real_, n-init+1, 3))
  mdr <- cbind(init:(n-1), matrix(0, n-init, 2))
  msr <- cbind(init:n, matrix(0, n-init+1, 2))
  coo <- cbind((init+1):n, matrix(NA_real_, n-init, 6))
  res <- matrix(NA_real_, n, n-init+1)
  bb <- matrix(NA_real_, n, n-init+1)
  LEV <- matrix(NA_real_, n, n-init+1)
  Un <- cbind((init+1):n, matrix(NA_real_, n-init, 10))
  WEI <- matrix(NA_real_, n, n-init+1)
  HET <- cbind(init:n, matrix(NA_real_, n-init+1, r))
  
  ncl <- setdiff(1:n, bsb)
  hhh <- 1
  if (qr(Xb)$rank != p) {
    print("O conjunto inicial fornecido não tem matriz de especificação de posto completo.")
  } else {
    
    for (m in m0:n) {
      
      if (n > 200) {
        if (length(intersect(m, seq100)) == 1) {
          cat("m = ", m, "\n")
        }
      }
      
      NoRankProblem <- (qr(Xb)$rank == p)
        
      if (NoRankProblem == TRUE) {
        Hart <- scoreH(yb, Xb, Zb, type="art", maxiter = 15000,
                      x_intercept=FALSE, z_intercept=FALSE)
        
        gamma <- Hart$het
        tau <- Z %*% gamma
        w <- h(tau)^(-0.5)
        Xw <- X * drop(w)
        
        if (m >= init) {
          WEI[, m-init+1] <- w
        }
        
        yw <- y * drop(w)
        Xb <- Xw[bsb,]
        yb <- yw[bsb]
        
        b <- Hart$b
        res_bsb <- yb - Xb %*% b
        blast <- b
      } else {
        print("A matriz de especificação não tem posto completo.")
        blast <- b
      }
      
      if (hhh == 1) {
        e <- yw - Xw %*% b
      } else {
        e <- y - X %*% b
      }
      
      if (m >= init) {
        bb[bsb, m-init+1] <- bsb
        
        if (NoRankProblem == TRUE) {
          bgls[m-init+1, 2:(p+1)] <- b
          HET[m-init+1, 2:(r+1)] <- gamma
          
          mAm <- t(Xb) %*% Xb
          mmX <- solve(mAm)
          dmmX <- diag(mmX)
          hi <- rowSums((Xw %*% mmX) * Xw)
          
          LEV[bsb, m-init+1] <- hi[bsb]
        }
      }
      
      if (m > p) {
        if (NoRankProblem == TRUE) {
          Sb <- (t(res_bsb) %*% res_bsb)/(m-p)
        }
      } else {
        Sb <- 0
      }
      
      if (m >= init) {
        res[, m-init+1] <- e
        
        if (NoRankProblem == TRUE) {
          S2[m-init+1, 2] <- Sb
          
          if (m < n) {
            a <- qnorm(0.5*(1+m/n))
            corr <- 1 - 2*(n/m)*a*dnorm(a)
          } else {
            corr <- 1
          }
          
          Sbrescaled <- Sb/corr
          S2[m-init+1, 4] <- Sbrescaled
          msrsel <- sort(abs(res_bsb)/sqrt(Sb*hi[bsb]))
          msr[m-init+1, 2] <- msrsel[m]
          
          S2[m-init+1, 3] <- 1 - var(res_bsb)/var(yb)
        }
      }
      
      R[, 2] <- e^2
      
      if (m > init) {
        if (NoRankProblem == TRUE) {
          bib <- bgls[m-init+1, 2:(p+1)] - bgls[m-init, 2:(p+1)]
          
          if (S2[m-init+1, 2] > 0) {
            #coo[m-init, 3:length(unit)+2] <- 1/(1-hi[unit]) * sqrt()
          }
          
          # TO-DO: FSRHeda.m 917-922
          
          
        }
      }
      
      if (m < n) {
        if (m >= init) {
          #browser()
          if (NoRankProblem == TRUE) {
            if (hhh == 1) {
              ord <- cbind(R[, 2]/(1+hi), e)
            } else {
              ord <- cbind(R[, 2]/(w+hi), e)
            }

            
            selmdr <- min(ord[ncl, 1]) 
            
            if (S2[m-init+1, 2] == 0) {
              warning("O valor de S2 é zero, mdr é NaN")
            } else {
              mdr[m-init+1, 2] <- sqrt(selmdr/Hart$sigma2)
            }
            
            selmdr <- ord[order(ord[, 1]) , ]
            mdr[m-init+1, 3] <- sign(selmdr[m+1, 2])*sqrt(selmdr[m+1, 1]/Hart$sigma2)
          }
        }
      }
      
    
    oldbsb <- bsb
      
    ord <- order(R[, 2])
    
    bsb <- ord[1:(m+1)]
    
    Xb <- X[bsb, ]
    yb <- y[bsb]
    Zb <- Z[bsb, ]
    
    if (m < (n-1)) {
      ncl <- ord[(m+2):n]
    }
    
   # if (m >= init) {
   #    unit <- setdiff(bsb, oldbsb)
   #   if (length(unit) <= 10) {
   #      Un[m-init+1, 2:(length(unit)+1)] <- unit
   #   } else {
   #      Un[m-init+1, 2:ncol(Un)] <- unit[1:10]
   #   }
   #  }
      
    } # end for
  }
  
  return(list(
    HET = HET,
    S2 = S2,
    WEI = WEI,
    mdr = mdr,
    bb = bb
  ))
}


