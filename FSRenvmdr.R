FSRenvmdr <- function(n, p, init, prob) {
  m0 <- init
  lp <- length(prob)
  probf <- 1-prob
  m <- matrix(m0:(n-1), n-m0, 1)
  mn <- nrow(m)
  mm <- matrix(rep(m, lp), n-m0, lp)
  pp <- matrix(rep(probf, each=mn), mn, lp)
  quantf <- qf(pp, 2*(n-mm), 2*(mm+1))
  quantb <- (mm+1)/(mm+1+(n-mm)*quantf)
  min_sca <- abs(qt(0.5*(1+quantb), mm-p))
  a <- qnorm(0.5*(1+mm/n))
  corr <- 1-2*(n/mm)*a*dnorm(a)
  MDRenv <- cbind(m, min_sca/sqrt(corr))
  
  return(MDRenv)
}

