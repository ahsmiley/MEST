#' EqT Power Calculator
#'
#' This function helps to conduct power analysis for Equivalence Testing with correlational data. To determine the necessary sample size, leave n=NULL and designate the desired power (i.e. power = .8). To calculate the power with a given sample size, leave power=NULL and specify a sample size (i.e. n = 150).
#' @param r is the predicted Pearson's r correlation (must be between -1 and 1).
#' @param low.sesoi is the lower smallest effect in Pearson's r correlation (must be between -1 and 1 and must be less than r).
#' @param up.sesoi is the upper smallest effect in Pearson's r correlation (must be between -1 and 1 and must be greater than r).
#' @param n is the sample size (# of pairs). Leave as NULL to calculate the necessary sample size with a given power.
#' @param power is the desired power (must be between 0 and 1).
#' @param int.size is the size of the two-tailed confidence interval for EqT (must be between 0 and 1). Default is .9.
#' @param n.sims is the number of simulations used if calculating power from a given n. Defaults to 1000. n.sims is not relevant in calculating necessary sample size. Increase n.sims for greater accuracy.
#' @keywords Equivalence
#' @export
#' @examples
#' EQT.power(r = 0, low.sesoi = -.2, up.sesoi = .2, n = NULL, power = NULL, int.size = .95)

EQT.power <- function(r, low.sesoi, up.sesoi, n = NULL, power = NULL, int.size = .9, n.sims = NULL){
  if (!require('DescTools')) install.packages('DescTools'); library('DescTools')
  ## Incorrect input
  if(is.null(n) & is.null(power)){
    stop("n and power cannot both be set to NULL")
  } else if(!is.null(n) & !is.null(power)){
    stop("either n or power must be set to NULL")
  } else if(r >= 1 | r <= -1 | low.sesoi >= 1 | low.sesoi <= -1 | up.sesoi >= 1 | low.sesoi <= -1){
    stop("r, low.sesoi, and up. sesoi must be between -1 and 1")
  } else if(int.size <= 0 | int.size >= 1){
    stop("int.size must be between 0 and 1")
  } else if(low.sesoi >= r){
    stop("low.sesoi must be less than r")
  } else if(up.sesoi <= r){
    stop("up.sesoi must be greater than r")
    ## Equivalence Testing power calculation
  } else if(is.null(n) & low.sesoi < r & r < up.sesoi){
    n <- 5
    upper <- 1
    lower <- -1
    if(power <= 0 | power >= 1){
    stop("power must be between 0 and 1")
    }
    while (upper > up.sesoi | lower < low.sesoi) {
      n <- n + 1
      cor.ci <- CorCI(rho = r, n = n, conf.level = int.size, alternative = "two.sided")
      CI.low <- CorCI(rho = unname(cor.ci[2]), n = n, conf.level = power, alternative = "greater")
      lower <- unname(CI.low[2])
      CI.up <- CorCI(rho = unname(cor.ci[3]), n = n, conf.level = power, alternative = "less")
      upper <- unname(CI.up[3])
    }
    cat("Equivalence Test", "\n",
        "Predicted correlation: ", r, "\n",
        "Upper SESOI: ", up.sesoi, "\n",
        "Lower SESOI: ", low.sesoi, "\n",
        "Power: ", power*100, "%", "\n",
        "CI: ", int.size*100, "%", "\n\n",
        "***Sample size Needed (# of pairs): ", n, "***",
        sep = "")
    invisible(list(r = r, up.sesoi = up.sesoi, low.sesoi = low.sesoi, n = n, power = power, int.size = int.size))
  } else if(is.null(power) & low.sesoi < r & r < up.sesoi) {
      if (!require('rockchalk')) install.packages('rockchalk'); library('rockchalk')
      if(is.null(n.sims)){
        n.sims <- 1000
      }
      if(n.sims < 1){
        stop("n.sims must be at least 1")
      } else {
        if(n.sims < 1000){
          warning("n.sims is less than 1000. Increase n.sims for greater accuracy.")
        }
          for(i in 1:n.sims){
            dat <- mvrnorm(n = n, Sigma = lazyCor(X = r, d = 2), mu = c(0,0))
            ## Calculate correlation between the two variables
            correl <- cor.test(dat[,1], dat[,2], conf.level = int.size)
            ## Lower CI of correlation
            sig[i] <- (correl$conf.int[1] > low.sesoi) & (correl$conf.int[2] < up.sesoi)
          }
          power <- mean(sig)
          cat("Equivalence Test", "\n",
              "Predicted correlation: ", r, "\n",
              "Upper SESOI: ", up.sesoi, "\n",
              "Lower SESOI: ", low.sesoi, "\n",
              "Sample size (# of pairs): ", n, "\n",
              "CI: ", int.size*100, "%", "\n\n",
              "***Power: ", 100*power, "%***", "\n",
              "Note: This was calculated using ", floor(n.sims), " simulations",
              sep = "")
          invisible(list(r = r, up.sesoi = up.sesoi, low.sesoi = low.sesoi, n = n, power = power, int.size = int.size, n.sims = n.sims))
  }
}
}

