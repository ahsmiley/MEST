#' MEST Power Calculator
#'
#' This function helps to determine the sample size needed for MEST. \cr\cr
#' See the article by Smiley, Glazier, and Shoda for details on how to conduct MEST: https://psyarxiv.com/465ha/
#' @param r is the predicted Pearson's r correlation (must be between -1 and 1).
#' @param geni is the Greatest Effect of No Interest in Pearson's r correlation (must be between -1 and 1).
#' @param n is the sample size (# of pairs). Leave as NULL to calculate the necessary sample size with a given power.
#' @param power is the statistical power (must be between 0 and 1). Leave as NULL to calculate the simulated power with a given sample size.
#' @param int.size is the size of the two-tailed confidence interval for MEST (must be between 0 and 1). Default is .9.
#' @param n.sims is the number of simulations used if calculating power from a given n. Defaults to 1000. n.sims is not relevant in calculating necessary sample size. Increase n.sims for greater accuracy.
#' @keywords MEST
#' @export
#' @examples
#' ## Calculate necessary sample size:
#' MEST.power(r = .5, geni = .3, n = NULL, power = .85, int.size = .95)
#'
#' ## Calculate statistical power:
#' MEST.power(r = .5, geni = .3, n = 100, power = NULL, int.size = .95, n.sims = 10000)

MEST.power <- function(r, geni, n = NULL, power = NULL, int.size = .9, n.sims = NULL){
  ## Install and require DescTools package (contains CorCI function)
  if (!require('DescTools')) install.packages('DescTools'); library('DescTools')
  ## Incorrect input
  if(is.null(n) & is.null(power)){
    stop("n and power cannot both be set to NULL")
  } else if(!is.null(n) & !is.null(power)){
    stop("either n or power must be set to NULL")
  } else if(r >= 1 | r <= -1 | geni >= 1 | geni <= -1){
    stop("r and geni must both be between -1 and 1")
  } else if(int.size <= 0 | int.size >= 1){
    stop("int.size must be between 0 and 1")
  } else if(geni == r){
    stop("GENI and r must be different values")
  ## NO n specified - calculate necessary n
  } else if(is.null(n)){
      n <- 5
      if(power <= 0 | power >= 1){
      stop("power must be NULL or between 0 and 1")
      }
      ## Calculate n when predicted correlation is greater than GENI
      if(r > geni){
        lower <- -1
        while (geni > lower) {
          n <- n + 1
          cor.ci <- CorCI(rho = r, n = n, conf.level = int.size, alternative = "two.sided")
          cor.ci2 <- CorCI(rho = unname(cor.ci[2]), n = n, conf.level = power, alternative = "greater")
          lower <- unname(cor.ci2[2])
        }
        cat("One-tailed MEST", "\n",
            "Predicted correlation: ", r, "\n",
            "GENI: ", geni, "\n",
            "Power: ", power*100, "%", "\n",
            "CI: ", int.size*100, "%", "\n\n",
            "***Sample size Needed: ", n, "***",
            sep = "")
      ## Calculate n when predicted correlation is lower than GENI
      } else if(geni > r){
        upper <- 1
        while (geni < upper) {
          n <- n + 1
          cor.ci <- CorCI(rho = r, n = n, conf.level = int.size, alternative = "two.sided")
          cor.ci2 <- CorCI(rho = unname(cor.ci[3]), n = n, conf.level = power, alternative = "less")
          upper <- unname(cor.ci2[3])
        }
        cat("One-tailed MEST", "\n",
            "Predicted correlation: ", r, "\n",
            "GENI: ", geni, "\n",
            "Power: ", power*100, "%", "\n",
            "CI: ", int.size*100, "%", "\n\n",
            "***Sample size Needed: ", n, "***",
            sep = "")
      }
      invisible(list(r = r, geni = geni, n = n, power = power, int.size = int.size))
  } else if(is.null(power)){
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
        sig <- rep(NA, n.sims)
        if(r > geni){
          for(i in 1:n.sims){
            dat <- mvrnorm(n = n, Sigma = lazyCor(X = r, d = 2), mu = c(0,0))
            ## Calculate correlation between the two variables
            correl <- cor.test(dat[,1], dat[,2], conf.level = int.size)
            ## Lower CI of correlation
            sig[i] <- correl$conf.int[1] > geni
          }
          power <- mean(sig)
          cat("One-tailed MEST", "\n",
              "Predicted correlation: ", r, "\n",
              "GENI: ", geni, "\n",
              "Sample size (# of pairs): ", n, "\n",
              "CI: ", int.size*100, "%", "\n\n",
              "***Power: ", 100*mean(sig), "%***", "\n",
              "Note: This was calculated using ", floor(n.sims), " simulations",
              sep = "")
        } else {
          for(i in 1:n.sims){
            dat <- mvrnorm(n = n, Sigma = lazyCor(X = r, d = 2), mu = c(0,0))
            ## Calculate correlation between the two variables
            correl <- cor.test(dat[,1], dat[,2], conf.level = int.size)
            ## Lower CI of correlation
            sig[i] <- correl$conf.int[2] < geni
          }
          power <- mean(sig)
          cat("One-tailed MEST", "\n",
              "Predicted correlation: ", r, "\n",
              "GENI: ", geni, "\n",
              "Sample size (# of pairs): ", n, "\n",
              "CI: ", int.size*100, "%", "\n\n",
              "***Power: ", 100*mean(sig), "%***", "\n",
              "Note: This was calculated using ", floor(n.sims), " simulations",
              sep = "")
          }
        }
    invisible(list(r = r, geni = geni, n = n, power = power, int.size = int.size, n.sims = n.sims))
      }
}


