#' Rank-based Tests for Multivariate Data in Nonparametric Factorial Designs
#'
#' The rankMANOVA function calculates an ANOVA-type
#' statistic (ATS) with (wild) bootstrap p-values for nonparametric factorial designs with
#' multivariate data.
#'
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variables and the right hand side contains the factor
#'   variables of interest. An interaction term must be specified. Data must be
#'   provided in wide format.
#' @param data A data.frame containing the variables in
#'   \code{formula}.
#' @param iter The number of iterations used for calculating the resampled
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "bootstrap"
#'   (sample-specific bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights). The default is "WildBS".
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no
#'   reproducible seed is set.
#'
#' @details Implemented is an ANOVA-type test statistic for testing hypotheses formulated in Mann-Whitney-type
#'  effects in nonparametric factorial designs. Statistical inference is based on a wild or a sample-specific
#'  bootstrap approach. The unweighted treatment effects considered do not depend on sample sizes and allow for
#'  transitive ordering. The package thus provides an extension of the univariate \code{\link[rankFD]{rankFD}}
#'  package to multivariate data.
#'
#' @return A \code{rankMANOVA} object containing the following components:
#' \item{Descriptive}{Some descriptive statistics of the data for all factor
#'   level combinations. Displayed are the number of individuals per factor
#'   level combination and the unweighted treatment effects for each dimension.}
#'   \item{Test}{The test statistic(s) and p-value(s) based on the
#'   chosen bootstrap approach.}
#'
#'@section NOTE: The number of bootstrap iterations has been set to 100 in the examples due
#' to runtime restrictions on CRAN. Usually it is recommended to use at least 1000 iterations.
#'
#' @examples
#'  library(ElemStatLearn)
#'  data("marketing")
#'  mymar <- marketing[, c("Sex", "Income", "Edu")]
#'  mymar2 <- na.omit(mymar)
#'  test <- rankMANOVA(cbind(Income, Edu) ~ Sex, data = mymar2, iter=100,
#'   resampling = "WildBS", CPU = 1)
#'  summary(test)
#'
#' @seealso \code{\link[rankFD]{rankFD}}
#'
#' @references
#' Dobler, D., Friedrich, S., and Pauly, M. (2017). Nonparametric MANOVA in Mann-Whitney effects.
#'
#' @importFrom graphics axis legend par plot title abline points
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom quantile
#' @importFrom methods hasArg
#' @importFrom parallel makeCluster parSapply detectCores
#'
#' @export

rankMANOVA <- function(formula, data,
                       iter = 10000, alpha = 0.05, CPU,
                       seed, resampling = "WildBS"){

  if (!(resampling %in% c("bootstrap", "WildBS"))){
    stop("Resampling must be one of 'bootstrap' and 'WildBS'!")
  }

  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha)

  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }

  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }

  dat <- model.frame(formula, data)
  nr_hypo <- attr(terms(formula), "factors")
  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  fac_names <- colnames(nr_hypo)

  outcome_names <- rownames(nr_hypo)[1]  # names of outcome variables
  # extract names of outcome variables
  if (grepl("cbind", outcome_names)){
  split1 <- strsplit(outcome_names, "(", fixed = TRUE)[[1]][-1]
  split2 <- strsplit(split1, ")", fixed = TRUE)[[1]]
  split3 <- strsplit(split2, ",")[[1]]
  } else {
    split3 <- outcome_names
  }

  EF <- rownames(nr_hypo)[-1]  # names of influencing factors
  nf <- length(EF)
  names(dat) <- c("response", EF)
  #no. dimensions
  d <- ncol(dat$response)
  if (!is.numeric(d)){
    d <- 1
  }
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, (aa + 1)]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, (jj + 1)]))
  }
  lev_names <- expand.grid(levels)

  if (nf == 1) {
    # one-way layout
    dat2 <- dat[order(dat[, 2]), ]
    fac.groups <- dat2[, 2]
    n.groups <- prod(fl)
    Y <- split(dat2, fac.groups)
    n <- sapply(Y, nrow)
    hypo <- (diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(d)
    statistic_out <- rep(NA, 2) # Test statistic, p-value
    # calculate results
    results <- rankbs(Y, n, hypo, d, iter, alpha, CPU, seed, resampling)
    statistic_out <- results$statistic
    names(statistic_out) <- c("Test statistic", "p-value")
  } else {
    # crossed design
    dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
    fac.groups <- do.call(list, dat2[, 2:(nf+1)])
    n.groups <- prod(fl)

    Y <- split(dat2, fac.groups)
    n <- sapply(Y, nrow)
    hypo_matrices <- HC_MANOVA(fl, perm_names, fac_names, d)[[1]]

    # ---------------------- error detection ------------------------------------
    # no factor combinations with less than 2 observations
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    #--------------------------------------------------------------------------

    statistic_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2) # Test statistic, p-value
    rownames(statistic_out) <- fac_names
    colnames(statistic_out) <- c("Test statistic", "p-value")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
      results <- rankbs(Y, n, hypo_matrices[[i]], d, iter, alpha, CPU, seed, resampling)
      statistic_out[i, ] <- results$statistic
    }
  }

  p_out <- results$p
  descriptive <- cbind(lev_names, n, p_out)
  colnames(descriptive) <- c(EF, "n", split3)

  # Output ------------------------------------------------------
  output <- list()
  output$input <- input_list
  output$Descriptive <- descriptive
  output$Test <- statistic_out
  class(output) <- "rankMANOVA"
  return(output)
}
