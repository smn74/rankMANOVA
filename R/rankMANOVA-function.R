#' Rank-based Tests for Multivariate Data in Semi-Parametric Factorial Designs
#'
#' The rankMANOVA function calculates the ANOVA-type
#' statistic (ATS) with bootstrap p-values.
#'
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in
#'   \code{formula}. Data must be in long format.
#' @param subject The column name of the subjects in the data.
#' @param iter The number of iterations used for calculating the resampled
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS"
#'   (parametric bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights). The Wild Bootstrap is calculated for all test statistics.
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no
#'   reproducible seed is set.
#' @param nested.levels.unique A logical specifying whether the levels of the nested factor(s)
#'   are labeled uniquely or not. Default is FALSE, i.e., the levels of the nested
#'   factor are the same for each level of the main factor.
#'
#'
#' @details
#'
#' @section
#'
#' @return A \code{MANOVA} object containing the following components:
#'   \item{Descriptive}{Some descriptive statistics of the data for all factor
#'   level combinations. Displayed are the number of individuals per factor
#'   level combination and the vector of means (one column per dimension).}
#'   \item{Covariance}{The estimated covariance matrix.}
#'   \item{WTS}{The value of the WTS along with degrees of freedom of the
#'   central chi-square distribution and p-value.}
#'   \item{ATS}{The value of the
#'   ATS, degrees of freedom of the central F distribution and the corresponding
#'   p-value.}
#'   \item{MATS}{The value of the MATS.}
#'   \item{resampling}{p-values for the test statistic based on the
#'   chosen resampling approach.}
#'
#'
#' @examples
#'
#' @seealso \code{\link{MANOVA.RM}}
#'
#' @references Konietschke, F., Bathke, A. C., Harrar, S. W. and Pauly, M. (2015).
#'   Parametric and nonparametric bootstrap methods for general MANOVA. Journal
#'   of Multivariate Analysis, 140, 291-301.
#'
#'   Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal data
#'   in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
#'
#'    Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and Hoeller, Y. (2016).
#'    Using EEG, SPECT, and Multivariate Resampling Methods
#'    to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
#'
#'   Friedrich, S., Konietschke, F., Pauly, M. (2016). GFD - An
#'   R-package for the Analysis of General Factorial Designs. Accepted for publication in
#'   Journal of Statistical Software.
#'
#'   Friedrich, S., and Pauly, M. (2017). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. arXiv preprint arXiv:1704.03731.
#'
#'
#' @importFrom graphics axis legend par plot title abline points
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom quantile
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster parSapply detectCores
#' @importFrom ellipse ellipse
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
  split1 <- strsplit(outcome_names, "(", fixed = TRUE)[[1]][-1]
  split2 <- strsplit(split1, ")", fixed = TRUE)[[1]]
  split3 <- strsplit(split2, ",")[[1]]

  EF <- rownames(nr_hypo)[-1]  # names of influencing factors
  nf <- length(EF)
  names(dat) <- c("response", EF)
  #no. dimensions
  d <- ncol(dat$response)
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
    names(statistic_out) <- c("Test statistic", "p-value")
    # calculate results
    results <- rankbs(Y, n, hypo, d, iter, alpha, CPU, seed, resampling)
    statistic_out <- results$statistic

  } else {
    # crossed design
    dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
    fac.groups <- do.call(list, dat2[, 2:(nf+1)])
    n.groups <- prod(fl)

    Y <- split(dat2, fac.groups)
    n <- sapply(Y, nrow)
    hypo_matrices <- HC_MANOVA(fl, perm_names, fac_names, d)[[1]]

    # ---------------------- error detection ------------------------------------

    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not implemented!")
    }
    # only 3-way nested designs are possible
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are
           not implemented!")
    }
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
  output$relEffect <- p_out
  output$Test <- statistic_out
  class(output) <- "rankMANOVA"
  return(output)
}
