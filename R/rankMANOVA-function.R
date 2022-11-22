#' Rank-based Tests for Multivariate Data in Nonparametric Factorial Designs
#'
#' The rankMANOVA function calculates an ANOVA-type
#' statistic (ATS) with (wild) bootstrap p-values for nonparametric factorial designs with
#' multivariate data.
#'
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variables and the right hand side contains the factor
#'   variables of interest. Data must be
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
#' @param nested.levels.unique A logical specifying whether the levels of the nested factor(s)
#'   are labeled uniquely or not. Default is FALSE, i.e., the levels of the nested
#'   factor are the same for each level of the main factor. For an example and more explanations
#'   see the GFD package and the corresponding vignette.
#' @param dec Number of decimals the results should be rounded to. Default is 3.
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
#'  data("marketing")
#'  mymar <- marketing[, c("Sex", "Income", "Edu")]
#'  mymar2 <- na.omit(mymar)
#'  test <- rankMANOVA(cbind(Income, Edu) ~ Sex, data = mymar2, iter=100,
#'   resampling = "WildBS", CPU = 1)
#'  summary(test)
#'
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
                       iter = 10000, alpha = 0.05, CPU, dec = 3,
                       seed, resampling = "bootstrap", nested.levels.unique = FALSE){

  if (!(resampling %in% c("bootstrap", "WildBS"))){
    stop("Resampling must be one of 'bootstrap' and 'WildBS'!")
  }

  output <- list()

  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }

  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }

  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha, resampling = resampling, dec = dec,
                     CPU = CPU, seed = seed)

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
    fl[aa] <- nlevels(as.factor(as.character(dat[, (aa + 1)])))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(as.character(dat[, (jj + 1)])))
  }
  lev_names <- expand.grid(levels)

  # number of hypotheses
  tmp <- 0
  for (i in 1:nf) {
    tmp <- c(tmp, choose(nf, i))
    nh <- sum(tmp)
  }
  # correct formula?
  if (length(fac_names) != nf && length(fac_names) != nh){
    stop("Something is wrong with the formula. Please specify all or no interactions in crossed designs.")
  }

  # mixture of nested and crossed designs is not possible
  if (length(fac_names) != nf && 2 %in% nr_hypo) {
    stop("A model involving both nested and crossed factors is
           not implemented!")
  }


  if (nf == 1) {
    # one-way layout
    nest <- FALSE
    dat2 <- dat[order(dat[, 2]), ]
    fac.groups <- dat2[, 2]
    hypo_matrices <- list((diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(d))

    #------------------------ end one-way layout -------------------------------------------------#

  } else {
    dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
    fac.groups <- do.call(list, dat2[, 2:(nf+1)])

  }
  Y <- split(dat2, fac.groups)
  n <- sapply(Y, nrow)

  nested <- grepl(":", formula)
  nested2 <- grepl("%in%", formula)

  if (sum(nested) > 0 || sum(nested2) > 0) {
    # nested
    nest <- TRUE

    # if nested factor is named uniquely
    if (nested.levels.unique){
      # delete factorcombinations which don't exist
      n <- n[n != 0]
      # create correct level combinations
      blev <- list()
      lev_names <- list()
      for (ii in 1:length(levels[[1]])) {
        blev[[ii]] <- droplevels(as.factor(as.character(dat[, 3][dat[, 2] == levels[[1]][ii]])))
        lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
      }
      if (nf == 2) {
        lev_names <- as.factor(unlist(lev_names))
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev)
      } else {
        lev_names <- lapply(lev_names, rep,
                            length(levels[[3]]) / length(levels[[2]]))
        lev_names <- lapply(lev_names, sort)
        lev_names <- as.factor(unlist(lev_names))
        blev <- lapply(blev, rep, length(levels[[3]]) / length(levels[[2]]))
        blev <- lapply(blev, sort)
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev, as.factor(levels[[3]]))
      }
      # correct for wrong counting of nested factors
      if (nf == 2) {
        fl[2] <- fl[2] / fl[1]
      } else if (nf == 3) {
        fl[3] <- fl[3] / fl[2]
        fl[2] <- fl[2] / fl[1]
      }
    }
    hypo_matrices <- HN_MANOVA(fl, d)
  } else {
    # crossed
    nest <- FALSE

    ## adapting formula argument, if interaction term missing
    if (nrow(perm_names) != nh) {
      #stop("For crossed designs, an interaction term must be specified in the formula.")
      form2 <- as.formula(paste(outcome_names, "~", paste(fac_names, collapse = "*")))
      perm_names2 <- t(attr(terms(form2), "factors")[-1, ])
      fac_names2 <- attr(terms(form2), "term.labels")
      hyps <- HC_MANOVA(fl, perm_names2, fac_names2, d, nh)
      hypo_matrices <- hyps[[1]]
      fac_names2 <- hyps[[2]]
      # choose only relevant entries of the hypo matrices
      indices <- grep(":", fac_names2, invert = T)
      hypo_matrices <- lapply(indices, function(x) hypo_matrices[[x]])

    } else if(nf !=1){
      hyps <- HC_MANOVA(fl, perm_names, fac_names, d, nh)
      hypo_matrices <- hyps[[1]]
      fac_names <- hyps[[2]]
    }
  }

  # correcting for "empty" combinations (if no interaction specified)
  n.groups <- prod(fl)
  if(nf != 1 & length(Y) != n.groups){
    index <- NULL
    for(i in 1:length(Y)){
      if(nrow(Y[[i]]) == 0){
        index <- c(index, i)
      }
    }
    Y <- Y[-index]
  }
  Y2 <- lapply(Y, function(x) x$response)
  if (d==1){
    Y2 <- lapply(Y2, function(x) as.matrix(x))
  }

  # ---------------------- error detection ------------------------------------

  # only 3-way nested designs are possible
  if (sum(nested) > 0 && nf >= 4) {
    stop("Four- and higher way nested designs are
           not implemented!")
  }
  # no factor combinations with less than 2 observations
  if (0 %in% n || 1 %in% n) {
    stop("There is at least one factor-level combination
           with less than 2 observations!")
  }

  #--------------------------------------------------------------------------#


  statistic_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2) # Test statistic, p-value
  rownames(statistic_out) <- fac_names
  colnames(statistic_out) <- c("Test statistic", "p-value")
  # calculate results
  for (i in 1:length(hypo_matrices)) {
    results <- rankbs(Y2, n, hypo_matrices[[i]], d, iter, alpha, CPU, seed, resampling)
    statistic_out[i, ] <- round(results$statistic, dec)
  }

  p_out <- round(results$p, dec)
  descriptive <- cbind(unique(lev_names), n, p_out)
  colnames(descriptive) <- c(EF, "n", split3)

  # other information needed, eg, for post-hoc tests-------------------#
  other <- list(dim=d, nf = nf, fl = fl, outcomes = split3, fac_names = fac_names,
                Y2 = Y2, quant = results$quant)


  # Output ------------------------------------------------------
  output$input <- input_list
  output$Descriptive <- descriptive
  output$Test <- statistic_out
  output$other <- other
  class(output) <- "rankMANOVA"
  return(output)
}
