#' Univariate post-hoc comparisons in nonparametric multivariate factorial designs
#'
#' @param object A \code{rankMANOVA} object.
#' @param factor The factor for which univariate comparisons are desired. Must be one of the factors
#' used in the main analysis, of course.
#' @param data The data set to be used for the analysis. If none is specified, the data used for fitting
#' \code{object} is re-used.
#' @param ... Not used yet.
#'
#' @details The univariate function computes p-values for univariate comparisons.
#' Details on the tests can be found in Dobler, Friedrich and Pauly (2019).
#' Note that due to the formulation of our effect size vectors, these tests can be performed by applying
#' the closed testing principle, i.e., no alpha-corretion is needed.
#'
#' @return p-values for the univariate post-hoc comparisons of the chosen factor.
#'
#' @references  Dobler, D., Friedrich, S., and Pauly, M. (2019).
#'              Nonparametric MANOVA in meaningful effects. Annals of the Institute of Statistical Mathematics.
#'
#' @export
univariate <- function(object, factor, data, ...){
  if(!(factor %in% object$other$fac_names)){
    stop("The specified factor is not found!")
  }

  test1 <- hasArg(data)
  if(!test1){
    data <- object$input$data
  }

  d <- object$other$dim
  # outcome variables
  outcome <- object$other$outcomes
  alpha <- object$input$alpha
  input <- object$input

  desc <- lapply(1:d, function(i){
    form <- as.formula(paste(outcome[i], "~", factor))
    fit <- rankMANOVA(form, data = data, iter = input$iter, alpha = alpha, dec = input$dec,
                      resampling = input$resampling, CPU = input$CPU)
    out <- list(de = fit$Descriptive[,ncol(fit$Descriptive)], pval = fit$Test)
    return(out)
  })

  de <- matrix(unlist(lapply(desc, function(x) x$de)), length(levels(data[, factor])), d)
  pval <- matrix(unlist(lapply(desc, function(x) x$pval)), d, 2, byrow = TRUE)

  colnames(de) <- outcome
  rownames(de) <- levels(data[, factor])

  rownames(pval) <- outcome
  colnames(pval) <- c("Test statistic", "p-value")

  # avoid printing zeros
  pval[pval[, "p-value"] == 0, "p-value"] <- "<0.001"

  out <- list(effects = de, pval = pval)
  return(out)
}


#' Pairwise post-hoc comparisons in nonparametric multivariate factorial designs
#'
#' @param object A \code{rankMANOVA} object.
#' @param type The type of the pairwise comparison must be specified here.
#' Calculation is based on the contrMat function in package multcomp, see the corresponding help page
#' for details on the types of contrasts available.
#' @param base an interger specifying which group is considered the baseline group
#' for Dunnett contrasts, see \code{\link[multcomp]{contrMat}}.
#' @param factor The factor for which pairwise comparisons are desired. Must be one of the factors
#' used in the main analysis, of course.
#' @param uni Logical: should univariate comparisons also be computed? Default is FALSE.
#' @param ... Not used yet.
#'
#' @details The pairwise function computes p-values for pairwise comparisons as provided by
#'  \code{\link[multcomp]{contrMat}}. For the nonparametric comparisons,
#'  only Tukey's all-pairwise and Dunnett's many-to-one comparisons are used.
#' Note that due to the formulation of our effect size vectors, these tests can be performed by applying
#' the closed testing principle, i.e., no alpha-corretion is needed, see Dobler, Friedrich and Pauly (2019)
#' for details.
#'
#' @return p-values for the multivariate and (if desired) univariate pairwise comparisons of the chosen factor.
#'
#' @references  Dobler, D., Friedrich, S., and Pauly, M. (2019).
#'              Nonparametric MANOVA in meaningful effects. Annals of the Institute of Statistical Mathematics.
#'
#' @seealso \code{\link[multcomp]{contrMat}}
#'
#' @importFrom multcomp contrMat
#' @importFrom stats as.formula
#'
#' @export

pairwise <- function(object, type = NULL, base = 1,
                     factor, uni = FALSE, ...){

  if(is.null(type)){
    stop("Please specify the type of pairwise comparison, see the multcomp-package for
           details.")
  }

  if(!(type %in% c("Dunnett", "Tukey"))){
    stop("Only Dunnett's many-to-one and Tukey's all pairwise comparisons can be used.")
  }

  if(!(factor %in% object$other$fac_names)){
    stop("The requested factor is not part of the model.")
  }

  if(grepl(":", factor)){
    stop("Pairwise comparisons cannot be computed for interaction effects.")
  }

  input <- object$input
  n <- object$Descriptive$n
  lev <- subset(object$Descriptive, select = 1:n)
  lev <- lev[-ncol(lev)]
  d <- object$other$dim
  nf <- object$other$nf

  # one-way
  if(nf == 1){
    names(n) <- lev[, 1]
  } else {
    names(n) <- sort(do.call(paste, c(lev, sep = " ")))
  }
  M <- contrMat(n, type = type, base)
  contmat <- M %x% diag(d)

  i <- 1
  cm <- list()
  while(i < nrow(contmat)){
    cm[[i]] <- contmat[i:(i+d-1), ]
    i <- i+d
  }
  cm <- cm[-which(sapply(cm, is.null))]

  pval <- t(sapply(cm, function(x){
    H <- t(x)%*%MASS::ginv(x%*%t(x))%*%x
    results <- rankbs(object$other$Y2, n, H, d, iter = input$iter, alpha = input$alpha,
                      resampling = input$resampling, CPU = input$CPU, seed = input$seed)
    statistic_out <- round(results$statistic, input$dec)
  }))

  rownames(pval) <- rownames(M)
  colnames(pval) <- c("Statistic", "multivariate p-value")

  # add univariate comparisons if desired
  if(uni){
    outcome <- object$other$outcomes

    pval2 <- lapply(1:d, function(i){
      form <- as.formula(paste(outcome[i], "~", factor))

      # necessary steps to get Y2
      dat <- model.frame(form, object$input$data)
      nr_hypo <- attr(terms(form), "factors")
      EF <- rownames(nr_hypo)[-1]  # names of influencing factors
      names(dat) <- c("response", EF)

      dat2 <- dat[order(dat[, 2]), ]
      fac.groups <- dat2[, 2]
      Y <- split(dat2, fac.groups)
      nuni <- sapply(Y, nrow)
      Y2 <- lapply(Y, function(x) x$response)
      Y2 <- lapply(Y2, function(x) as.matrix(x))

      M <- contrMat(n, type = type, base)

      out <- t(apply(M, 1, function(x) {
      # x is a row vector, hence t() must be switched in the computation of H
      H <- x%*%MASS::ginv(t(x)%*%x)%*%t(x)
      results <- rankbs(Y2, nuni, H, d = 1, iter = input$iter, alpha = input$alpha,
                        resampling = input$resampling, CPU = input$CPU, seed = input$seed)
      statistic_out <- round(results$statistic, input$dec)
      }))

      colnames(out) <- c("Statistic", outcome[i])
      return(out)
    })

    pval2 <- do.call(cbind, pval2)
    pval <- cbind(pval, pval2)

    pval <- as.data.frame(pval)
    pval <-  pval[, names(pval) != "Statistic"]
  }
  # avoid printing zeros
  pval[cbind(row(pval)[which(pval == 0)], col(pval)[which(pval == 0)])] <- "<0.001"

  return(pval)
}




