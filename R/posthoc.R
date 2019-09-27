#' Univariate post-hoc comparisons in nonparametric multivariate factorial designs
#'
#' @param object A \code{rankMANOVA} object.
#' @param factor The factor for which univariate comparisons are desired. Must be one of the factors
#' used in the main analysis, of course. Defaults to all factors in the model (without interactions).
#' @param data The data set to be used for the analysis. If none is specified, the data used for fitting
#' \code{object} is re-used.
#' @param ... Not used yet.
#'
#' @details The univariate function computes p-values for univariate comparisons.
#' If no factor is specified, all factors in the model are used. In this case, the unweighted effects
#' don't change compared to the multivariate model and are thus not returned, only p-values and test
#' statistics for the univariate tests.
#' Details on the tests can be found in Dobler, Friedrich and Pauly (2019).
#' Note that due to the formulation of our effect size vectors, these tests can be performed by applying
#' the closed testing principle, i.e., no alpha-correction is needed.
#'
#' NOTE: If an interaction is significant in the main model, the data set should be split accordingly
#' for the post-hoc analyses. Thus, the univariate comparisons are not computed for interaction effects.
#'
#' @return p-values for the univariate post-hoc comparisons of the chosen factor.
#'
#' @references  Dobler, D., Friedrich, S., and Pauly, M. (2019).
#'              Nonparametric MANOVA in meaningful effects. Annals of the Institute of Statistical Mathematics.
#'
#' @export
univariate <- function(object, factor = NULL, data, ...){

  if(!is.null(factor) && !(factor %in% object$other$fac_names)){
    stop("The specified factor is not found!")
  }

  nf <- object$other$nf
  fac.spec <- TRUE

  if(is.null(factor)){
    factor <- object$other$fac_names
    if(nf != 1){
      factor <- paste(object$other$fac_names[!grepl(":", object$other$fac_names)], collapse = "+")
    }
    fac.spec <- FALSE
  }

  if(grepl(":", factor)){
    stop("Univariate comparisons cannot be computed for interaction effects.")
  }

  test1 <- hasArg(data)
  if(!test1){
    data <- object$input$data
  }

  if((nf == 1 || fac.spec ) && !is.factor(data[, factor])){
    data[, factor] <- as.factor(as.character(data[, factor]))
    #warning(paste("The variable", factor, "was transformed into a factor. Please check for correct labeling!"))
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

  if(fac.spec){
    de <- matrix(unlist(lapply(desc, function(x) x$de)), length(levels(data[, factor])), d)
    pval <- matrix(unlist(lapply(desc, function(x) x$pval)), d, 2, byrow = TRUE)

    colnames(de) <- outcome
    rownames(de) <- levels(data[, factor])

    rownames(pval) <- outcome
    colnames(pval) <- c("Test statistic", "p-value")
    # avoid printing zeros
    pval[pval[, "p-value"] == 0, "p-value"] <- "<0.001"

    out <- list(effects = de, pval = pval)
  } else {
    pval <- matrix(unlist(lapply(desc, function(x) x$pval)), d, 2)
    rownames(pval) <- outcome
    colnames(pval) <- c("Test statistic", "p-value")
    # avoid printing zeros
    pval[pval[, "p-value"] == 0, "p-value"] <- "<0.001"

    out <- list(pval = pval)
  }



  return(out)
}


#' Pairwise post-hoc comparisons in nonparametric multivariate factorial designs
#'
#' @param object A \code{rankMANOVA} object. The pairwise comparisons can only be
#' performed for one factor at a time, so this model must be a one-way model.
#' @param type The type of the pairwise comparison must be specified here.
#' Calculation is based on the contrMat function in package multcomp, see the corresponding help page
#' for details on the types of contrasts available.
#' @param base An integer specifying which group is considered the baseline group
#' for Dunnett contrasts, see \code{\link[multcomp]{contrMat}}.
#' @param factor If desired, the specific factor for which the calculations shall be computed.
#' Defaults to all factors present in the model (without interactions, though).
#' @param uni Logical: should univariate comparisons also be computed? Default is FALSE.
#' @param ... Not used yet.
#'
#' @details The pairwise function computes p-values for pairwise comparisons as provided by
#'  \code{\link[multcomp]{contrMat}}. For the nonparametric comparisons,
#'  only Tukey's all-pairwise and Dunnett's many-to-one comparisons are used.
#' Note that due to the formulation of our effect size vectors, these tests can be performed by applying
#' the closed testing principle, i.e., no alpha-correction is needed, see Dobler, Friedrich and Pauly (2019)
#' for details.
#'
#' NOTE: If an interaction is significant in the main model, the data set should be split accordingly
#' for the post-hoc analyses. Thus, the pairwise comparisons are not computed for interaction effects.
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
                     factor = NULL, uni = FALSE, ...){

  if(is.null(type)){
    stop("Please specify the type of pairwise comparison, see the multcomp-package for
           details.")
  }

  if(!(type %in% c("Dunnett", "Tukey"))){
    stop("Only Dunnett's many-to-one and Tukey's all pairwise comparisons can be used.")
  }

  input <- object$input
  d <- object$other$dim
  nf <- object$other$nf

  # for a specific factor only
  if(! is.null(factor) || nf == 1){

    if(is.null(factor) && nf == 1){
      factor <- object$other$fac_names
    }

    if(!(factor %in% object$other$fac_names)){
      stop("The requested factor is not part of the model.")
    }

    if(grepl(":", factor)){
      stop("Pairwise comparisons cannot be computed for interaction effects.")
    }

    # number of levels of the relevant factor:
    n <- object$Descriptive$n
    n.fac <- 1:object$other$fl[which(object$other$fac_names == factor)]
    lev <- unique(object$Descriptive[, factor])
    names(n.fac) <- lev
    M <- contrMat(n.fac, type = type, base)
    # additionally inflate the contrast matrix by the number of other factor levels
    if(nf != 1){
      p.n.fl <- prod(object$other$fl[-which(object$other$fac_names == factor)])
    } else {
      p.n.fl <- 1
      #contmat <- M %x% diag(d)
    }
  } else {
    # for all factors in the model
    n <- object$Descriptive$n
    lev <- subset(object$Descriptive, select = 1:n)
    lev <- lev[-ncol(lev)]
    names(n) <- do.call(paste, c(lev, sep = " "))
    M <- contrMat(n, type = type, base)
    p.n.fl <- 1
  }

  contmat <- M %x% diag(d) %x% diag(p.n.fl)
  i <- 1
  cm <- list()
  while(i < nrow(contmat)){
    cm[[i]] <- contmat[i:(i+d*p.n.fl-1), ]
    i <- i+d*p.n.fl
  }
  if(length(cm)!=1){
    cm <- cm[-which(sapply(cm, is.null))]
  }

  pval <- t(sapply(cm, function(x){
    H <- t(x)%*%MASS::ginv(x%*%t(x))%*%x
    results <- rankbs(object$other$Y2, n, H, d, iter = input$iter, alpha = input$alpha,
                      resampling = input$resampling, CPU = input$CPU, seed = input$seed)
    statistic_out <- round(results$statistic, input$dec)
  }))



  # contmat <- M %x% diag(d)
  #
  # pval <- t(apply(contmat, 1, function(x){
  #   # x is a row vector, hence t() must be switched in the computation of H
  #   H <- x%*%MASS::ginv(t(x)%*%x)%*%t(x)
  #   results <- rankbs(object$other$Y2, n, H, d, iter = input$iter, alpha = input$alpha,
  #                     resampling = input$resampling, CPU = input$CPU, seed = input$seed)
  #   statistic_out <- round(results$statistic, input$dec)
  # }))


  rownames(pval) <- rownames(M)
  colnames(pval) <- c("Statistic", "multivariate p-value")

  # add univariate comparisons if desired
  if(uni){
    outcome <- object$other$outcomes

    if(is.null(factor)){
      factor <- object$other$fac_names
      if(nf != 1){
        factor <- paste(object$other$fac_names[!grepl(":", object$other$fac_names)], collapse = "+")
      }
    }

    pval2 <- lapply(1:d, function(i){
      form <- as.formula(paste(outcome[i], "~", factor))

      # necessary steps to get Y2
      dat <- model.frame(form, object$input$data)
      nr_hypo <- attr(terms(form), "factors")
      EF <- rownames(nr_hypo)[-1]  # names of influencing factors
      names(dat) <- c("response", EF)
      if (!grepl("+", factor, fixed = TRUE)) {
        dat2 <- dat[order(dat[, 2]), ]
        fac.groups <- dat2[, 2]
      } else {
        dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
        fac.groups <- do.call(list, dat2[, 2:(nf+1)])
      }
      Y <- split(dat2, fac.groups)
      nuni <- sapply(Y, nrow)
      Y2 <- lapply(Y, function(x) x$response)
      Y2 <- lapply(Y2, function(x) as.matrix(x))

      M <- contrMat(nuni, type = type, base)

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




