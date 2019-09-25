#' Display rankMANOVA object
#'
#' Returns a short summary of the results (test statistic with p-values)
#'
#' @param x A rankMANOVA object
#' @param ... Additional parameters (currently not used)
#'
#' @export
print.rankMANOVA <- function(x, ...) {
  object <- x
  a <- object$input
  # avoid printing zeros
  Test <- object$Test
  Test[Test[, "p-value"] == 0, "p-value"] <- "<0.001"
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Test:", "\n", sep = "")
  print(Test)
}


#' Summarizing a rankMANOVA object
#'
#' Returns a summary of the results including sample sizes and relative treatment effects for all groups as well
#' as the test statistic with resampling-based p-values
#'
#' @param object A rankMANOVA object
#' @param ... Additional parameters (currently not used)
#'
#' @export
summary.rankMANOVA <- function (object, ...) {
  a <- object$input
  Test <- object$Test
  Test[Test[, "p-value"] == 0, "p-value"] <- "<0.001"
bs <- switch(object$input$resampling,
             bootstrap = "nonparametric bootstrap",
             WildBS = "wild bootstrap")

  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Test: p-values based on the ", bs, "\n", sep = "")
  print(Test)
}
