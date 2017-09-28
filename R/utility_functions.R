#' @export
print.rankMANOVA <- function(x, ...) {
  a <- x$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Test:", "\n", sep = "")
  print(x$Test)
}

#' @export
summary.rankMANOVA <- function (object, ...) {
  a <- object$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Test:", "\n", sep = "")
  print(object$Test)
}
