context("Correct labelling of factors")

data("marketing")
mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
mymar2 <- na.omit(mymar)

test_that("Correct labelling of factors",{
  test1 <- rankMANOVA(cbind(Income, Edu) ~ Sex*Language , data = mymar2, iter=10, resampling = "bootstrap")
  mymar2$Sex <- factor(mymar2$Sex, labels = c("M", "F"))
  mymar2$Language <- factor(mymar2$Language, labels = c("English", "Spanish", "Other"))
  test2 <- rankMANOVA(cbind(Income, Edu) ~ Sex*Language , data = mymar2, iter=10, resampling = "bootstrap")
  rownames(test1$Descriptive) <- NULL
  rownames(test2$Descriptive) <- NULL
  expect_equal(test1$Descriptive[, 4:5], test2$Descriptive[, 4:5])
})


test_that("Correct labelling of factors - One-way",{
  test1 <- rankMANOVA(cbind(Income, Edu) ~ Language , data = mymar2, iter=10, resampling = "bootstrap")
  mymar2$Language <- factor(mymar2$Language, labels = c("English", "Spanish", "Other"))
  test2 <- rankMANOVA(cbind(Income, Edu) ~ Language , data = mymar2, iter=10, resampling = "bootstrap")
  rownames(test1$Descriptive) <- NULL
  rownames(test2$Descriptive) <- NULL
  expect_equal(test1$Descriptive[, 3:4], test2$Descriptive[, 3:4])
})
