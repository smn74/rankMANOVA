context("Overall functioning")

library(ElemStatLearn)
data("marketing")

test_that("One-way layout", {
  mymar <- marketing[, c("Sex", "Income", "Edu")]
  mymar2 <- na.omit(mymar)
  test.oneway <- rankMANOVA(cbind(Income, Edu) ~ Sex, data = mymar2, iter=10,
                            resampling = "WildBS", CPU = 1, seed = 987)
  expect_equal_to_reference(summary(test.oneway), file = "oneway.rds")
})


test_that("Two-way interaction",{
  mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
  mymar2 <- na.omit(mymar)
  test1 <- rankMANOVA(cbind(Income, Edu) ~ Sex*Language , data = mymar2, iter=10, resampling = "bootstrap", seed = 123)
  expect_equal_to_reference(summary(test1), file = "interaction.rds")
})

test_that("Two-way no interaction", {
  mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
  mymar2 <- na.omit(mymar)
  test.add <- rankMANOVA(cbind(Income, Edu) ~ Sex + Language, data = mymar2, iter=10, resampling = "bootstrap", seed = 456)
  expect_equal_to_reference(summary(test.add), file = "no-interaction.rds")
})

