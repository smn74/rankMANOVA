context("Post-hoc comparisons")

library(ElemStatLearn)
data("marketing")

test_that("Univariate comparisons - subset of data", {
  mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
  mymar2 <- na.omit(mymar)
  mymar2$Sex <- factor(mymar2$Sex, labels = c("M", "F"))
  mymar2$Language <- factor(mymar2$Language, labels = c("English", "Spanish", "Other"))
  test <- rankMANOVA(cbind(Income, Edu) ~ Sex*Language, data = mymar2, iter=10,
                            resampling = "WildBS", CPU = 1, seed = 987)
  male <- mymar2[mymar2$Sex == "M", ]
  m1 <- rankMANOVA(cbind(Income, Edu) ~ Language, data = male, iter = 1000)

  expect_equal(univariate(m1, data = male, factor = "Language"),
               univariate(test, data = male, factor = "Language"))
})


test_that("Pairwise comparisons one-way", {
  mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
  mymar2 <- na.omit(mymar)
  mymar2$Sex <- factor(mymar2$Sex, labels = c("M", "F"))
  mymar2$Language <- factor(mymar2$Language, labels = c("English", "Spanish", "Other"))
  Mult2<-rankMANOVA(cbind(Income, Edu) ~ Language, data = mymar2, iter = 1000)
  expect_equal(pairwise(Mult2, type="Tukey"), pairwise(Mult2, type="Tukey", factor = "Language"))
})

test_that("Pairwise comparisons one-way II", {
  mymar <- marketing[, c("Sex", "Income", "Edu", "Language")]
  mymar2 <- na.omit(mymar)
  mymar2$Sex <- factor(mymar2$Sex, labels = c("M", "F"))
  mymar2$Language <- factor(mymar2$Language, labels = c("English", "Spanish", "Other"))
  Mult1<-rankMANOVA(cbind(Income, Edu)  ~ Sex, data = mymar2, iter = 1000, seed=123)
  expect_equal(as.character(pairwise(Mult1, type="Tukey")[, "Statistic"]), as.character(Mult1$Test[, "Test statistic"]))
})

> summary(Mult1)
