test_that("fit works", {

  x <- rnorm(1000)
  y <- rnorm(1000)
  data <- as.data.frame(cbind(x,y))
  names(data) <- c('x', 'y')

  fit <- blblm(y~x, data, logit = FALSE, m = 3, B = 100)
  expect_s3_class(fit, 'blblm')

  co <- coef(fit)
  expect_equal(length(co), 2)

  ci <- confint(fit)
  expect_equal(dim(ci), c(1,2))

  sigma_ci <- sigma(fit, confidence = TRUE)
  expect_equal(length(sigma_ci), 3)

  pred <- predict(fit, data.frame(x = c(2.5, 3)))
  expect_equal(length(pred), 2)

  pred_ci <- predict(fit, data.frame(x = c(2.5, 3)), confidence = TRUE)
  expect_equal(dim(pred_ci), c(2,3))


  z <- as.integer(rnorm(1000) > 1)
  data1 <- as.data.frame(cbind(x,z))
  names(data1) <- c('x', 'z')

  fit1 <- blblm(z~x, data1, logit = TRUE, m = 3, B = 100)
  expect_s3_class(fit1, 'blblm')

  co1 <- coef(fit1)
  expect_equal(length(co1), 2)

  ci1 <- confint(fit1)
  expect_equal(dim(ci1), c(1,2))

  sigma_ci1 <- sigma(fit1, confidence = TRUE)
  expect_equal(length(sigma_ci1), 3)

  pred1 <- predict(fit1, data.frame(x = c(2.5, 3)), type="response")
  expect_equal(length(pred1), 2)

  pred_ci1 <- predict(fit1, data.frame(x = c(2.5, 3)), , type="response", confidence = TRUE)
  expect_equal(dim(pred_ci1), c(2,3))

})
