#' @importFrom utils capture.output
#' @import purrr
#' @import stats
#' @import parallel
#' @import readr
#' @importFrom magrittr %>%
#' @details
#' Bag of little bootstraps works by averaging the results of bootstrapping multiple subsets. BLB is well suited to modern parallel architectures.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param file_list a character vector. file names that containing data.
#' @param path character. the path of the file list.
#' @param m the number of subsamples.
#' @param B the repeat time of bootstrapping for each subsample.
#' @param parallel logical. If TRUE the function will use parallelization.
#' @param num_cl numeric. The number of clusters to be created for parallelization.
#' @param logit logical. If TRUE the function will do logistic regression.
#'
#' @return an object of class "blblm"
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @export
blblm <- function(formula, data = NULL, file_list = NULL, path= '', m = 10, B = 5000, parallel = FALSE, num_cl = 0, logit = FALSE) {

  # check if a file list exists
  if (length(file_list) > 0) {

    # check if parallelization is needed
    if (parallel) {
      cl <- makeCluster(num_cl) ## create clusters
      estimates <- parLapply(cl,
                             file_list,
                             function(x) {lm_each_subsample(formula = formula, data = readr::read_csv(paste(path, x, sep = '')), logit = logit, n = nrow(readr::read_csv(paste(path, x, sep = ''))), B = B)})
      stopCluster(cl)
    } else {
      estimates <- map(
        file_list,
        ~ lm_each_subsample(formula = formula, data = readr::read_csv(paste(path, x, sep = '')), logit = logit, n = nrow(readr::read_csv(paste(path, x, sep = ''))), B = B))
    }

  } else {

    data_list <- split_data(data, m)

    # check if parallelization is needed
    if (parallel) {
      cl <- makeCluster(num_cl)
      estimates <- parLapply(cl,
                             data_list,
                             function(x) {lm_each_subsample(formula = formula, data = x, logit = logit, n = nrow(data), B = B)})
      stopCluster(cl)
    } else {
      estimates <- map(
        data_list,
        ~ lm_each_subsample(formula = formula, data = ., logit = logit, n = nrow(data), B = B))
    }
  }

  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


# split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


# compute the estimates
lm_each_subsample <- function(formula, data, logit, n, B) {
  replicate(B, lm_each_boot(formula, data, logit, n), simplify = FALSE)
}


# compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, logit, n) {
  # generate the weights for data
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, logit, freqs)
}


# estimate the regression estimates based on given the number of repetitions
lm1 <- function(formula, data, logit, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  # check if logistic regression is needed.
  if (logit) {
    fit <- glm(formula, data, family = binomial, weights = freqs)
  } else {
    fit <- lm(formula, data, weights = freqs)
  }
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


# compute the coefficients from fit
blbcoef <- function(object) {
  coef(object)
}


# compute sigma from fit
blbsigma <- function(object) {
  p <- object$rank
  y <- model.extract(object$model, "response")
  e <- fitted(object) - y
  w <- object$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
