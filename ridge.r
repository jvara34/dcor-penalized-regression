# Installing necessary packages 
library(ggplot2)
library(dplyr)
library(energy)
library(glmnet)
library(stats)
library(AICcmodavg)


# New Dataset for gene expression 
bcTCGA <- readRDS("~/Documents/Research/R Files/bcTCGA.rds")

summary(bcTCGA)

# Synthetic dataset
# Only working with the synthetic dataset will work on the real dataset later 
set.seed(123)
n <- 100
p <- 500
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)
y <- rnorm(n)

compute_dcor <- function(X, y) {
  dcor_values <- apply(X, 2, function(xj) dcor(xj, y))
  names(dcor_values) <- paste0("X", 1:ncol(X))
  return(sort(dcor_values, decreasing = TRUE))
}

subset_by_dcor <- function(X, dcor_values, threshold) {
  keep_vars <- names(dcor_values[dcor_values > threshold])
  if (length(keep_vars) == 0) {
    stop("No predictors meet the distance correlation threshold.")
  }
  return(X[, keep_vars, drop = FALSE])
}

# Ridge regression with AIC and BIC
fit_ridge_with_aic_bic <- function(x_sub, y) {
  x_mat <- as.matrix(x_sub)
  ridge_fit <- cv.glmnet(x_mat, y, alpha = 0)
  lambda_min <- ridge_fit$lambda.min
  final_fit <- glmnet(x_mat, y, alpha = 0, lambda = lambda_min)
  
  y_pred <- predict(final_fit, newx = x_mat, s = lambda_min)
  residuals <- y - y_pred
  rss <- sum(residuals^2)
  
  df <- final_fit$df
  n <- length(y)
  
  # Compute AIC and BIC
  aic <- n * log(rss / n) + 2 * df
  bic <- n * log(rss / n) + log(n) * df
  # Returning both AIC and BIC 
  return(list(fit = final_fit, AIC = aic, BIC = bic, MSE = mean(residuals^2)))
}

# Thresholds and results
res <- compute_dcor(X, y)
thresholds <- seq(0.1, 1.0, by = 0.1)

ridge_results <- data.frame(threshold = thresholds, MSE = NA, AIC = NA, BIC = NA)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  cat("\nThreshold:", t, "\n")
  
  x_sub <- tryCatch(
    subset_by_dcor(X, res, t),
    error = function(e) {
      cat("Threshold", t, "skipped -", e$message, "\n")
      return(NULL)
    }
  )
  if (is.null(x_sub)) next
  
  result <- fit_ridge_with_aic_bic(x_sub, y)
  ridge_results$MSE[i] <- result$MSE
  ridge_results$AIC[i] <- result$AIC
  ridge_results$BIC[i] <- result$BIC
}
# Printing the Threshold, MSE, AIC, and the BIC results 
print(ridge_results)

