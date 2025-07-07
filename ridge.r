# Installing necessary packages 
library(ggplot2)
library(dplyr)
library(energy)
library(glmnet)


# New Dataset for gene expression 
bcTCGA <- readRDS("~/Documents/Research/R Files/bcTCGA.rds")

summary(bcTCGA)

# Create Random High dimensional dataset  
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


# Call function compute_dcor and save the results  
get_mse <- function(fit) {
  min(fit$cvm)
}

fit_ridge <- function(x_sub, y) {
  x_mat <- as.matrix(x_sub)
  fit <- cv.glmnet(x_mat, y, alpha = 0)
  #best_lamb <- fit$lambda.min
 # min_mse <- min(fit$cvm)
  #cat("Mini CV MSE:", min_mse, "\n")
  
  return(fit)
}

subset_by_dcor <- function(X, dcor_values, threshold) {
  keep_vars <- names(dcor_values[dcor_values > threshold])
  if (length(keep_vars) == 0) {
    stop("No predictors meet the distance correlation threshold.")
  }
  return(X[, keep_vars, drop = FALSE])
}

# Fit Several Linear Models for various threshold values of distance correlation
res <- compute_dcor(X, y)
thresholds <- seq(0.1, 1.0, by = 0.1)
mses <- numeric(length(thresholds))

# This function will use get_mse() and fit_ridge() to compute all of the dcor values and give MSE 
compute_lm_dcor <- function(X, y, res, thresholds){
  # for linear models 
  models <- list()
  # compare the res with the different thresholds 
  for (i in seq_along(thresholds)) {
    t <- thresholds[i]
    cat("\nThreshold:", t, "\n")
    
    # 1) find the variables that fit within the thresholds   
    x_sub <- tryCatch(subset_by_dcor(X, res, t), error = function(e) {
      cat("Threshold", t, "was skipped -", e$message, "\n")
      mses[i] <- NA 
      models[[paste0("Thresh", t)]] <- NULL 
      return(NULL)
    }
    )
    if (is.null(x_sub)) next
    # 2) Fit the linear model with the subset of dcor variables 
    fit <- fit_ridge (x_sub, y)
    mse <- get_mse(fit)
    mses[i] <- mse
    models[[paste0("thresh_", t)]] <- fit
    # Mean Square Error 
    mses[i] <- get_mse(fit)
  }
  # Trying to figure out the best lambda to use for this dataset 
  mse_results <- data.frame(threshold = thresholds, MSE = mses)
  
  # Return the values of mse_results 
  return(mse_results)
}

L <- compute_lm_dcor(X, y, res, thresholds)
print(L)

