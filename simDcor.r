library(energy)
library(glmnet)
#is legit just the profs code 
set.seed(123)
n <- 100
p <- 500
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)#name columns for model later
y <- rnorm(n)

#Compute distance correlation
compute_dcor <- function(X, y) {
  dcor_values <- apply(X, 2, function(xj) dcor(xj, y))
  names(dcor_values) <- colnames(X)
  return(sort(dcor_values, decreasing = TRUE))
}

# Makes subset by dCor threshold passed in
subset_by_dcor <- function(X, dcor_results, threshold) {
  keep_vars <- names(dcor_results[dcor_results > threshold])
  if (length(keep_vars) == 0) {
    stop("No predictors meet the threshold.")
  }
  return(X[, keep_vars, drop = FALSE])
}

# 4) Function for LASSO 
fit_lasso <- function(X_sub, y) {
  X_mat <- as.matrix(X_sub)
  fit <- cv.glmnet(X_mat, y, alpha = 1)
  return(fit)
}

#Gets minimum MSE 
get_mse <- function(fit) {
  min(fit$cvm)
}

#loops through threshholds in .1 increases 
thresholds <- seq(0.1, 1.0, by = 0.1)
models <- list()
mses <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  cat("\nThreshold:", t, "\n")
  tryCatch(
    {
      X_sub <- subset_by_dcor(X, dcor_results, threshold = t)
      cat("Subset dims:", dim(X_sub), "\n") #displays how much fits in the threshold requirment 
      
      fit <- fit_lasso(X_sub, y)
      models[[paste0("thresh_", t)]] <- fit
      mses[i] <- get_mse(fit)
      
      # I still cant get it to plot it how it wants will fix later... maybe idek
      # plot(fit, main = paste("LASSO CV Plot - Threshold:", t))
    },
    error = function(e) {
      cat("No predictors meet threshold.\n")
      mses[i] <- NA
    }
  )
}

#prints the accuracy for each threshold 
mse_results <- data.frame(threshold = thresholds, MSE = mses)
print(mse_results)
