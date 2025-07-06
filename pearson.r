library(energy)
library(glmnet)

set.seed(123)
n <- 100
p <- 500
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- rnorm(n)

# Compute Pearson correlation 
compute_pearson <- function(X, y) {
  pearson_values <- cor(X, y, method = 'pearson')
  names(pearson_values) <- paste0("X", 1:p)#name columns for model later
  return(sort(pearson_values, decreasing = TRUE))
}

pearson_values <- compute_pearson(X, y)
print(pearson_values)