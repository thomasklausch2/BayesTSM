# Transforms correlation matrix to variance-covariance matrix using SDs
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}