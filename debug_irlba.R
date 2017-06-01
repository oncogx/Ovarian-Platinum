# Made this into its own file because you cant make breakpoints in RMDs 
install.packages("irlba")
library(irlba)
sparse_pca <- function(x, n_pcs, mu=NULL, s=NULL, center_scale=TRUE) {
  if (is.null(mu) && center_scale) mu <- colMeans(x)
  if (is.null(s) && center_scale) s <- apply(x, 2, sd, na.rm=TRUE)
  
  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba::irlba(x, n_pcs, center=mu, scale=s) # Error occurs on this line
  } else {
    svd_res <- irlba::irlba(x, n_pcs)
  }
  
  # compute explained variance
  n <- dim(x)[1]
  variance_sum <- sum(apply(x,2,var,na.rm=TRUE)/(s^2)) # sample variance sum
  var_pcs <- svd_res$d^2/(n-1)/variance_sum
  
  return(list(x=svd_res$u %*% diag(svd_res$d), rotation=svd_res$v, sdev=svd_res$d/sqrt(n-1),
              tot_var=variance_sum, var_pcs=var_pcs))
}

pcaSRS <- sparse_pca(t(exprs(gbmSRS_log)), 10)
