experiment1 <- function(n){
  
  library(knockoff)
  library(tictoc)
  library(distr)
  
  set.seed(1234)
  # model parameters
  p = 200
  k = 20
  amp = 5
  B = 1 # normalization of rows of X
  
  # model sampling
  rho = 0.25
  Sigma <- toeplitz(rho^(0:(p-1)))
  X <- matrix(rnorm(n*p),n)%*%chol(Sigma)
  X <- X/sqrt(rowSums(X^2))*B
  nonzero <- sample(p,k)
  beta <- amp* (1:p %in% nonzero)
  
  # privacy parameters
  epsilon = 0.5
  delta = 1e-4
  r = 200 # JLT dimension
  w_sq <- 4*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
  Lap <- DExp(4*B^2/epsilon)
  Z <- r(Lap)(1)
  
  
  num_iter = 50 # number of iterations to approximate FDR and power
  
  Xbeta = X %*% beta
  fdps = rep(0,num_iter)
  powers = rep(0,num_iter)
  power = function(selected) sum(beta[selected] != 0)/k
  fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
  for (i in 1:num_iter){
    set.seed(Sys.time())
    y <- Xbeta +rnorm(n) # resample the noise here for the FDR computation
    
    # lower bounding singular values of [X;y]
    A <- cbind(X,y)
    UDV <- svd(A)
    min_sin <- min(UDV$d)
    
    # Johnson-Lindenstrauss transform
    R = 1/sqrt(r)*matrix(rnorm(r*n),r)
    if (min_sin^2 > w_sq+Z+4*B^2*log(1/delta)/epsilon) {
      RX <- R%*%X
      Ry <- R%*%y
      matrix_altered <- FALSE
    } else {
      RA <- rbind(R%*%A,diag(rep(sqrt(w_sq),p+1)))
      RX <- RA[,1:p]
      Ry <- RA[,p+1]
      matrix_altered <- TRUE
    } 
    
    result = knockoff.filter(RX,Ry,fdr=0.2)
    fdps[i] <- fdp(result$selected)
    powers[i] <- power(result$selected)
    cat(rep("-",40),"\n")
    cat(i,"/",num_iter, " iteration for n =",n,"\nDP Knockoffs selects:", result$selected, "\nTrue model:", sort(nonzero),"\nFDP =",fdps[i],"Power =", powers[i],"\n")
  }
  
  cat("For n =",n, "FDR=",mean(fdps),"Power=",mean(powers))
  result_list <- list(fdr = mean(fdps), power = mean(powers))
  return(result_list)
}

fdrs <- rep(0,15)
powers <- rep(0,15)

for (i in 1:10) {
  n = 50000+10000*i
  result <- experiment1(n)
  fdrs[i] = result$fdp
  powers[i] = result$power
}