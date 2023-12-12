experiment1 <- function(n){
  
  library(knockoff)
  library(tictoc)
  library(distr)
  set.seed(1919)
  # model parameters
  p = 50
  k = 15
  amp = 4.5
  B = 1 # normalization of rows of X
  
  # model sampling
  rho = 0.25
  Sigma <- toeplitz(rho^(0:(p-1)))
  X <- matrix(rnorm(n*p),n)%*%chol(Sigma)
  X <- X/sqrt(rowSums(X^2))*B
  nonzero <- sample(p,k)
  beta <- amp* (1:p %in% nonzero)
  
  # privacy parameters
  epsilon = 0.1
  delta = 1/n
  r = 200 # JLT dimension
  w_sq <- 8*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
  Lap <- DExp(4*B^2/epsilon)
  Z <- r(Lap)(1)
  
  num_iter = 20 # number of iterations to approximate FDR and power
  
  Xbeta = X %*% beta
  fdps = rep(0,num_iter)
  powers = rep(0,num_iter)
  power = function(selected) sum(beta[selected] != 0)/k
  fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
  for (i in 1:num_iter){
    y <- Xbeta + rnorm(n) # resample the noise here for the FDR computation
    
    # lower bounding singular values of [X;y]
    A <- cbind(X,y)
    UDV <- svd(A)
    min_sin <- min(UDV$d)
    
    # Johnson-Lindenstrauss transform
    
    if (min_sin^2 > w_sq+Z+4*B^2*log(1/delta)/epsilon) {
      R = 1/sqrt(r)*matrix(rnorm(r*n),r)
      RX <- R%*%X
      Ry <- R%*%y
      matrix_altered <- FALSE
    } else {
      R = 1/sqrt(r)*matrix(rnorm(r*(n+p+1)),r)
      A_alt <- rbind(A,diag(rep(sqrt(w_sq),p+1)))
      RA <- R%*%A_alt
      RX <- RA[,1:p]
      Ry <- RA[,p+1]
      matrix_altered <- TRUE
    } 
    
    result = knockoff.filter(RX,Ry,fdr=0.2)
    fdps[i] <- fdp(result$selected)
    powers[i] <- power(result$selected)
    cat(rep("-",50),"\n")
    cat(i,"/",num_iter, " iteration for n =",n,
        "\nw^2 =",w_sq,"square of the minimum singular value =", min_sin^2,
        "matrix altered?:", matrix_altered,
        "\nDP Knockoffs select:", result$selected, 
        "\nTrue model:", sort(nonzero),
        "\nFDP =",fdps[i],
        "Power =", powers[i],"\n")
  }
  
  cat(rep("=",50),"\n")
  cat("For n =",n, "FDR=",mean(fdps),"Power=",mean(powers),"\n")
  result_list <- list(fdr = mean(fdps), power = mean(powers))
  return(result_list)
}

fdrs <- rep(0,10)
powers <- rep(0,10)

for (i in 0:10) {
  n = 600000+10000*i
  result <- experiment1(n)
  fdrs[i] = result$fdr
  powers[i] = result$power
}