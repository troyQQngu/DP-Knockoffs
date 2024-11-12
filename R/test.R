# this script is used for small tests
library(knockoff)
library(tictoc)
library(distr)
library(ggplot2)

# set seed
set.seed(1919)
# model parameters
n = as.integer(10^5)
p = 50
k = 15
amp = 4.5
B = 1 # normalization of rows of X
r = 3000

# model sampling 
# independent features
X <- matrix(rnorm(n*p),n)
# X <- X/sqrt(max(rowSums(X^2)))*B # X is normalized in row norm <= 1
nonzero <- sample(p,k)
beta <- amp* (1:p %in% nonzero)

# privacy parameters
epsilon = 0.2
delta = 2*p/n

Lap <- DExp(4*B^2/epsilon)
Z <- r(Lap)(1)

num_iter = 20 # number of iterations to approximate FDR and power

Xbeta = X %*% beta
fdps = rep(0,num_iter)
powers = rep(0,num_iter)
power = function(selected) sum(beta[selected] != 0)/k
fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))

for (i in 1:num_iter){
  E <- rnorm(n)
  y <- Xbeta + E # resample the noise here for the FDR computation
  y <- y/sqrt(k*amp^2+1) # normalize y so that it's about the same size as X
  
  # lower bounding singular values of [X;y]
  A <- cbind(X,y)
  A <- A/sqrt(max(rowSums(A^2)))*B # normalize A so that the row norm of A is bounded by B
  # B <- max(sqrt(rowSums(A^2)))
  w_sq <- B^2*(1+(4/epsilon+2/log(8/delta))*(sqrt(2*r*log(8/delta))+log(8/delta)))
  UDV <- svd(t(A)%*%A)
  min_sin_sq <- min(UDV$d)
  
  s <- min_sin_sq+Z-4*B^2*log(1/delta)/epsilon # private release of minimum singular value
  
  # Johnson-Lindenstrauss transform
  w_adj <- sqrt(max(0,w_sq-s))
  R = 1/sqrt(r)*matrix(rnorm(r*(n+p+1)),r)
  A_alt <- rbind(A,diag(rep(w_adj,p+1)))
  RA <- R %*% A_alt
  RX <- RA[,1:p]
  Ry <- RA[,p+1]
  
  result = knockoff.filter(RX,Ry,fdr=0.2)
  fdps[i] <- fdp(result$selected)
  powers[i] <- power(result$selected)
  cat(rep("-",50),"\n")
  cat(i,"/",num_iter, " iteration for n =",n,"r=",r,
      "\nw^2 =",w_sq,"min_sin^2 =", min_sin_sq,"adjusted w^2=",w_adj^2,
      "\nDP Knockoffs select:", result$selected, 
      "\nTrue model:", sort(nonzero),
      "\nFDP =",fdps[i],
      "Power =", powers[i],"\n")
}

cat(rep("=",50),"\n")
cat("For n =",n, "FDR=",mean(fdps),"Power=",mean(powers),"\n")
result_list <- list(fdr = mean(fdps), power = mean(powers))

