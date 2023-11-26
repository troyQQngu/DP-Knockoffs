library(knockoff)
library(tictoc)
library(distr)

set.seed(1234)
# model parameters
n = 20000
p = 200
k = 20
amp = 5
B = 1 # normalization of rows of X

# model setup
rho = 0.25
Sigma <- toeplitz(rho^(0:(p-1)))
X <- matrix(rnorm(n*p),n)%*%chol(Sigma)
X <- X/sqrt(rowSums(X^2))*B
nonzero <- sample(p,k)
beta <- amp* (1:p %in% nonzero)
y.sample <- function(X) X %*% beta +rnorm(n)
y <- y.sample(X)

# privacy parameters
epsilon = 2.5
delta = 1e-4
r = 200 # JLT dimension


# lower bounding singular values of [X;y]
A <- cbind(X,y)
w_sq <- 4*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
Lap <- DExp(4*B^2/epsilon)
Z <- r(Lap)(1) 
UDV <- svd(A)
min_sin <- min(UDV$d)

# Johnson-Lindenstrauss transform
R = 1/sqrt(r)*matrix(rnorm(r*n),r)
if (min_sin^2 > w_sq+Z+4*B^2*log(1/delta)/epsilon) {
  RX <- R%*%X
  Ry <- R%*%y
  matrix_altered <- FALSE
} else {
  RX <- rbind(R%*%X,diag(rep(sqrt(w_sq),p)))
  Ry <- rbind(R%*%y,matrix(0,nrow=p,ncol=1))
  matrix_altered <- TRUE
} 

set.seed(Sys.time())
tic()
result = knockoff.filter(RX,Ry,fdr=0.2)
toc()

print(result)

fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
fdp(result$selected)
power = function(selected) length(selected)/k
power(result$selected)

