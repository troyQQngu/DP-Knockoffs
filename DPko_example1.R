library(knockoff)
library(tictoc)
library(distr)

set.seed(1234)
# model parameters
n = 500
p = 500
k = 30
amp = 100

# model setup
mu <- rep(0,p)
rho = 0.25
Sigma <- toeplitz(rho^(0:(p-1)))
X <- matrix(rnorm(n*p),n)%*%chol(Sigma)

nonzero <- sample(p,k)
beta <- amp* (1:p %in% nonzero)/sqrt(n)
y.sample <- function(X) X %*% beta +1/4*rnorm(n)
y <- y.sample(X)

# privacy parameters
epsilon = 0.5
delta = 1e-1
r = 400 # JLT dimension 

# lower bounding singular values of [X;y]
B <- max(sqrt(rowSums(X^2)))
w_sq <- 8*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
Lap <- DExp(4*B^2/epsilon)
Z <- r(Lap)(1) 
A <- cbind(X,y)
UDV <- svd(A)
min_sin <- min(UDV$d)

# Johnson-Lindenstrauss transform
R = 1/sqrt(r)*matrix(rnorm(r*n),r)
if (TRUE) { #(min_sin^2 > w_sq+Z+4*B^2*log(1/delta)/epsilon) {
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
result = knockoff.filter(RX,Ry)
toc()

print(result)

fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
fdp(result$selected)

