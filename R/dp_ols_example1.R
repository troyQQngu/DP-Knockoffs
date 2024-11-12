## Compute the minimum singular value of matrix A and threshold w for altering 
## matrix A from the DP OLS paper.
##
## In this example we find that using n = 1000 and n = 10000 and define 
## everything else as the DP OLS paper describes, the w^2 is approximately 
## 20,000, the minimum singular value is around 80, which makes the square 1600.
## This is when we set r=0, which means there it algorithm will always alter the
## matrix regardless of our choice of r. And the regression of the altered matrix
## is equivalent to a ridge regression with a penalty constant around 150!

# parameters
set.seed(1234)
n = 100000
p = 3
beta = matrix(c(0.5,-0.25,0),3)

X <- matrix(rnorm(n*p),n)
y.sample <- function(X) X %*% beta + rnorm(n)/sqrt(1-sum(beta^2))
y <- y.sample(X)

# privacy parameters
epsilon = 1
delta = 1e-6
r = 1500# JLT dimension 

# lower bounding singular values of [X;y]
A <- cbind(X,y)
B <- max(sqrt(rowSums(A^2)))
w_sq <- 4*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
Lap <- DExp(4*B^2/epsilon)
Z <- r(Lap)(1) 
A <- cbind(X,y)
UDV <- svd(A)
min_sin <- min(UDV$d)

