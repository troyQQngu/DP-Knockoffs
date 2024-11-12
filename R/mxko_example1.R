library(knockoff)
library(tictoc)

# parameters
n = 200
p = 50
k = 15
amp = 4.5
B = 1 # normalization of rows of X

rho = 0
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n)%*%chol(Sigma)
X <- X/sqrt(max(rowSums(X^2)))*B


nonzero = sample(p,k)
beta = amp* (1:p %in% nonzero)
y.sample = function(X) X %*% beta +rnorm(n)
y = y.sample(X)

tic()
result = knockoff.filter(X,y)
toc()

print(result)

fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
fdp(result$selected)
power = function(selected) sum(beta[selected]!=0)/k
power(result$selected)

