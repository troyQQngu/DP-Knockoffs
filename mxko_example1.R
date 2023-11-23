library(knockoff)
library(tictoc)

set.seed(1234)
# parameters
n = 1000
p = 1000
k = 60
amp = 4.5

mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n)%*%chol(Sigma)

nonzero = sample(p,k)
beta = amp* (1:p %in% nonzero)/sqrt(n)
y.sample = function(X) X %*% beta +rnorm(n)
y = y.sample(X)

tic()
result = knockoff.filter(X,y)
toc()

print(result)

fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
fdp(result$selected)

