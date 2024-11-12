
experiment4 <- function(n,r) {
  library(knockoff)
  library(tictoc)
  library(distr)
  # set seed
  set.seed(1919)
  # model parameters
  p = 50
  k = 15
  amp = 4.5
  B = 1 # normalization of rows of X
  
  # model sampling 
  # notice here we use independent features
  X <- matrix(rnorm(n*p),n)
  X <- X/sqrt(max(rowSums(X^2)))*B
  nonzero <- sample(p,k)
  beta <- amp* (1:p %in% nonzero)
  
  # privacy parameters
  epsilon = 0.2
  delta = 2*p/n
  w_sq <- 8*B^2/epsilon*(sqrt(2*r*log(8/delta))+2*log(8/delta))
  Lap <- DExp(4*B^2/epsilon)
  Z <- r(Lap)(1)
  
  num_iter = 20# number of iterations to approximate FDR and power
  
  Xbeta = X %*% beta
  fdps = rep(0,num_iter)
  powers = rep(0,num_iter)
  power = function(selected) sum(beta[selected] != 0)/k
  fdp = function(selected) sum(beta[selected] == 0)/max(1,length(selected))
  for (i in 1:num_iter){
    y <- Xbeta + rnorm(n) # resample the noise here for the FDR computation
    
    # lower bounding singular values of [X;y]
    A <- cbind(X,y)
    UDV <- svd(X)
    min_sin <- min(UDV$d)
    s <- min_sin^2+Z-4*B^2*log(1/delta)/epsilon
    # Johnson-Lindenstrauss transform
    
    w_adj <- sqrt(max(0,w_sq-s))
    
    R = 1/sqrt(r)*matrix(rnorm(r*(n+p+1)),r)
    A_alt <- rbind(A,diag(rep(w_adj,p+1)))
    RA <- R%*%A_alt
    RX <- RA[,1:p]
    Ry <- RA[,p+1]
    matrix_altered <- TRUE

    
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

n_iter = 13

fdrs <- rep(0,n_iter)
powers <- rep(0,n_iter)

exp <- seq(3,6,length.out = n_iter)
r = 400

for (i in 1:n_iter) {
  n <- as.integer(10^exp[i])
  result <- experiment4(n,r)
  fdrs[i] = result$fdr
  powers[i] = result$power
}

sample_size <- exp
xrange = c(exp[1],exp[n_iter])
yrange = c(0,1)

png("experiment_4_r400_plot.png")

plot(sample_size,fdrs,
     main = "Sample Size vs FDR & Power",
     xlim = xrange,
     ylim = yrange,
     xlab = "log10(n)",
     ylab = "",
     type="b",
     col = "blue")

lines(sample_size, rep(0.2,n_iter),type = "l", col="green",lty=2)

lines(sample_size,powers,
      xlab="",
      ylab="",
      col="red",
      type="b",
      lty=1,
      xlim=xrange,
      ylim=yrange)

legend("topleft",legend = c("FDR","FDR=0.2","Power"),
       col=c("blue","green","red"),
       lty =c(1,2,1),
       cex=0.8)
dev.off()

