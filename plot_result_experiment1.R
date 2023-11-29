load("results_experiment1.RData")

sample_size = 50000+10000*c(1:10)
xrange = c(60000,150000)
yrange = c(0,1)

plot(sample_size,fdrs,
     main = "Sample Size vs FDR & Power",
     xlim = xrange,
     ylim = yrange,
     xlab = "Sample Size",
     ylab = "",
     type="b",
     col = "blue")

lines(sample_size, rep(0.2,10),type = "l", col="green",lty=2)

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
