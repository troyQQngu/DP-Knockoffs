# plot convergence in r with different n

sample_size <- exp
xrange = c(exp[1],exp[r_iter])
yrange = c(0,1)

png("experiment_5_fdrs.png")

plot(sample_size,fdrs_ne4_5,
     main = "Projected Sample Size vs FDR",
     xlim = xrange,
     ylim = yrange,
     xlab = "log10(r)",
     ylab = "",
     type="b",
     col = "blue")

lines(sample_size,fdrs_ne5,
      xlab="",
      ylab="",
      col="red",
      type="b",
      lty=1,
      xlim=xrange,
      ylim=yrange)

lines(sample_size, rep(0.2,r_iter),type = "l", col="green",lty=2)

legend("topleft",legend = c("n=10^4.5","n=10^5","FDR=0.2"),
       col=c("blue","green","red"),
       lty =c(1,2,1),
       cex=0.8)
dev.off()


png("experiment_5_powers.png")

plot(sample_size,powers_ne4_5,
     main = "Projected Sample Size vs Power",
     xlim = xrange,
     ylim = yrange,
     xlab = "log10(r)",
     ylab = "",
     type="b",
     col = "blue")

lines(sample_size,powers_ne5,
      xlab="",
      ylab="",
      col="red",
      type="b",
      lty=1,
      xlim=xrange,
      ylim=yrange)


legend("topleft",legend = c("n=10^4.5","n=10^5"),
       col=c("blue","red","green"),
       lty =c(1,1,2),
       cex=0.8)
dev.off()
