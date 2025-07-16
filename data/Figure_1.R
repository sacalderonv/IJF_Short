 #Results are very similar to those reported in the paper, but not identical, due
 #to the lack of seeds specified prior to random number generation.

 library(mtarm)
 data(riverflows)
 fit3.5 <- mtar( ~ Bedon + LaPlata | Rainfall, row.names=Date, dist="Laplace",
	               data=riverflows, subset={Date<"2009-04-05"}, 
	               ars=list(p=c(5,5,5)), n.burnin=1000, n.sim=3000, n.thin=2)
 fitGau <- update(fit3.5, dist="Gaussian")
 res <- residuals(fit3.5)
 resGau <- residuals(fitGau)

 dev.new()
 par(mfrow=c(2,2))
 acf(res$by.component[,1],main="(a)",ylim=c(-0.1,0.1),xlim=c(1.5,31))
 acf(res$by.component[,2],main="(b)",ylim=c(-0.1,0.1),xlim=c(1.5,31))
 qqnorm(res$full,pch=20,main="(c)")
 abline(0,1,lty=3)
 qqnorm(resGau$full,pch=20,main="(d)")
 abline(0,1,lty=3)
