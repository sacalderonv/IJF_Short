 #Results are very similar to those reported in the paper, but not identical, due
 #to the lack of seeds specified prior to random number generation.

 library(mtarm)
 data(riverflows)
 fit3.5 <- mtar( ~ Bedon + LaPlata | Rainfall, row.names=Date, dist="Laplace",
	               data=riverflows, subset={Date<"2009-04-05"}, 
	               ars=list(p=c(5,5,5)), n.burnin=1000, n.sim=3000, n.thin=2)
 future <- subset(riverflows, Date>="2009-04-05")                  
 forecast <- forecasting(fit3.5,data=future,row.names=Date)

 dev.new()
 plot(forecast,n=300,historical=list(type="l",col="black",lty=1),
               forecasts=list(type="l",col="black",lty=3,ylab="River's flow"),
               forecasts.PI=list(col="light gray",border=NULL))
