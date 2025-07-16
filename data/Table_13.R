 #Results are very similar to those reported in the paper, but not identical, due
 #to the lack of seeds specified prior to random number generation.

 library(mtarm)
 data(riverflows)
 fit3.5 <- mtar( ~ Bedon + LaPlata | Rainfall, row.names=Date, dist="Laplace",
	               data=riverflows, subset={Date<"2009-04-05"}, 
	               ars=list(p=c(5,5,5)), n.burnin=1000, n.sim=3000, n.thin=2)
 future <- subset(riverflows, Date>="2009-04-05")                  
 forecast <- forecasting(fit3.5,data=future,row.names=Date,credible=0.95)

 a <- cbind(future[,c("Date","Bedon","LaPlata")],forecast$summary)
 colnames(a)[c(2,3,4,7)] <- c("True","True","Forecast","Forecast")
 b <- a[,c(1,2,4,5,6)]
 c <- a[,c(1,3,7,8,9)]
 d <- rbind(b,c)
 table.out <- d[order(d$Date),]
 table.out[,-1] <- round(table.out[,-1],digits=3)
 rownames(table.out) <- NULL
 
 print(table.out)
