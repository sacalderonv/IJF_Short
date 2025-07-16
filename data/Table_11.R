 #Results are very similar to those reported in the paper, but not identical, due
 #to the lack of seeds specified prior to random number generation.
 #The processing time is approximately 75 minutes.

 library(mtarm)
 data(riverflows)
 DIC.all <- WAIC.all <- RMSE.all <- Log.score.all <- vector()
 future <- subset(riverflows, Date>"2009-03-25" & Date<"2009-04-05")                  
 myfunc <- function(x){
	  temp <- forecasting(eval(parse(text=x)),data=future,row.names=Date,out.of.sample=TRUE)
	  return(c(sqrt(mean(temp$SE)),mean(temp$log.score)))
 }
 fits <- matrix(paste0(paste0("fit", rep(1:3, each=5)), ".", 1:5))
 dists <- c("Gaussian","Student-t","Slash","Contaminated normal","Hyperbolic","Laplace")
 
 for(dist in dists){
	fit1.1 <- mtar( ~ Bedon + LaPlata | Rainfall, row.names=Date, dist=dist,
	                  data=riverflows, subset={Date>="2006-01-05" & Date<="2009-03-25"}, 
	                  ars=list(p=1), n.burnin=1000, n.sim=3000, n.thin=2)
	fit1.2 <- update(fit1.1, ars=list(p=2), subset={Date>="2006-01-04" & Date<="2009-03-25"})
	fit1.3 <- update(fit1.1, ars=list(p=3), subset={Date>="2006-01-03" & Date<="2009-03-25"})
	fit1.4 <- update(fit1.1, ars=list(p=4), subset={Date>="2006-01-02" & Date<="2009-03-25"})
	fit1.5 <- update(fit1.1, ars=list(p=5), subset={Date<="2009-03-25"})
	fit2.1 <- update(fit1.3, ars=list(p=c(1,1)))
	fit2.2 <- update(fit1.3, ars=list(p=c(2,2)))
	fit2.3 <- update(fit1.3, ars=list(p=c(3,3)))
	fit2.4 <- update(fit1.4, ars=list(p=c(4,4)))
	fit2.5 <- update(fit1.5, ars=list(p=c(5,5)))
	fit3.1 <- update(fit1.3, ars=list(p=c(1,1,1)))
	fit3.2 <- update(fit1.3, ars=list(p=c(2,2,2)))
	fit3.3 <- update(fit1.3, ars=list(p=c(3,3,3)))
	fit3.4 <- update(fit1.4, ars=list(p=c(4,4,4)))
	fit3.5 <- update(fit1.5, ars=list(p=c(5,5,5)))
	DICs <- c(DIC(fit1.1,fit1.2,fit1.3,fit1.4,fit1.5,verbose=FALSE),
	          DIC(fit2.1,fit2.2,fit2.3,fit2.4,fit2.5,verbose=FALSE),
	          DIC(fit3.1,fit3.2,fit3.3,fit3.4,fit3.5,verbose=FALSE))
    DIC.all <- cbind(DIC.all,DICs)			   
	WAICs <- c(WAIC(fit1.1,fit1.2,fit1.3,fit1.4,fit1.5,verbose=FALSE),
	           WAIC(fit2.1,fit2.2,fit2.3,fit2.4,fit2.5,verbose=FALSE),
	           WAIC(fit3.1,fit3.2,fit3.3,fit3.4,fit3.5,verbose=FALSE))
    WAIC.all <- cbind(WAIC.all,WAICs)				
	oos <- apply(fits, 1, myfunc)
	RMSE.all <- cbind(RMSE.all,oos[1,])
	Log.score.all <- cbind(Log.score.all,oos[2,])
 }
 table.out <- round(rbind(DIC.all,WAIC.all,RMSE.all,Log.score.all),digits=2)
 table.out <- cbind(rep(1:3,each=5,times=4),rep(1:5,times=12),table.out)
 colnames(table.out) <- c("l","p*",dists)
 rownames(table.out) <- rep(c("DIC","WAIC","RMSE","log-score"),each=15)

 print(table.out)
