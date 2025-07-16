 #Results are very similar to those reported in the paper, but not identical, due
 #to the lack of seeds specified prior to random number generation.

 library(mtarm)
 data(riverflows)
 fit3.5 <- mtar( ~ Bedon + LaPlata | Rainfall, row.names=Date, dist="Laplace",
	               data=riverflows, subset={Date<"2009-04-05"}, 
	               ars=list(p=c(5,5,5)), n.burnin=1000, n.sim=3000, n.thin=2)
 a <- coef(fit3.5,FUN=mean)
 table.out <- vector()
 p <- max(fit3.5$ars$p)
 k <- ncol(fit3.5$data[[1]]$y)
 l <- fit3.5$regim
 for(i in 1:l){
	temp <- vector()
	if(p > fit3.5$ars$p[i]) b <- cbind(matrix(NA,k,k),t(a[[i]]$location),matrix(NA,k,k*(p-fit3.5$ars$p[i])),a[[i]]$scale)
	else b <- cbind(matrix(NA,k,k),t(a[[i]]$location),a[[i]]$scale)
	if(fit3.5$Intercept) b <- b[,-1]
	for(ii in 1:(fit3.5$Intercept+p+1)) temp <- rbind(temp,b[,c(1:k)+k*(ii-1)])
	table.out <- cbind(table.out,temp)
 }
 colnames(table.out) <- rep("",l*k)
 colnames(table.out)[seq(1,l*k,k)] <- paste(rep("Regime",l),1:l)
 rows <- rep("",p*k+k)
 rows[seq(1,p*k+k,k)] <- c(paste0("lag(",1:p,")"),"Scale")
 if(fit3.5$Intercept){
   rows <- c(rep("",k),rows)
   rows[1] <- "Intercept" 
 }  
 if(l > 1){
   table.out <- rbind(table.out,cbind(rbind(NA,a$thresholds,a$delay),matrix(NA,nrow(a$thresholds)+2,ncol(table.out)-1)))
   rows <- c(rows,"",rownames(a$thresholds),"delay") 
 }
 rownames(table.out) <- rows
 table.out <- round(table.out,digits=4)
 
 print(table.out, na.print="")
