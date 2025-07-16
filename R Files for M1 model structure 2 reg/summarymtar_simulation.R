summarymtar_simulation <- function(object, credible=0.95, digits=5,Get.results=FALSE,Print.Results=TRUE,...){
  k <- ncol(object$data[[1]]$y)
  n.sim <- object$n.sim
  
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic"))
  {
    
    myvec <- c("thresholds", "d","extra", "non_structural")
    salida_estimation<-setNames(vector("list", length = length(myvec)), myvec)
    salida_estimation$non_structural<-vector("list",object$regim)
    for(regime in 1:object$regim)
    {
      myvec1 <- c("Autoregressive", "Scale") 
      salida_estimation$non_structural[[regime]]<-setNames(vector("list", length = length(myvec1)), myvec1)
    }
    
  }
  else{
    myvec <- c("thresholds", "d", "non_structural")
    salida_estimation<-setNames(vector("list", length = length(myvec)), myvec)
    salida_estimation$non_structural<-vector("list",object$regim)
    for(regime in 1:object$regim)
    {
      myvec1 <- c("Autoregressive", "Scale") 
      salida_estimation$non_structural[[regime]]<-setNames(vector("list", length = length(myvec1)), myvec1)
    }
  }
  resumen <- function(x){
    x <- matrix(x,ifelse(is.null(nrow(x)),1,nrow(x)),ifelse(is.null(ncol(x)),length(x),ncol(x)))
    y <- matrix(0,nrow(x),5) 
    y[,1] <- apply(x,1,mean)
    y[,2] <- apply(x,1,sd)
    y[,3] <- apply(x,1,function(x) min(mean(sign(median(x))*x > 0),1-1/200000))
    y[,3] <- 2*(1 - y[,3])
    ks <- seq(credible,1,length=n.sim*(1-credible))
    lis <- t(apply(x,1,quantile,probs=ks-credible))
    lss <- t(apply(x,1,quantile,probs=ks))
    dif <- apply(abs(lss-lis),1,which.min)
    y[,4] <- lis[cbind(1:nrow(x),dif)]
    y[,5] <- lss[cbind(1:nrow(x),dif)]
    colnames(y) <- c("   Mean"," Std.Dev"," 2(1-PD) ","HDI_Low","HDI_high")
    return(y)
  }
  
  if(object$regim > 1){
    thresholds <- matrix(round(resumen(matrix(object$chains$thresholds,nrow=object$regim-1)),digits=digits)[,c(1,4,5)],ncol=3)
    h <- round(mean(object$chains$h),digits=0)
    thresholds1 <- paste0(c("(-Inf",paste0("(",round(thresholds[,1],digits=digits))),",",c(paste0(round(thresholds[,1],digits=digits),"]"),"Inf)"))
    thresholds2 <- paste0(c("(-Inf",paste0("(",round(thresholds[,2],digits=digits))),",",c(paste0(round(thresholds[,2],digits=digits),"]"),"Inf)"))
    thresholds3 <- paste0(c("(-Inf",paste0("(",round(thresholds[,3],digits=digits))),",",c(paste0(round(thresholds[,3],digits=digits),"]"),"Inf)"))
    d <- data.frame(cbind(thresholds1,thresholds2,thresholds3))
    rownames(d) <- paste("Regime",1:nrow(d))
    colnames(d) <- c("Thresholds (mean,","HDI_low,","HDI_high)")
    salida_estimation$thresholds<-thresholds
    salida_estimation$d<-h
  }
  if(Print.Results)
  {
    cat("\nResponse          :",ifelse(length(colnames(object$data[[1]]$y))==1,colnames(object$data[[1]]$y),paste(colnames(object$data[[1]]$y),collapse="    |    ")))
    if(object$regim > 1) cat("\nThreshold series  :",object$ts,"(mean)")	  
    cat("\nError distribution:",object$dist)
    cat("\n\n")
  }
  if(object$regim > 1) if(Print.Results)
  {print(d)}
  for(i in 1:object$regim){
    out <- outs <- vector()
    for(j in 1:k){
      temp <- object$chains[[i]]$location[,seq(j,n.sim*k,k)]
      temps <- matrix(matrix(object$chains[[i]]$scale,k,n.sim*k)[,seq(j,n.sim*k,k)],nrow=k)
      if(j > 1){
        out <- cbind(out,matrix(0,nrow(out),1),round(resumen(temp),digits=digits))
        outs <- cbind(outs,round(resumen(temps)[,c(1,4,5)],digits=digits+1))
      }  
      else{
        out <- round(resumen(temp),digits=digits)
        outs <- round(resumen(temps)[,c(1,4,5)],digits=digits+1)
      }
    }
    ########Tomar Salidas####
    rownames(out) <- object$name[[i]]
    salida_estimation$non_structural[[i]]$Autoregressive<-out
    salida_estimation$non_structural[[i]]$Scale<-outs
    outs <- matrix(cbind(outs,0),k,3*k+1)
    outs <- outs[,c(seq(1,3*k,3),3*k+1,seq(2,3*k,3),3*k+1,seq(3,3*k,3))]
    rownames(outs) <- colnames(object$data[[1]]$y)
    colnames(outs) <- c(rownames(outs),"",rownames(outs),"",rownames(outs))
    if(Print.Results)
    {
      cat("\n\nRegime",i,":")
      cat("\nAutoregressive coefficients\n")
      print(format(out, justify = "right", format = "+/-", zero.print="   |   "), quote=FALSE)
      cat("\nScale parameter (mean, HDI_low, HDI_high)\n")
      print(format(outs, justify = "right", format = "+/-", zero.print="   ."), quote=FALSE)
    }
  }
  if(object$dist %in% c("Slash","Contaminated normal","Student-t","Hyperbolic")){
    out <- round(resumen(object$chains$extra),digits=digits+1)
    salida_estimation$extra<-out
    out[,3] <- 0
    if(object$dist %in% c("Slash","Student-t","Hyperbolic")) rownames(out)[nrow(out)] <- paste0("nu",paste0(rep("",max(nchar(object$name[[1]]))-1),collapse=" "))
    else rownames(out)[nrow(out):(nrow(out)-1)] <- paste0(c("nu2","nu1"),paste0(rep("",max(nchar(object$name[[1]]))-2),collapse=" "))
    if(Print.Results)
    {
      cat("\n\nExtra parameter","\n")
      print(format(out, justify = "right", flag="+", zero.print="   .   "), quote=FALSE)
    }
  }
  if(Print.Results)
  {
    cat("\n\n")
  }
  if(Get.results){return(salida_estimation)}
}

