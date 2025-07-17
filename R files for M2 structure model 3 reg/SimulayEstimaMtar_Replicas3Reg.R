####Program for simulating MTAR process and estimating its parameters
###This only includes M2 structure model, that is bivariate a MTAR with 3 regimes
###This an example for simulating and estimating a MTAR model with distribution error Slash with parameter extra 4.


# Required packages====

library(GIGrvg)
library(Formula)
library(Rfast)
library(mtarm)
library(stats)

####Defining parameters of the model and  lists for storing the results====
n_rep=1000 ##number of replications

#repl_student_estimation<-vector("list", n_rep)
#repl_gaussian_estimation<-vector("list", n_rep)
#repl_laplace_estimation<-vector("list", n_rep)
#repl_hyperbolic_estimation<-vector("list", n_rep)
repl_slash_estimation<-vector("list", n_rep)
#repl_contaminated_estimation<-vector("list", n_rep)


#repl_student<-vector("list", n_rep)
#repl_gaussian<-vector("list", n_rep)
#repl_laplace<-vector("list", n_rep)
#repl_hyperbolic<-vector("list", n_rep)
repl_slash<-vector("list", n_rep)
#repl_contaminated<-vector("list", n_rep)

calen=100 ## Starting values
long=1000 ##Long of Series
h.ahead=10 ## forecast horizon
Tlen = long+calen+h.ahead




##Parameters of the  Process Yt ====
k = 2 #number of components of the output vector
v=0


###Function for simulating arma model with constant for threshold process
simulate_arma <- function(ar_params = NULL, ma_params = NULL, constant = 0, n = 100, sigma2 = 1){ 
  
  # Load necessary library
  library(stats)
  
  # Validate input
  if (is.null(ar_params) && is.null(ma_params)) {
    stop("Please specify either AR or MA parameters.")
  }
  
  # Generate innovation series
  innovations <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # Initialize series
  arma_series <- numeric(n)
  
  # Handle initial values for AR components
  if (!is.null(ar_params)) {
    arma_series[1:max(length(ar_params), 1)] <- innovations[1:max(length(ar_params), 1)]
  }
  
  # Simulate the ARMA process
  for (t in (max(length(ar_params), 1) + 1):n) {
    ar_part <- ifelse(is.null(ar_params), 0, sum(ar_params * arma_series[(t - 1):(t - length(ar_params))]))
    ma_part <- ifelse(is.null(ma_params), 0, sum(ma_params * innovations[(t - 1):(t - length(ma_params))]))
    arma_series[t] <- constant + ar_part + ma_part + innovations[t]
  }
  return(arma_series)
}



###Structural parameters and parameters for the exogenous and threshold processes
dist <- "Slash"
#extra<-c(0.1,0.15)
extra=4
delay <- 1 #Delay parameter
Intercept <- TRUE
ars <- list(p=c(1,1,1),q=c(0,0,0),d=c(0,0,0))  #Autoregressive orders for each regime
ars_aux<-list(p=c(1,1,1)) 

inic<-calen

###Parameters for the threshold process
phi_ar=0.6
sigma2_inn_ar<-1
constant<-1


### Paramaters of the R1 regime ====
Location_R1 = list(phi1 = matrix(c(0.8,0,-0.2,0.5),k,k,byrow = TRUE))

Sigma_R1 = matrix(c(1,0,0,4),k,k,byrow = TRUE)

cs_1=matrix(c(2,1),nrow=k)

R1 = list(orders = list(p = 1,q = 0,d = 0),Location = Location_R1,Sigma = Sigma_R1,cs=cs_1)

###Parameters of the R3 regimen
Location_R2 = list(phi1 = matrix(c(0.3,0,0,-0.6),k,k,byrow = TRUE))

Sigma_R2 = matrix(c(1,0,0,1),k,k,byrow = TRUE)

cs_2=matrix(c(0.4,-2),nrow=k)

R2 = list(orders = list(p = 1,q = 0,d = 0),
          Location = Location_R2 ,Sigma = Sigma_R2,cs=cs_2)

###Parameters of the R3 regimen

Location_R3 = list(phi1 = matrix(c(0.6,0,-0.2,0.8),k,k,byrow = TRUE))



Sigma_R3 = matrix(c(2,0,0,1),k,k,byrow = TRUE)

cs_3=matrix(c(-3,0),nrow=k)




R3 = list(orders = list(p = 1,q = 0,d = 0),
          Location_ = Location_R3  ,Sigma = Sigma_R3,cs=cs_3)




## list of the 
Rg = list(R1 = R1,R2 = R2,R3=R3)


####Mean and Variace ARMA
mean_arma=constant/(1-sum(phi_ar))

var_arma<-ltsa::tacvfARMA(phi=phi_ar,sigma2=sigma2_inn_ar)[1]

params <- list()
umbrales<-qnorm(p=c(0.33,0.66),mean =mean_arma,sd = sqrt(var_arma) ) ##Thoeretical thresholds 


###Creates a list of structural parameters to be used for the simulation 

for(i in 1:length(ars$p)){
  np <- Intercept + ars$p[i]*k + ars$q[i]*v + ars$d[i]
  params[[i]] <- list()
  params[[i]]$location <-rbind(t(Rg[[i]][[4]]),matrix(unlist(Rg[[i]][[2]]),ncol=k,byrow=TRUE))
  #params[[i]]$location <- matrix(c(rbeta(np*k,shape1=4,shape2=16)),np,k)
  # params[[i]]$scale <- diag(rgamma(k,shape=1,scale=1))
  params[[i]]$scale <- Rg[[i]][[3]]
  params[[i]]$scale2 <- chol(params[[i]]$scale)
}
params




start.time <- Sys.time()
bar <- txtProgressBar(min=0, max=n_rep, initial=0, char="=", style=3)
for(replicas in 1:n_rep){
  ## 3 dimensions ====
  
  ## Simulating the process Zt
  
  
  Zt<-tseries  <- simulate_arma(ar_params = phi_ar, constant = constant, sigma2 = sigma2_inn_ar,n = Tlen+max(ars$p,ars$q,ars$d,delay))
  
  myseries <- matrix(0,Tlen+max(ars$p,ars$q,ars$d,delay),k)
  
  regimen <- cut(Zt[(max(ars$p,ars$q,ars$d,delay)+1-delay):(length(Zt)-delay)],breaks=c(-Inf,umbrales,Inf),labels=1:length(ars$p))
  myseries[1:max(ars$p,ars$q,ars$d),] <- rnorm(max(ars$p,ars$q,ars$d)*k)
  
  for(i in 1:Tlen){
    current <- max(ars$p,ars$q,ars$d,delay) + i
    regimeni <- regimen[i]
    if(Intercept) X <- 1 else X <-  vector()
    for(j in 1:ars$p[regimeni]) X <- c(X,myseries[current-j,])
    if(ars$q[regimeni] > 0) for(j in 1:ars$q[regimeni]) X <- c(X,Xt[current-j,])
    if(ars$d[regimeni] > 0) for(j in 1:ars$d[regimeni]) X <- c(X,Zt[current-j,])
    Theta <- params[[regimeni]]$location   
    mu <- apply(matrix(X,nrow(Theta),ncol(Theta))*Theta,2,sum)
    u <- 1
    if(dist=="Student-t")  u <- 1/rgamma(1,shape=extra/2,rate=extra/2)
    if(dist=="Slash")  u <- 1/rbeta(1,shape1=extra/2,shape2=1)
    if(dist=="Contaminated normal")  if(runif(1)<=extra[1]) u <- 1/extra[2]
    if(dist=="Laplace")  u <- rexp(1,rate=1/8)
    if(dist=="Hyperbolic") u <- rgig(n=1,lambda=1,chi=1,psi=extra^2)
    myseries[current,] <- t(params[[regimeni]]$scale2)%*%matrix(rnorm(k,mean=0,sd=sqrt(u)),k,1) + matrix(mu,k,1)
  }
  
  datos <- data.frame(myseries[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt),])
  colnames(datos) <- paste("X",1:k,sep="")
  datos <- data.frame(datos,threshold=Zt[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt)])
  if(max(ars$q)>0){
    colnames(Xt) <- paste("X",(k+1):(k+v),sep="")
    datos <- data.frame(datos,Xt[(max(ars$p,ars$q,ars$d,delay)+inic+1):length(Zt),])
  }
  
  
  ### Creating dates for simulating time serie data====
  Fechas=seq(as.Date("2000/1/1"), by = "day", length.out = (Tlen-calen))
  
  
  #### 3 reg defining long of series ====
  
  datos1=data.frame(datos,Fecha=Fechas)
  
  
    fecha_final<-Fechas[long]
    
    ##  Procedure for estimating parameters using function mtar of the mtarm R package====

  ##########
  
  fit1_Bayes <- mtar(~X1+X2|threshold, data=datos1, row.names=Fecha,subset={Fecha<=fecha_final} ,ars=ars_aux, dist=dist, n.burnin=1000, n.sim=1500, n.thin=1)
  ###Function summarymtar_simulation must be ran previously
  summary_Bayes<-summarymtar_simulation(fit1_Bayes,Get.results = TRUE,Print.Results = FALSE)
  
  
  nano <- forecasting(fit1_Bayes,subset(datos1,Fecha > fecha_final),row.names=Fecha)
  forecasting_salida<-nano$summary
  repl_slash[[replicas]]=list(forecasting=forecasting_salida,True_Values=subset(datos1,Fecha > fecha_final))
  repl_slash_estimation[[replicas]]<-summary_Bayes
  
  
  setTxtProgressBar(bar,replicas)
  }

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
####Save Replications====

save(repl_slash_estimation,repl_slash,extra,umbrales,Rg,file="replicas_slash_1000_3reg1.rds")

