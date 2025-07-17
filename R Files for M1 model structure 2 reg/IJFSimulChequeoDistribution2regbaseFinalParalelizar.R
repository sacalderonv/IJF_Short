####Program for simulating, and estimating a MTAR with 2 regimes model
####also Compare with other distribution for errors via DIC and WAIC 


# Packages====

library(GIGrvg)
library(Formula)
library(Rfast)
library(tsDyn)
library(mtarm)
library(foreach)
library(doParallel)

# Detect the number of available cores
n_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)


####Defining parameters of the model and  lists for storing the results====
n_rep=1000



####Storing forecast
repl_student<-vector("list", n_rep)  #Store forcasting for TRUE distribution
#repl_gaussian<-vector("list", n_rep)
#repl_laplace<-vector("list", n_rep)
#repl_hyperbolic<-vector("list", n_rep)
#repl_slash<-vector("list", n_rep)
#repl_contaminated<-vector("list", n_rep)


#Store estimation for TRUE and alternative distributions
repl_student_estimation<-vector("list", n_rep)
repl_gaussian_estimation<-vector("list", n_rep)
repl_laplace_estimation<-vector("list", n_rep)
repl_hyperbolic_estimation<-vector("list", n_rep)
repl_slash_estimation<-vector("list", n_rep)
repl_contaminated_estimation<-vector("list", n_rep)

##Store WAIC and DIC for comparing distribution in missespecification
repl_student_compare_dist<-vector("list", n_rep) 
#repl_gaussian_compare_dist<-vector("list", n_rep)
#repl_laplace_compare_dist<-vector("list", n_rep)
#repl_hyperbolic_compare_dist<-vector("list", n_rep)
#repl_slash_compare<-vector("list", n_rep)
#repl_contaminated_compare_dist<-vector("list", n_rep)
calen=100
long=1000
h.ahead=10
Tlen = long+calen+h.ahead
inic<-calen
######


dimeUt=3
Sigma_ut = diag(2,dimeUt,dimeUt) 
Phi_ut = list(phi1 = matrix(c(0.24,0.48,-0.12,0.46,-0.36,0.1,-0.12,-0.47,0.58),dimeUt,dimeUt,byrow = TRUE))
cs_var<-matrix(c(0,0,2),nrow=dimeUt)
R_ut = list(orders = list(p = 1),Phi = Phi_ut,Sigma = Sigma_ut,cs =cs_var )
mean_VAR=solve(diag(3)-R_ut$Phi$phi1)%*%R_ut$cs



##Parameters of the process Yt ====
k = 3
v=2


dist <- "Student-t"
#extra <- c(0.05,0.1)
#extra=0.11
extra=3
delay <- 0

ars <- list(p=c(1,2),q=c(1,0),d=c(1,0))


Intercept <- TRUE
### R1 regime ====
#Location_R1 = list(phi1 = matrix(c(0.1,0.6,0.4,-0.4,0.5,-0.7,0.2,0.6,-0.3),k,k,byrow = TRUE)) 
Location_R1 = list(phi1 = matrix(c(0.1,0.6,0.4,-0.4,0.5,-0.7,0.2,0.6,-0.3),k,k,byrow = TRUE),beta1=matrix(c(0.6,-0.5,-0.4,0.6,0.1,0.3),k,v,byrow = TRUE),delta1=matrix(c(0.6,1,-0.4),k,1,byrow = TRUE)) 


Sigma_R1 = matrix(c(1,0.3,0.4,0.3,1,-0.5,0.4,-0.5,1),k,k,byrow = TRUE) 

cs_1=matrix(c(1,-2,6),nrow=k) 

R1 = list(orders = list(p = 1,q = 1,d = 1),Location = Location_R1,Sigma = Sigma_R1,cs=cs_1)

### R2 regime ====
#Location_R2 = list(phi1 = matrix(c(0.3,0.5,-0.5,0.2,0.7,-0.1,0.3,-0.4,0.6),k,k,byrow = TRUE))
Location_R2 = list(phi1 = matrix(c(0.3,0.5,-0.5,0.2,0.7,-0.1,0.3,-0.4,0.6),k,k,byrow = TRUE),phi2 = matrix(c(0.3,0.1,0.2,-0.2,-0.6,0.4,0.3,-0.1,0.5),k,k,byrow = TRUE))

Sigma_R2 = matrix(c(1.5,0.2,-0.4,0.2,1,0.7,-0.4,0.7,2),k,k,byrow = TRUE)

cs_2=matrix(c(5,-3,-1),nrow=k)

R2 = list(orders = list(p = 2,q = 0,d = 0),
          Location = Location_R2,Sigma = Sigma_R2,cs=cs_2)






Rg = list(R1 = R1,R2 = R2) # 2 reg
umbrales = mean_VAR[3] # 2 reg
params <- list()

for(i in 1:length(ars$p)){
  np <- Intercept + ars$p[i]*k + ars$q[i]*v + ars$d[i]
  params[[i]] <- list()
  params[[i]]$location <-rbind(t(Rg[[i]][[4]]),matrix(unlist(Rg[[i]][[2]]),ncol=k,byrow=TRUE))

  params[[i]]$scale <- Rg[[i]][[3]]
  params[[i]]$scale2 <- chol(params[[i]]$scale)
}


params
start.time <- Sys.time()
bar <- txtProgressBar(min=0, max=n_rep, initial=0, char="=", style=3)
for(replicas in 1:n_rep){

  
  ## Obtain process Ut=(Xt,Zt)
  
  Ut = tsDyn::VAR.sim(B=cbind(R_ut$cs,R_ut$Phi$phi1), n = Tlen+max(ars$p,ars$q,ars$d,delay) ,lag=R_ut$orders$p,include="const" ,varcov  = R_ut$Sigma)
  Zt = as.matrix(Ut[,(v+1)])
  Xt= Ut[,1:v] # Procesos Bidimensional v=2
  myseries <- matrix(0,Tlen+max(ars$p,ars$q,ars$d,delay),k)
  
  
  regimen <- cut(Zt[(max(ars$p,ars$q,ars$d,delay)+1-delay):(length(Zt)-delay)],breaks=c(-Inf,umbrales,Inf),labels=1:length(ars$p))
  myseries[1:max(ars$p,ars$q,ars$d),] <- rnorm(max(ars$p,ars$q,ars$d)*k)
  
  ##Simulate output process Y_t
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
  
  
  #plot(as.ts(datos[,1:k]))
  
  
  Fechas=seq(as.Date("2000/1/1"), by = "day", length.out = (Tlen-calen))
  
  
 
  
  datos1=data.frame(datos,Fecha=Fechas)
  
  
  

  fecha_final<-Fechas[long]
  
  
  
  
  
  
  
  ####Paralelizing
  
  
  
  results <- foreach(dist = c("Gaussian", "Student-t", "Slash", "Contaminated normal", "Laplace", "Hyperbolic")) %dopar% {
                       # Fit the model
                       fit <- mtarm::mtar(~X1+X2+X3|threshold|X4+X5, 
                                   data=datos1, 
                                   ars=ars, 
                                   dist=dist, 
                                   row.names=Fecha, 
                                   subset={Fecha<=fecha_final}, 
                                   n.burnin=500, 
                                   n.sim=1500, 
                                   n.thin=1, 
                                   Intercept=Intercept)
                       
                       # Calculate DIC and WAIC
                       dic_value <- mtarm::DIC(fit)
                       waic_value <- mtarm::WAIC(fit)
                       
                       # Return a list containing the fit, DIC, and WAIC
                       lista<-list(fit = fit, dic = setNames(dic_value, dist), waic = setNames(waic_value, dist))
                     }
 
  
  
  fits <- lapply(results, function(x) x$fit)
  dic_values <- sapply(results, function(x) x$dic)
  waic_values <- sapply(results, function(x) x$waic)
  
  fit1_2reg_Gaussian <- fits[[1]]
  fit1_2reg_Student <- fits[[2]]
  fit1_2reg_Slash <- fits[[3]]
  fit1_2reg_Contaminated <- fits[[4]]
  fit1_2reg_Laplace <- fits[[5]]
  fit1_2reg_Hyperbolic <- fits[[6]]
  
  
  
  DIC_models_dist<-as.matrix(t(dic_values))
  WAIC_models_dist<-as.matrix(t(waic_values))
  
  
  min_index_DIC_dist <- which.min(DIC_models_dist)
  min_index_WAIC_dist <- which.min(WAIC_models_dist)
  ######
  min_name_DIC_dist <- colnames(DIC_models_dist)[min_index_DIC_dist]
  min_value_DIC_dist <- DIC_models_dist[min_index_DIC_dist]
  
  min_name_WAIC_dist <- colnames(WAIC_models_dist)[min_index_WAIC_dist]
  min_value_WAIC_dist <- WAIC_models_dist[min_index_WAIC_dist]
  
  repl_student_compare_dist[[replicas]]<-list(DIC_models_dist=DIC_models_dist,WAIC_models_dist=WAIC_models_dist,Min_DIC_dist=c(min_name_DIC_dist,min_value_DIC_dist), Min_WAIC_dist=c(min_name_WAIC_dist,min_value_WAIC_dist) )
  #####
  
  
  
  summary_Bayes_student<-summarymtar_simulation(fit1_2reg_Student,Get.results = TRUE,Print.Results = FALSE)
  summary_Bayes_gaussian<-summarymtar_simulation(fit1_2reg_Gaussian,Get.results = TRUE,Print.Results = FALSE)
  summary_Bayes_slash<-summarymtar_simulation(fit1_2reg_Slash,Get.results = TRUE,Print.Results = FALSE)
  summary_Bayes_laplace<-summarymtar_simulation(fit1_2reg_Laplace,Get.results = TRUE,Print.Results = FALSE)
  summary_Bayes_contaminated<-summarymtar_simulation(fit1_2reg_Contaminated,Get.results = TRUE,Print.Results = FALSE)
  summary_Bayes_hyperbolic<-summarymtar_simulation(fit1_2reg_Hyperbolic,Get.results = TRUE,Print.Results = FALSE)
  
  
  nano <- forecasting(fit1_2reg_Student,subset(datos1,Fecha > fecha_final),row.names=Fecha)
  forecasting_salida<-nano$summary
  repl_student[[replicas]]=list(forecasting=forecasting_salida,True_Values=subset(datos1,Fecha > fecha_final)) ##Get the forecasting
  repl_student_estimation[[replicas]]<-summary_Bayes_student
  repl_gaussian_estimation[[replicas]]<-summary_Bayes_gaussian
  repl_slash_estimation[[replicas]]<-summary_Bayes_slash
  repl_laplace_estimation[[replicas]]<-summary_Bayes_laplace
  repl_contaminated_estimation[[replicas]]<-summary_Bayes_contaminated
  repl_hyperbolic_estimation[[replicas]]<-summary_Bayes_hyperbolic
  
  setTxtProgressBar(bar,replicas)
  #print(replicas)
}

stopCluster(cl)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken



###True distribution for errors ----
repl_estimation<-repl_student_estimation
repl_compare_dist_todo<-repl_student_compare_dist
repl<-repl_student
#### Compare distributions ----
repl_estimation_oth_dist_Gaussian<-repl_gaussian_estimation
repl_estimation_oth_dist_slash<-repl_slash_estimation
repl_estimation_oth_dist_contaminated<-repl_contaminated_estimation
repl_estimation_oth_dist_hyperbolic<-repl_hyperbolic_estimation
repl_estimation_oth_dist_laplace<-repl_laplace_estimation


param.extra=TRUE  #put TRUE or FALSE depends if the true distribution error has or not extra parameter
if(param.extra){save(repl,repl_estimation,repl_compare_dist_todo,repl,repl_estimation_oth_dist_Gaussian,repl_estimation_oth_dist_slash,repl_estimation_oth_dist_contaminated,repl_estimation_oth_dist_hyperbolic,repl_estimation_oth_dist_laplace,umbrales,extra,Rg,ars,h.ahead,Rg,file="replicas_Identify_Dist_RelativeBias_Forecasting_1000_2reg_ForStudent_1.rds")}else
{save(repl,repl_estimation,repl_compare_dist_todo,repl,repl_estimation_oth_dist_Gaussian,repl_estimation_oth_dist_slash,repl_estimation_oth_dist_contaminated,repl_estimation_oth_dist_hyperbolic,repl_estimation_oth_dist_laplace,umbrales,Rg,ars,h.ahead,file="replicas_Identify_Dist_RelativeBias_Forecasting_1000_2reg_ForStudent_1.rds")}

