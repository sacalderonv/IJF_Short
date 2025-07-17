####Resume 1000 replications
#####
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#####
#####Part 1:
#####Configure
#####
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
rm(list = ls())



regimes=3
n_rep=1000

####

###Load information, you should activate one of the following 6 lines and update the names in lines from 23 to 28 for the corresponding distribution error
#load("replicas_gaussian_1000_3reg1.rds")
load("replicas_student_1000_3reg1.rds")
#load("replicas_slash_1000_3reg1.rds")
#load("replicas_contaminated_1000_3reg1.rds")
#load("replicas_hyperbolic_1000_3reg1.rds")
#load("replicas_laplace_1000_3reg1.rds")


## list for Estimation of the selected distribution, you must change the name of the distribution in the next line in "repl_hyperbolic_estimation"
repl_estimation<-repl_student_estimation

##Flist for forecasting n.ahead with selected distribution, you must change the name of the distribution in the next line in "repl_hyperbolic"
repl<-repl_student



###Setting Structural Parameters. 
h.ahead=10
k = 2 ##Dimension of the output vector
v=0 ## Dimension of the exogenous variable vector
ars <- list(p=c(1,1,1),q=c(0,0,0),d=c(0,0,0))
delay <- 1
para.extra=TRUE #Depend if the distribution error has or nor extra parameter tou have to change to TRUE or FALSE



######non-structural Parameters. do not change 

### Regime 1
Location_R1 = list(phi1 = matrix(c(0.8,0,-0.2,0.5),k,k,byrow = TRUE))
Sigma_R1 = matrix(c(1,0,0,4),k,k,byrow = TRUE)

cs_1=matrix(c(2,1),nrow=k)

R1 = list(orders = list(p = 1,q = 0,d = 0),Location = Location_R1,Sigma = Sigma_R1,cs=cs_1)

###Regime 2
#Location_R2 = list(phi1 = matrix(c(0.3,0,0,-0.6),k,k,byrow = TRUE),phi2 = matrix(c(-0.1,0,0,0.3),k,k,byrow = TRUE))
Location_R2 = list(phi1 = matrix(c(0.3,0,0,-0.6),k,k,byrow = TRUE))
Sigma_R2 = matrix(c(1,0,0,1),k,k,byrow = TRUE)

cs_2=matrix(c(0.4,-2),nrow=k)

R2 = list(orders = list(p = 1,q = 0,d = 0),
          Location = Location_R2 ,Sigma = Sigma_R2,cs=cs_2)

###Regime 3

#Location_R3 = list(phi1 = matrix(c(0.6,0,-0.2,0.8),k,k,byrow = TRUE),phi2 = matrix(c(0.2,0,0,-0.1),k,k,byrow = TRUE))
Location_R3 = list(phi1 = matrix(c(0.6,0,-0.2,0.8),k,k,byrow = TRUE))


Sigma_R3 = matrix(c(2,0,0,1),k,k,byrow = TRUE)

cs_3=matrix(c(-3,0),nrow=k)




R3 = list(orders = list(p = 1,q = 0,d = 0),
          Location_ = Location_R3  ,Sigma = Sigma_R3,cs=cs_3)




## crea lista de objeto tipo Regime
Rg = list(R1 = R1,R2 = R2,R3=R3)

###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#####
#####Part 2:
#####Estimation
#####
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###Concatene if it is necessary. Use only in the the case when you have many files to get the total of replication 
#repl_estimation<-c(repl_estimation,repl_slash_estimation)
#repl<-c(repl,repl_slash)
###Lists to store the results for forecasting
fun = function(x){x/n_rep}
veces<-vector("list",h.ahead)
longitud<-vector("list",h.ahead)
sd<-vector("list",h.ahead)

#####Lists to store the results for estimation of the parameters

veces_est<-vector("list",regimes)
longitud_est<-vector("list",regimes)
sd_est<-vector("list",regimes)


###Estimación
for(reg in 1:regimes){
  dime_inter<-dim(Rg[[reg]]$cs)
  dime_phi<-dim(Rg[[reg]]$Location$phi1)
  dime_sigma<-dim(Rg[[reg]]$Sigma)
  suma_intercept<-matrix(rep(0,prod(dime_inter)),dime_inter[1],dime_inter[2])
  suma_phi<-matrix(rep(0,prod(dime_phi)),dime_phi[1],dime_phi[2])
  suma_sigma<-matrix(rep(0,prod(dime_sigma)),dime_sigma[1],dime_sigma[2])
  for(rep in 1:n_rep){
    
    ###Autoregressive
    row.names.obj_Reg_autoreg<-rownames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)
    col.names.obj_Reg_autoreg<-colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)
    pos.mean<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[1])
    pos.low<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[4])
    pos.high<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[5])
    pos.sd<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive)==col.names.obj_Reg_autoreg[2])
    
    ####IC
    ###low
    int.low_Reg<-repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive[,pos.low]
    int.low.intercep_Reg<-int.low_Reg[row.names.obj_Reg_autoreg[1],]
    int.low_Phi_Reg<-t(int.low_Reg[row.names.obj_Reg_autoreg[2:(k+1)],])
    
    ###high
    int.high_Reg<-repl_estimation[[rep]]$non_structural[[reg]]$Autoregressive[,pos.high]
    int.high.intercep_Reg<-int.high_Reg[row.names.obj_Reg_autoreg[1],]
    int.high_Phi_Reg<-t(int.high_Reg[row.names.obj_Reg_autoreg[2:(k+1)],])
    ####Comparación
    Intercept_TF<-(int.low.intercep_Reg<Rg[[reg]]$cs) & (Rg[[reg]]$cs<int.high.intercep_Reg)
    Include_Intercept<- matrix(as.integer(as.logical(Intercept_TF)),dime_inter[1],dime_inter[2])
    suma_intercept<-suma_intercept+Include_Intercept
    
    Phi_TF<-(int.low_Phi_Reg<Rg[[reg]]$Location$phi1) & (Rg[[reg]]$Location$phi1<int.high_Phi_Reg)
    Include_Phi<- matrix(as.integer(as.logical(Phi_TF)),dime_phi[1],dime_phi[2])
    suma_phi<-suma_phi+Include_Phi
    
    ###Sigma
    row.names.obj_Reg_scale<-rownames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)
    col.names.obj_Reg_scale<-colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)
    pos.mean_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[1])
    pos.low_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[2])
    pos.high_scale<-which(colnames(repl_estimation[[rep]]$non_structural[[reg]]$Scale)==col.names.obj_Reg_scale[3])
    
    ####IC
    ###low
    int.low_Reg_scale<-repl_estimation[[rep]]$non_structural[[reg]]$Scale[,pos.low_scale]
    
    int.low_Sigma_Reg<-t(int.low_Reg_scale)
    
    ###high
    int.high_Reg_scale<-repl_estimation[[rep]]$non_structural[[reg]]$Scale[,pos.high_scale]
    
    int.high_Sigma_Reg<-t(int.high_Reg_scale)
    
    
    ####Comparación
    
    Sigma_TF<-(int.low_Sigma_Reg<(Rg[[reg]]$Sigma)) & ((Rg[[reg]]$Sigma)<int.high_Sigma_Reg)
    Include_Sigma<- matrix(as.integer(as.logical(Sigma_TF)),dime_sigma[1],dime_sigma[2])
    suma_sigma<-suma_sigma+Include_Sigma
  }
  veces_est[[reg]]<-list(Intercept=suma_intercept,Phi=suma_phi,Sigma=suma_sigma)
}

if(isTRUE(para.extra)){
extra_1<-extra
suma_extra<-rep(0,length(extra_1))
names.extra<-colnames(repl_estimation[[1]]$extra)
pos.low.extra<-which(colnames(repl_estimation[[1]]$extra)==names.extra[4])
pos.high.extra<-which(colnames(repl_estimation[[1]]$extra)==names.extra[5])
sesgo_extra<-rep(0,length(extra_1))
}


suma_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
dife_thresholds<-matrix(rep(0,regimes-1),regimes-1,1)
r_true<-as.matrix(umbrales)
delay<-rep(0,n_rep)

for(rep in 1:n_rep){
  dife_thresholds=dife_thresholds+(repl_estimation[[rep]]$thresholds[,1]-r_true)
 
  delay[rep]<-repl_estimation[[rep]]$d
  
  if(isTRUE(para.extra)){
    for(l in 1:length(extra_1)){
  	sesgo_extra=sesgo_extra+(repl_estimation[[rep]]$extra[,l]-extra_1)
  
  if((repl_estimation[[rep]]$extra[l,pos.low.extra]<extra_1[l]) & ((repl_estimation[[rep]]$extra[l,pos.high.extra]>extra_1[l])) )
  {suma_extra[l]=suma_extra[l]+1}
    }
  }
  
  Threshold_TF<-r_true>repl_estimation[[rep]]$thresholds[,2] & r_true<repl_estimation[[rep]]$thresholds[,3]
  Include_Threshold<-matrix(as.integer(as.logical(Threshold_TF)),regimes-1,1)
  suma_thresholds=suma_thresholds+Include_Threshold
}



sesgo_thresholds<-dife_thresholds/n_rep




divided<-function(x){x/10}


prop_Reg1<-lapply(veces_est[[1]],fun)
prop_Reg2<-lapply(veces_est[[2]],fun)
prop_Reg3<-lapply(veces_est[[3]],fun)

###Percentage of time credible intervals captures real individual parameters table 1 respective column 
print("Table 3: Percentage of times that the true parameter values lie on the 95% credible intervals considering M2 with sample size 1000")

print(list("Regime 1", lapply(veces_est[[1]], divided)))
print(list("Regime 2", lapply(veces_est[[2]], divided)))
print(list("Regime 3", lapply(veces_est[[3]], divided)))



#suma_extra/n_rep
prop.thresholds<-suma_thresholds/n_rep
prop.delay<-table(delay)/n_rep

cat("Table 3 threshold parameter: ",prop.thresholds*100,"\n")
cat("Table 3 delay parameter:",prop.delay*100,"\n")
###Relative Bias 

if(isTRUE(para.extra)){
sesgo_extra_def<-sesgo_extra/n_rep
sesgo_relativo_extra<-(sesgo_extra_def/abs(extra))*100
prop.extra=suma_extra/n_rep
cat("Table 3 Percentage of times that extra parameter values lie on the 95% credible intervals considering M2 with sample size 1000,:",prop.extra*100,"\n")
cat("Table 4 Relative bias for extra parameter:",sesgo_relativo_extra,"\n")
}

sesgo_relativo_thresholds<-(sesgo_thresholds/abs(umbrales))*100
cat("Table 4 Relative bias for threshold parameter:",sesgo_relativo_thresholds,"\n")
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#####
#####Part 3:
#####Forecasting
#####
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################

for(pasos in 1:h.ahead){
  suma_pasos=rep(0,k)
  suma_longitud=rep(0,k)
  #suma_sd=rep(0,k)
  for(rep in 1:n_rep){
    int.low.pos=which(names(repl[[rep]]$forecasting[pasos,])=="HDI_Low")
    int.high.pos=which(names(repl[[rep]]$forecasting[pasos,])=="HDI_high")
    #sd.pos=which(names(repl[[rep]]$forecasting[pasos,])=="SD")
    for(entrada in 1:k){
      if((repl[[rep]]$True_Values[pasos,entrada]>repl[[rep]]$forecasting[pasos,int.low.pos[entrada]])&(repl[[rep]]$True_Values[pasos,entrada]<repl[[rep]]$forecasting[pasos,int.high.pos[entrada]])){suma_pasos[entrada]=suma_pasos[entrada]+1}
      suma_longitud[entrada]<-suma_longitud[entrada]+(repl[[rep]]$forecasting[pasos,int.high.pos[entrada]]-repl[[rep]]$forecasting[pasos,int.low.pos[entrada]])
      #suma_sd[entrada]<-suma_sd[entrada]+repl[[rep]]$forecasting[pasos,sd.pos[entrada]]
    }
    #print(suma_longitud_gauss)
    #print(suma_sd_gauss)
  }
  veces[[pasos]]=suma_pasos
  longitud[[pasos]]=suma_longitud
  #sd[[pasos]]=suma_sd
}


prop_pred_int<-lapply(veces,fun)

mult<-function(x){x*100}


###Corresponding column of table 6
print("Table 6: Percentage of times that all components of the vector y_{T+h} for h=1,...,10 lie in their 95% prediction intervals considering M2 with sample size 1000.")
print(lapply(prop_pred_int,mult),digits=3)


