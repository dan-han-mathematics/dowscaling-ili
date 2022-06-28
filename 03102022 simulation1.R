####Dynamic downscaling and daily nowcasting from influenza surveillance data 
##############################
##########################################################################################################
#                               Dan Han, University of Louisiville, USA                         #######
#                                   Email: dan.han @louisville.edu
#                                         Version 1, March 10th 2022                              #######
##########################################################################################################
##########################################################################################################
# Purpose: This file uses dataset #1 that contains simulated CDC weekly ili rates, constitutional daily rates and respiratory
# daily rates to estimate downscaled daily ili rates based on the model in reference (1). And calculate Mean squred error (MSE)
# between estimates and the true value
##########################################################################################################

# Preinstalled packages: bsts, matlib, Semipar, Hmisc, MASS, mvtnorm, ggplot2, corpcor, RGeode, gridExtra, autostsm, lubridate, stats
# If you don't have these packages installed, please use install.packages("file_path\\package_file_name.extension",
# repos = NULL, type = "source") to install these packages in advance.
##########################################################################################################

# License: GPL (GNU Publice License),https://www.gnu.org/licenses/gpl-3.0.en.html
##########################################################################################################

# references:
# 1. Paul R, Dan H, DeDoncker E, Prieto D. Dynamic downscaling and daily nowcasting from influenza surveillance data. Statistics in Medicine. 2022; doi: 10.1002/sim.9502
# 2. Durbin, J. and Koopman, S.J., 2012. Time series analysis by state space methods (Vol. 38). OUP Oxford.
##########################################################################################################

# load packages: bsts, matlib, Semipar, Hmisc, MASS, mvtnorm, ggplot2, corpcor, RGeode, gridExtra, autostsm, lubridate, stats
library("bsts")     
library("matlib")
library(SemiPar)
library(Hmisc)
library(MASS)
library(mvtnorm)
library(ggplot2)
library("corpcor")
library("RGeode")
library(gridExtra)
library(autostsm)
library(lubridate)
# set up working directory
setwd("/path/to/my/directory")


#Read Constitutional Rate and Respiratory rate
Data_Neat <- read.csv("sim_new_data_1.csv", header = T)
Constitutional<-Data_Neat$xc_sim[1:1820]
Respiratory<-Data_Neat$xr_sim[1:1820]
CDC_daily_rates<-Data_Neat$CDC_daily[1:1820]
CDC_week_full=matrix(CDC_daily_rates,7,260)[1,]
dates<-as.Date(Data_Neat$Date[1:1820])

#Data frame for constitutional rates
df_c <- data.frame(date=dates, val=Constitutional)
#Data frame fro respiratory rates
df_r <- data.frame(date=dates, val=Respiratory)

#data frame for Constitutional rates 
Constitutional_data<-data.frame(
  day=dates,
  value=Constitutional)

#data frame for Respiratory rates 
Respiratory_data<-data.frame(
  day=dates,
  value=Respiratory)

#data frame for CDC daily rates 
CDC_daily_data<-data.frame(
  day=dates,
  value=CDC_daily_rates)


## set the missing date from 
head(CDC_daily_data)
#############################################################################
# BSTS Model for only Respiratory rates daily in 1820 days
niter_0=5000
ss_R<- AddLocalLinearTrend(list(), Respiratory)
ss_R <- AddSeasonal(ss_R, Respiratory, nseasons =2, season.duration = 180)  
model11 <- bsts(Respiratory, state.specification = ss_R,niter = niter_0,model.options = BstsOptions(save.state.contributions = TRUE),seed=1000)
summary(model11)
### Get a suggested number of burn-ins
burn1 <- SuggestBurn(0.1, model11)
pdf("Respiratory trend and sesonal effects 1.pdf", # name of pdf (need to include .pdf)
    width = 10, # width of resulting pdf in inches
    height = 10 # height of resulting pdf in inches
) 
plot(model11, "components")
dev.off()
##############################################################################
# BSTS Model for only Constutional rates daily
ss_C<- AddLocalLinearTrend(list(), Constitutional)
ss_C<- AddSeasonal(ss_C, Constitutional, nseasons =2, season.duration = 180)  
model22 <- bsts(Constitutional, state.specification = ss_C,niter = niter_0,model.options = BstsOptions(save.state.contributions = TRUE),seed=2000)


### Get a suggested number of burn-ins
burn2<- SuggestBurn(0.1, model22)
pdf("Constitutional trend and sesonal effects 1.pdf", # name of pdf (need to include .pdf)
    width = 10, # width of resulting pdf in inches
    height = 10 # height of resulting pdf in inches
) 
plot(model22, "components")
dev.off()
### Get a suggested number of burn-ins
burn=max(burn1,burn2)

summary1=summary(model11,burn = burn)
summary2=summary(model22,burn = burn)


sigma.obs_r_range=quantile(model11$sigma.obs[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.obs_r_range
sigma.level_r_range=quantile(model11$sigma.trend.level[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.level_r_range
sigma.slope_r_range=quantile(model11$sigma.trend.slope[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.slope_r_range
sigma.seasonal_r_range=quantile(model11$sigma.seasonal.2.180[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.seasonal_r_range


sigma.obs_c=quantile(model22$sigma.obs[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.obs_c
sigma.level_c=quantile(model22$sigma.trend.level[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.level_c
sigma.slope_c=quantile(model22$sigma.trend.slope[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.slope_c
sigma.seasonal_c=quantile(model22$sigma.seasonal.2.180[(burn+1):niter_0],probs=c(0.025,0.5,0.975))
sigma.seasonal_c

#############################################################################################################
################################################################################################################################
######################################full Downscaled Ili Daily #####################################################################################

# this is the downscaled loop number
downscaled_loop_number=niter_0-burn
cat("downscaled loop number")
downscaled_loop_number



# the initial values of sigma_w2 and sigma_z2
A=t(as.matrix(c(1,1,1,1,1,1,1)))*1/7
sigma_w2<-var(CDC_daily_rates)
sigma_z2<-1


myknots=read.csv("knots_sim_1.txt",sep = "",header=FALSE)
knots_dim=dim(myknots)
#set up the prior of theta, the variance is denoted as sigma_theta
sigma_theta=1000*diag(knots_dim[1]+3)
mu_theta=10*rep(1,knots_dim[1]+3)
thetar=rep(1,knots_dim[1]+3)

w_k=list()
w_k[[1]]=-1.15*10^{13}
xc_list=list()
xr_list=list()
h_loop_full=list()
#x1,x2,x3 are fixed coeffcicents vector
x1=matrix(0,downscaled_loop_number,1)
x2=matrix(0,downscaled_loop_number,1)
x3=matrix(0,downscaled_loop_number,1)
fixed_coef_full=list()
## downscaled ili loop is the loop of downscaled ili daily
downscaled_ili_loop_full=list()
spm_residual=list()

#This loop only downscaled the place other than the missing time
for (k in 1:downscaled_loop_number) {
  ### Extract the components for Constitutional rate 
  components_c3<- cbind.data.frame(
    model22$state.contributions[burn+k,"trend",],                               
    model22$state.contributions[burn+k,"seasonal.2.180",],dates)
  names(components_c3) <- c("C Trend", "C Seasonality","Date")  
  
  ### Extract the components for Respiratory rate
  components_R3 <- cbind.data.frame(
    model11$state.contributions[burn+k,"trend",],                               
    model11$state.contributions[burn+k,"seasonal.2.180",],dates)
  names(components_R3) <- c("R Trend", "R Seasonality","Date")
  #######################################################################
  #get h function by semipar
  # first define x1 (xc, stands for consitutional) and x2(xr,stands for respiratory)
  xc3<-components_c3$`C Trend`+components_c3$`C Seasonality`
  xr3<-components_R3$`R Trend`+components_R3$`R Seasonality`
  
  
  if (k==1){
    p=1
  } else {
    rate=min(1,exp(w_k[[k]]-w_k[[k-1]]))
    if (is.nan(rate)){p=1}
    else{
      p=rbinom(1,1,rate)
    }
  }
  
  if (p==1){
    xc=xc3
    xr=xr3
  } else {
    xc=xc
    xr=xr
  }
  
  xc_list[[k]]=xc
  xr_list[[k]]=xr
  
  # run the following code, if the knots are not fixed:
  # fitspm.z_new<-spm(CDC_daily_rates~f(xc,xr)), then save the knots for the first time
  fitspm.z_new<-spm(CDC_daily_rates~f(xc,xr,knots=myknots))

  
  
  fit_new<-fitspm.z_new$fit
  fittedvalue_new<-fitted.spm(fitspm.z_new)
  h_loop_full[[k]]=fittedvalue_new
  spm_residual[[k]]=fit_new$residuals
  number_of_weeks=length(fittedvalue_new)/7
  X_Z=as.matrix(fit_new$data[,2:(knots_dim[1]+4)])
  
  h_w=matrix(fittedvalue_new, nrow = 7, ncol=number_of_weeks)
  
  ILi_Daily_new=matrix(0,7,number_of_weeks)
  for (i in 1:number_of_weeks){
    Sigma_yz=solve(diag(7)/sigma_w2+t(A)%*%A/sigma_z2)
    Mu_yz=Sigma_yz%*%(diag(7)%*%as.matrix(h_w[,i])/sigma_w2+t(A)%*%CDC_week_full[i]/sigma_z2)
    ILi_Daily_new[,i]=mvrnorm(1,Mu_yz,Sigma_yz)}
  
  downscaled_ili_loop_full[[k]]= ILi_Daily_new
  
  vec_ILi_Daily=as.vector(downscaled_ili_loop_full[[k]])
  
  #get the scale for the truncated inverse gamma distribution of sigma_z2
  scale_sum=0
  for (i in 1:number_of_weeks){
    scale_sum=scale_sum+(CDC_week_full[i]-A%*% ILi_Daily_new[,i])^2
  }
  
  #simulate sigma_z2
  sigma_z2=1/rgammatr(1,number_of_weeks/2+1,2/scale_sum,range=c(0.01,10^9))
  
  # generate weights
  sigma_w2=1/rgammatr(1,1820/2+1,2/(t(vec_ILi_Daily-fittedvalue_new)%*%(vec_ILi_Daily-fittedvalue_new)),range=c(0.01,10^9))
  
  
  spm.coef=fit_new$coefficients
  fixed_coef=spm.coef$fixed
  fixed_coef_full[[k]]=fixed_coef
  x1[k]=fixed_coef[[1]]
  x2[k]=fixed_coef[[2]]
  x3[k]=fixed_coef[[3]]
  random_coef=spm.coef$random
  
  mu_theta=c(unname(fixed_coef),random_coef$dummy.group.vec.Handan)
  fit_aux=fitspm.z_new$aux
  sigma_theta=fit_aux$cov.mat
  
  Sigma_thetar=solve(t(X_Z)%*%X_Z/sigma_w2+solve(sigma_theta))
  Mu_thetar=Sigma_thetar%*%(t(X_Z)%*%vec_ILi_Daily/(sigma_w2)+solve(sigma_theta)%*%mu_theta)
  thetar=mvrnorm(1,Mu_thetar,Sigma_thetar)
  
  
  
  w_k[[k+1]]=dmvnorm(vec_ILi_Daily,X_Z%*%thetar,sigma_w2*diag(1820),log=TRUE)
}

# create a matrix according to the list of downscaled_ili_loop
downscaled_ili_loop_matrix_full=matrix(0,nrow=downscaled_loop_number,ncol=1820)
for (i in 1:downscaled_loop_number){
  downscaled_ili_loop_matrix_full[i,]=downscaled_ili_loop_full[[i]]
}

# create a matrix according to the list of h
h_loop_matrix_full=matrix(0,nrow=downscaled_loop_number,ncol=1820)
for (i in 1:downscaled_loop_number){
  h_loop_matrix_full[i,]=h_loop_full[[i]]
}


# #Find its 25 percentile and 97.5 percentile
downscaled_Ili_CI_full=matrix(0,1820,2)
downscaled_Ili_mid_full=matrix(0,1820,1)
for (i in 1:1820){
  downscaled_Ili_CI_full[i,]=quantile (as.vector(downscaled_ili_loop_matrix_full[,i]), probs=c(0.025,0.975))
  downscaled_Ili_mid_full[i,]=quantile (as.vector(downscaled_ili_loop_matrix_full[,i]), probs=0.5)
}


# Find the median of h
h_mid_full=matrix(0,1820,1)
for (i in 1:1820){
  h_mid_full[i,]=quantile (as.vector(h_loop_matrix_full[,i]), probs=0.5)
}

# use downscaled_ili_loop to calculate the error between estimated weekly CDC rate and true weekly CDC rate
estimated_CDC_weekly=list()
error_z=list()
sigma_error_z=rep(0,downscaled_loop_number)
for (i in 1:downscaled_loop_number){
  estimated_CDC_weekly[[i]]=unname(tapply(downscaled_ili_loop_full[[i]], (seq_along(downscaled_ili_loop_full[[i]])-1) %/% 7, sum))/7
  error_z[[i]]=estimated_CDC_weekly[[i]]-CDC_week_full
  sigma_error_z[i]=sd(error_z[[i]])
}
mean(sigma_error_z)
median(sigma_error_z)
quantile(sigma_error_z,0.025)
quantile(sigma_error_z,0.975)

# use spm_residual to calculate the standard deviation of error during spm
sigma_error_y=rep(0,downscaled_loop_number)
for (i in 1:downscaled_loop_number){
  sigma_error_y[i]=sd(spm_residual[[i]])
}
mean(sigma_error_y)
median(sigma_error_y)
quantile(sigma_error_y,0.025)
quantile(sigma_error_y,0.975)

#find the confidence interval for fixed coefficients
beta_0_range=quantile(x1,c(0.025,0.5,0.975))
beta_0_range

beta_1_range=quantile(x2,c(0.025,0.5,0.975))
beta_1_range

beta_2_range=quantile(x3,c(0.025,0.5,0.975))
beta_2_range



#make a data frame for h 
df.downscaled_Ili_h_full<-data.frame(
  day=dates,
  value= h_mid_full
)

downscaled_Ili_CI_fullup=downscaled_Ili_CI_full[,2]
downscaled_Ili_CI_fulllow=downscaled_Ili_CI_full[,1]

df.downscaled_Ili_CI_full<-data.frame(
  day=dates,
  value=downscaled_Ili_CI_fulllow,
  value2=downscaled_Ili_CI_fullup
)



df.downscaled_ILi_full<-data.frame(
  day=dates,
  value=downscaled_Ili_mid_full
)


#error for CDC daily between estimate and truth
error_estimate_truth=downscaled_Ili_mid_full-CDC_daily_rates

error_estimate_truth_range=quantile(error_estimate_truth,c(0.025,0.5,0.975))
error_estimate_truth_range
error_estimate_truth_mean=mean(error_estimate_truth)
error_estimate_truth_mean

error_estimate_truth_CDC_week=colSums(matrix(downscaled_Ili_mid_full,7,260))-colSums(matrix(CDC_daily_rates,7,260))
error_estimate_truth_CDC_week_range=quantile(error_estimate_truth_CDC_week,c(0.025,0.5,0.975))
error_estimate_truth_CDC_week_range
error_estimate_truth_CDC_week_mean=mean(error_estimate_truth_CDC_week)
error_estimate_truth_CDC_week_mean

# the error is the dowscaled ili estimate minus the true CDC daily rate
error_truth_estimate=downscaled_Ili_mid_full-CDC_daily_rates
MSE=mean(error_truth_estimate^2)
MSE_record=data.frame(number=1,MSE=MSE)
write.table(MSE_record, file = "MSE.csv",append = TRUE,sep=",",col.names = FALSE) 

write.table(downscaled_Ili_mid_full, file = "downscaled_Ili_mid_full.csv",append = TRUE,col.names = FALSE,sep=",") 

pdf("all rates 1820 with CI_new_sim_data_1.pdf", # name of pdf (need to include .pdf)
    width = 10, # width of resulting pdf in inches
    height = 10 # height of resulting pdf in inches
)    
ggplot()+
  geom_line(data=df.downscaled_ILi_full,aes(x=day,y=value,color="a"))+
  geom_ribbon(data=df.downscaled_Ili_CI_full,aes(x=day,ymin=value,ymax=value2),alpha=0.2)+
  geom_line(data=Constitutional_data,aes(x=day,y=value,color="b"))+
  geom_line(data=Respiratory_data,aes(x=day,y=value,color="c"))+
  geom_point(data=CDC_daily_data,aes(x=day,y=value,color="d"),size=0.1)+
  scale_color_manual(labels = c("Downscaled ILi Daily", "Cons Rate","Resp Rates","CDC", "h"),values = c("a"="black", "b"="seagreen4","c"="dodgerblue4","d"="red"))+
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none",text=element_text(size=18))+
  labs( y="log rates per million population", x = "Days")+
  scale_x_date(date_labels = "%Y %b",date_breaks  ="3 month")


dev.off()

save.image(file="new_sim_data_1.RData") 
