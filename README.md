# dowscaling-ili
Bayesian Downscaling Ili project
####Dynamic downscaling and daily nowcasting from influenza surveillance data 
##############################
##########################################################################################################
#                               Dan Han, University of Louisiville, USA                         #######
#                                   Email: dan.han @louisville.edu
#                                         Version 1, March 10th 2022                              #######
##########################################################################################################
##########################################################################################################
# Purpose: 
The article [1] used publicly available data from the Centers for Disease Control and Prevention (CDC) and data obtained from the Michigan Department of Health and Human Services (MDHHS). 
The data from the MDHHS cannot be shared. This coding project uses the dataset sim_new_data_1.csv that contains simulated CDC weekly ili rates, constitutional daily rates and respiratory
# daily rates to estimate downscaled daily ili rates based on the model in the paper [1]. And calculate Mean squred error (MSE)
# between estimates and the true value.
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
# How to use the code (use "03102022 simulation1.R" as an example): 
1. The user needs to replace the dataset "sim_new_data_1.csv" with their own datasets.
2. If the user has presettings for knots in semiparametric regression part, replace the knots file.
   If the user has no presettings for knots in semiparametric regression part, aftet line 146, set k=1, and run the following code:

 k=1

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
  fitspm.z_new<-spm(CDC_daily_rates~f(xc,xr)), then save the knots for the first time

Then rerun the loop from line 148 to the end of the code in the first simulation example "03102022 simulation1.R"
