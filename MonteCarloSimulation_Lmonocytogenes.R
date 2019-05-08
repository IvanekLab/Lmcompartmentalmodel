#=====================================================================================================
# Title: Public health impact of foodborne exposure to naturally occuring virulence attenuated Listeria monocytogenes:
# inference from mouse and mathematical models
# Authors: Stout, A. E., Van Stelten, A., Heiden, J., Marquis, H., Ballou, M., Nightingale, K., Ivanek, R. 
# Correspondence: Alison Stout - aek68@cornell.edu 
# R version 3.5.0 (2018-04-23)
# Code last updated
#=====================================================================================================

setwd("C:/Users/aek68/Desktop/Listeria Priming/Journal of Royal Interface/Final Codes and output")

#Packages needed
library(deSolve)

#=====================================================================================================
#Model Parameters - VA(1) and FV(2)
#===========================================================================
US_population <- 300e6 #US Population 

Day <- 365 #Scaling parameter to convert from days to years
illnes.P_VA <- 1-exp((-10^(-13.76))*2.00*10^6) # 3.475602e-08 #Last term is the dose   
illnes.P_FV <- 1-exp((-10^(-10.44))*1.76*10^4) # 6.390172e-0  #Last term is the dose    
kappa_value_VA <- illnes.P_VA  
kappa_value_FV <- illnes.P_FV   
kappa_value_sc <- 0.0989 
kappa_value_pc <- 0.0992 
lambda_value_VA <- -log(1-10/365)  #0.02777956 - Rate of exposure to attenuated strain      
lambda_value_FV <- -log(1-12/365)  #0.0334293 - Rate of exposure to fully virulent strain 
alpha_value <- 1/14 #Colonized to protected
gamma_value <- 0.0359 #Loss of immunity 
births_value <- (1/(70*365))*US_population  #US births
deaths_value <- (1/(70*365))  #Natural deaths 
mu_value <-  0.0114    #Deaths rate from Lm
time_max <- 10000 # Days to run model for 

#====================================================================================================
# Model Initial Values 
#===========================================================================
S_init <- US_population
I_VA_init <- 0 #Ill VA
C_VA_init <- 0 #Colonized VA
I_FV_init <- 0  #Ill FV
C_FV_init <- 0  #Colonized FV
P_VA_init <- 0  #Protected
P_FV_init <- 0  #Protected 
D_VA_init <- 0 #Deaths from VA
D_FV_init <- 0 #Deaths from FV 
D_other_init <- 0 #Natural deaths 
I_VAyear_init <- 0 #Cumulative VA infections
I_FVyear_init <- 0 #Cumulative FV infections
C_VAyear_init <- 0 #Cumulative VA colonized
C_FVyear_init <- 0 #Cumulative FV colonized
P_VAyear_init <- 0 #Cumulative VA protected
P_FVyear_init <- 0 #Cumulative FV protected 

#=====================================================================================================
#Model Equations - 
#=====================================================================================================

basic_SP.dyn <- function(t,var,par) {
  
  # Parameters 
  lambda_1 <- par[1]     
  lambda_2 <- par[2] 
  kappa_VA <- par[3]  
  kappa_FV <- par[4] 
  kappa_sc <- par[5]  
  kappa_pc <- par[6] 
  alpha <- par[7] 
  gamma <- par[8] 
  births <- par[9] 
  death_rate <- par[10] 
  mu <- par[11] 
  
  #State values
  S <- var[1]
  I_VA <- var[2]
  C_VA <- var[3]
  I_FV <- var[4]
  C_FV <- var[5]
  P_VA <- var[6]
  P_FV <- var[7]
  D_VA <- var[8]
  D_FV <- var[9]
  D_other <- var[10]
  I_VAyear <- var[11]
  I_FVyear <- var[12]
  C_VAyear <- var[13]
  C_FVyear <- var[14]
  P_VAyear <- var[15]
  P_FVyear <- var[16]
  
  # Calculate the derivatives
  
  dS <- births + mu*(I_VA+I_FV) - (lambda_1*(kappa_VA+kappa_sc)+lambda_2*(kappa_FV+kappa_sc)+death_rate)*S + gamma*(P_VA + P_FV) 
  dI_VA <- (lambda_1*kappa_VA)*S - (alpha+mu+death_rate)*I_VA
  dC_VA <- lambda_1*((kappa_sc*S)+(kappa_pc*P_FV)+(kappa_pc*P_VA))-(alpha+death_rate)*C_VA
  dI_FV <- lambda_2*kappa_FV*S - (alpha+mu+death_rate)*I_FV
  dC_FV <- lambda_2*((kappa_sc*S)+(kappa_pc*P_VA)+(kappa_pc*P_FV)) - (alpha+death_rate)*C_FV
  dP_VA <- alpha*C_VA + alpha*I_VA - ((lambda_1*kappa_pc)+(lambda_2*kappa_pc)+gamma+death_rate)*P_VA
  dP_FV <- alpha*C_FV + alpha*I_FV - ((lambda_2*kappa_pc)+(lambda_1*kappa_pc)+gamma+death_rate)*P_FV
  dD_VA <- mu*I_VA
  dD_FV <- mu*I_FV
  dD_other <- death_rate*(S+I_VA+C_VA+I_FV+C_FV+P_VA+P_FV)
  dI_VAyear <- (lambda_1*kappa_VA)*S
  dI_FVyear <- (lambda_2*kappa_FV)*S
  dC_VAyear <- lambda_1*((kappa_sc*S)+(kappa_pc*P_FV)+(kappa_pc*P_VA))
  dC_FVyear <- lambda_2*((kappa_sc*S)+(kappa_pc*P_VA)+(kappa_pc*P_FV)) 
  dP_VAyear <- alpha*C_VA + alpha*I_VA 
  dP_FVyear <- alpha*C_FV + alpha*I_FV
  
  
  # Last instruction: return a list 
  
  return(list(c(dS, dI_VA, dC_VA, dI_FV, dC_FV, dP_VA, dP_FV, dD_VA, dD_FV, dD_other, dI_VAyear, dI_FVyear, dC_VAyear, dC_FVyear, dP_VAyear, dP_FVyear)))
}

#=====================================================================================================
# Run the model - 
#===========================================================================

basic_SP.sol=NULL
S=NULL
I_VA=NULL
C_VA=NULL
I_FV=NULL
C_FV=NULL
P_VA=NULL
P_FV=NULL
D_VA=NULL
D_FV=NULL
D_other=NULL
I_VAyear=NULL
I_FVyear=NULL
C_VAyear=NULL
C_FVyear=NULL
P_VAyear=NULL
P_FVyear=NULL
Total_population=NULL
Total_deaths=NULL
Total_illness=NULL
Total_colonized=NULL
Total_protected=NULL
I_all=NULL
C_year=NULL
P_Year=NULL

#=====================================================================================================
#Parameters to apply distributions and sample from (Monte Carlo)
#===========================================================================
n <- 10000 #Number of iterations

set.seed(145689)
kappa_value_sc <- numeric(n)
kappa_value_pc <- numeric(n)
I_all_1825 <- numeric(n) #Cases in year five 
gamma_value <- numeric(n) #Loss of immunity 
Total_colonized_1825 <- numeric(n)

Factor.x <- runif(n, min=0.5, max = 1.5)
kappa_value_sc <-  runif(n, min = 0, max = kappa_value_FV*10^5.5) 
kappa_value_pc <-  Factor.x * kappa_value_sc    
gamma_value <- runif(n, min=0.00006, max = 0.07) #Loss of immunity

upar<-cbind(kappa_value_sc, kappa_value_pc, gamma_value, Factor.x)

for(i in  1:n)
{
  basic_SP.par <- c(lambda_value_VA, lambda_value_FV, kappa_value_VA, kappa_value_FV, 
                    kappa_value_sc[i], kappa_value_pc[i], alpha_value, gamma_value[i], births_value, deaths_value, mu_value) 
  basic_SP.t <- seq(0, time_max, 1)
  basic_SP.init <- c(S_init, I_VA_init, C_VA_init, I_FV_init, C_FV_init, P_VA_init, P_FV_init, D_VA_init, D_FV_init, D_other_init, I_VAyear_init, I_FVyear_init,  C_VAyear_init, C_FVyear_init, P_VAyear_init, P_FVyear_init) 
  basic_SP.sol[[i]]  <- lsoda(basic_SP.init, basic_SP.t, basic_SP.dyn,basic_SP.par)
}

#=====================================================================================================
# Extract variables
#===========================================================================

for (i in 1:n) {
  
  I_VAyear[[i]] <- basic_SP.sol[[i]][,12] #Cumulative attenuated cases
  I_FVyear[[i]] <- basic_SP.sol[[i]][,13] #Cumulative fully virulent cases
  
#=====================================================================================================
# Calculated values
#===========================================================================
  
  I_all[[i]] <- I_VAyear[[i]] + I_FVyear[[i]]
  I_all_1825[[i]] <- I_all[[i]][1825] - I_all[[i]][1460]
  
}

data <- as.data.frame(cbind(I_all_1825, kappa_value_sc, gamma_value, Factor.x, kappa_value_pc))

View(data)

#=====================================================================================================
# Dataset of those iterations that are between 517 - 3161 
#===========================================================================
data <- data[ which(data$I_all_1825>557 & data$I_all<3161), ]
dim(data) 

#Revisualize
x11();par(mfrow=c(2,4))
hist(data$I_all_1825, main = "Cases")
hist(data$kappa_value_sc, main = "Kappa SC")
hist(data$kappa_value_pc, main = "Kappa PC")
hist(data$gamma_value, main = "Gamma")
boxplot(data$I_all_1825, main = "Cases")
boxplot(data$kappa_value_sc, main = "Kappa SC")
boxplot(data$kappa_value_pc, main = "Kappa PC")
boxplot(data$gamma_value, main = "Gamma")

median(data$I_all_1825)  
median(data$kappa_value_sc) 
median(data$kappa_value_pc) 
median(data$gamma_value) 

mean(data$I_all_1825)  
mean(data$kappa_value_sc) 
mean(data$kappa_value_pc) 
mean(data$gamma_value) 


#=====================================================================================================
# Dataset of those iterations that are between 1500-1700
#===========================================================================
data1 <- data[ which(data$I_all_1825>1500 & data$I_all<1700), ]
dim(data1) 

#Revisualize
x11();par(mfrow=c(2,4))
hist(data1$I_all_1825, main = "Cases")
hist(data1$kappa_value_sc, main = "Kappa SC")
hist(data1$kappa_value_pc, main = "Kappa PC")
hist(data1$gamma_value, main = "Gamma")
boxplot(data1$I_all_1825, main = "Cases")
boxplot(data1$kappa_value_sc, main = "Kappa SC")
boxplot(data1$kappa_value_pc, main = "Kappa PC")
boxplot(data1$gamma_value, main = "Gamma")

median(data1$I_all_1825)  
median(data1$kappa_value_sc) 
median(data1$kappa_value_pc) 
median(data1$gamma_value) 

mean(data1$I_all_1825)  
mean(data1$kappa_value_sc) 
mean(data1$kappa_value_pc) 
mean(data1$gamma_value) 




