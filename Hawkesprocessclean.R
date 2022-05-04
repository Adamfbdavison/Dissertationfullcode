
######## Modelling COVID-19 using discrete-time Hawkes process ######## 
setwd("C:/Users/afbda/Desktop/Hawkes processes clean")
library(lme4)
library(ggplot2)
library(dplyr)
library(stringr)
library(zoo)
library(magrittr)
library(tidyverse)
library(tidyr)
library(reshape2)
library(gridExtra)
library(pracma)
library(zoo)
library(mvtnorm)
library(GGally)
library(xtable)
library(scales)
library(ggpubr)
# Re-initialise seed
set.seed(NULL)

#######################################Data prep ####################################
######################################################################################



Englanddeaths <- read.csv("Nationcase.csv")
Englanddeaths <- Englanddeaths[,c(2,4,6)]
Englanddeaths <- na.omit(Englanddeaths)
Englanddeaths <- Englanddeaths %>% 
  rename(
    Date = date,
    Country = areaName,
    n = newDeaths28DaysByDeathDate
  )

Englanddeaths <- Englanddeaths[,c(2,1,3)]
Englanddeaths <- Englanddeaths %>% map_df(rev)
Englanddeaths$n <- rollmean(Englanddeaths$n, 5, na.pad = TRUE , align = "center")
Englanddeaths <- na.omit(Englanddeaths)
Englanddeaths$n <- as.integer(Englanddeaths$n)
Englanddeaths$n <- as.numeric(Englanddeaths$n)
Englanddeaths$Date <- as.Date(Englanddeaths$Date)
Englanddeaths <- as.data.frame(Englanddeaths)

Englanddeaths[c(1:150),2] <- rep("England1",150)
Firstwaveup <- Englanddeaths[c(1:37),]
Firstwavedown <- Englanddeaths[c(38:150),]

Englanddeaths[c(195:400),2] <- rep("England2",206)
Secwaveup <- Englanddeaths[c(195:324),]
Secwavedown <- Englanddeaths[c(325:400),]

Englanddeaths[c(460:580),2] <- rep("England3",121)
Thiwaveup <- Englanddeaths[c(460:600),]


#################################functions###########################################
#######################################################################################



exp_trig = function(t, beta){
  beta * (1-beta)^(t-1)
}

calc_decay_vals = function(beta, decay_times, decay_times_counts, calc_grad=FALSE){
  
  # Calculations for decay function
  t = 1:max(unlist(decay_times))
  calc_decay = sapply(t, function(t) exp_trig(t, beta))
  decay_vals = lapply(decay_times, function(x) calc_decay[x])
  decay_vals_counts = Map('*',decay_vals,decay_times_counts)
  sum_decay = sapply(1:length(decay_vals_counts), function(x) sum(decay_vals_counts[[x]]))
  return(list(sum_decay))
}

lambda_cond_int = function(mu, alpha, beta, decay_times, decay_times_counts){
  
  sum_decay = calc_decay_vals(beta, decay_times, decay_times_counts)[[1]]
  lambda = mu + alpha * sum_decay
  
  return(lambda)
}

log_likelihood = function(country, mu_1, alpha_1, beta_1, mu_2, alpha_2, beta_2){
  
  eval_lambda_up = lambda_cond_int(mu_1, alpha_1, beta_1, decay_times_upward, decay_times_counts_up)
  eval_lambda_down = lambda_cond_int(mu_2, alpha_2, beta_2, decay_times_down, decay_times_counts_down)
  eval_lambda = c(eval_lambda_up, eval_lambda_down)
  event_counts = c(event_count_up[,2], event_count_down[,2])
  log_lik = sum(dpois(event_counts, lambda=eval_lambda, log=TRUE))
  return(log_lik)
}


loglikelihood = function(p0){
  u1 <- p0[1]
  a1 <- p0[2]
  b1 <- p0[3]
  u2 <- p0[4]
  a2 <- p0[5]
  b2 <- p0[6]
  
  eval_lambda_up = lambda_cond_int(u1, a1, b1, decay_times_upward, decay_times_counts_upward)
  eval_lambda_down = lambda_cond_int(u2, a2, b2, decay_times_downward, decay_times_counts_downward)
  eval_lambda = c(eval_lambda_up, eval_lambda_down)
  event_counts = c(event_count_upward[,2], event_count_downward[,2])
  deviation = (sqrt(10)/10)*(sqrt(event_counts))
  
  log_lik = -sum(dnorm(ceil(eval_lambda),mean  = event_counts, sd = deviation,  log=TRUE))
  return(log_lik)
}
# Function to filter data 
filter_data_cp_peak = function(data, country, cps){
  
  # Filter data to peak of interest
  modelling_data.df = data %>%
    filter(Country == country) 
  
  # Start observation period once there is 10 events 
  t1_ind = as.numeric(which(cumsum(modelling_data.df$n) >=10)[1])
  
  # Calculate number of days in sample
  start_date = modelling_data.df$Date[t1_ind]
  
  end_date = max(modelling_data.df$Date)
 
  max_T = as.numeric(difftime(end_date,  start_date, units="days")) + 1
  
  # Modelling data
  modelling_data.df %<>%
    filter(Date <= end_date) %>%
    mutate(time = as.numeric(difftime(Date, start_date, units="days") + 1))  %>%
    select(time, n) 
  
  # Calculate individual events and their previous events
  event_count = as.matrix(sapply(modelling_data.df, as.numeric))
  decay_times = lapply(1:max_T, function(t) as.numeric(t- event_count[event_count[,1] <t & event_count[,2] >0 ,1]))
  decay_times_counts = lapply(1:max_T, function(t) as.numeric(event_count[event_count [,1] <t & event_count[,2] >0 ,2]))
  
  # Create global variables
  assign("start_date", start_date, envir = .GlobalEnv)
  assign("event_count", event_count, envir = .GlobalEnv)
  
    
  # Upward trajectory    
  max_T_upward = cps[1]
  event_count_upward = event_count[t1_ind:(max_T_upward+t1_ind-1),]
  decay_times_upward = decay_times[1:max_T_upward]
  decay_times_counts_upward = decay_times_counts[1:max_T_upward]
    
  # Create global variables
  assign("event_count_upward", event_count_upward, envir = .GlobalEnv)
  assign("decay_times_upward", decay_times_upward, envir = .GlobalEnv)
  assign("decay_times_counts_upward", decay_times_counts_upward, envir = .GlobalEnv)
    
  # Downward trajectory   
    
  event_count_downward = event_count[(max_T_upward+t1_ind):(max_T+t1_ind-1),]
  decay_times_downward = decay_times[(max_T_upward+1):max_T]
  decay_times_counts_downward = decay_times_counts[(max_T_upward+1):max_T]
  assign("decay_times_downward", decay_times_downward, envir = .GlobalEnv)
  assign("decay_times_counts_downward", decay_times_counts_downward, envir = .GlobalEnv)
  assign("event_count_downward", event_count_downward, envir = .GlobalEnv)
      
    
  
  return()
}



########################################Nealder mead######################################################


###################################Wave 1 deaths ####################################
#####################################################################################


filter_data_cp_peak(Englanddeaths,"England1", 31)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



p0 <- c(10,10,0.1,10,10,0.1)
lower = c(0,0,0,0,0,0)
upper = c(100000,100,6,100000,100,6)


Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,10,2,20,10,2)))



Englamdaup1 <- lambda_cond_int(3.0973876 ,67.0352427,  4.1613673 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup1, type = "l")
plot(event_count_upward, type ="l")

Englamdadown1 <- lambda_cond_int(  0.3608873 ,63.8277412  ,4.2291585, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown1, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England1") 

Engfulllambda1 <- as.data.frame(c(Englamdaup1, Englamdadown1))
Engfulllambda1$Deaths <- modelling_data.df[c(6:150),3]
Engfulllambda1$Date <- (modelling_data.df[c(6:150),1])

par(mfrow = c(1,1))
p3 <- ggplot(data = Engfulllambda1, aes(y= Engfulllambda1[,1], x = Engfulllambda1[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda1[,2],x=Engfulllambda1[,3], colour = "Daily counts"), size = 1.4) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 1 Deaths")

p3

calc_decay = sapply(t, function(t) exp_trig(t, 4.2223))
decay_vals = lapply(decay_times_upward, function(x) calc_decay[x])
decay_vals_counts = Map('*',decay_vals,decay_times_counts_upward)
decay_vals_counts[[31]]

calc_decay = sapply(t, function(t) exp_trig(t, 2.27960733))
decay_vals = lapply(decay_times_downward, function(x) calc_decay[x])
decay_vals_counts = Map('*',decay_vals,decay_times_counts_downward)
decay_vals_counts[[31]]
#############################Wave 2 deaths ###########################################
#######################################################################################



filter_data_cp_peak(Englanddeaths,"England2", 127)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")


p0 <- c(10,20,0.1,10,20,0.1)
lower = c(0,0,0,0,0,0)
upper = c(100000,100000,5,100000,100000,5)


Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,10,2,20,10,2)))



Englamdaup2 <- lambda_cond_int( 0.03498618, 9.34338224, 2.31694503   , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2 <- lambda_cond_int(0.00000000, 7.88523028, 2.22620579, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England2") 

Engfulllambda2 <- as.data.frame(c(Englamdaup2, Englamdadown2))
Engfulllambda2$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda2$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda2, aes(y= Engfulllambda2[,1], x = Engfulllambda2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda2[,2],x=Engfulllambda2[,3]))


################################Wave 3 deaths #######################################
#######################################################################################



filter_data_cp_peak(Englanddeaths, "England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")


p0 <- c(10,10,0.1,10,10,0.1)
lower = c(0,0,0,0,0,0)
upper = c(100000,100,5,100000,100,6)


Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,10,2,20,10,2)))



Englamdaup3 <- lambda_cond_int(0.1519286 ,4.8385041, 1.7613092 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3 <- lambda_cond_int( 7.7268413 ,1.9213014 ,1.1497254, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England3") 

Engfulllambda3 <- as.data.frame(c(Englamdaup3, Englamdadown3))
Engfulllambda3$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3, aes(y= Engfulllambda3[,1], x = Engfulllambda3[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3[,2],x=Engfulllambda3[,3], colour = "Daily deaths")) + xlab("Date(2021)") + ylab("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths")
 

############################## PREPARE CASES DATA ##################################
######################################################################################



Englandcases <- read.csv("Nationcase.csv")
Englandcases <- Englandcases[,c(2,4,5)]
Englandcases <- na.omit(Englandcases)
Englandcases <- Englandcases %>% 
  rename(
    Date = date,
    Country = areaName,
    n = newCasesByPublishDate
  )

Englandcases <- Englandcases[,c(2,1,3)]
Englandcases <- Englandcases %>% map_df(rev)
Englandcases$n <- rollmean(Englandcases$n, 5, na.pad = TRUE , align = "center")
Englandcases <- na.omit(Englandcases)
Englandcases$n <- as.integer(Englandcases$n)
Englandcases$n <- as.numeric(Englandcases$n)
Englandcases$Date <- as.Date(Englandcases$Date)
Englandcases <- as.data.frame(Englandcases)

Englandcases[c(1:150),2] <- rep("England1",150)
Firstwaveup <- Englandcases[c(1:37),]
Firstwavedown <- Englandcases[c(38:150),]

Englandcases[c(200:420),2] <- rep("England2",221)
Secwaveup <- Englandcases[c(200:340),]
Secwavedown <- Englandcases[c(340:400),]

Englandcases[c(480:547),2] <- rep("England3",68)
Thiwaveup <- Englandcases[c(480:547),]


#####################WAVE 2 CASES (IN 2 SECTIONS) #####################################
#######################################################################################


filter_data_cp_peak(Englandcases,"England2", 142)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")


p0 <- c(300,20,1,300,20,1)
lower = c(0,0,0,0,0,0)
upper = c(100000,100000,6,100000,100000,6)


Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,10,2,20,10,2)))



Englamdaup2c <- lambda_cond_int(269.480452 , 53.672138  , 3.999651, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2c <- lambda_cond_int(350.666053 , 48.399510  , 3.968782, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England2") 

Engfulllambda2c <- as.data.frame(c(Englamdaup2c, Englamdadown2c))
Engfulllambda2c$Deaths <- modelling_data.df[c(1:221),3]
Engfulllambda2c$Date <- (modelling_data.df[c(1:221),1])

par(mfrow = c(1,1))
p4 <- ggplot(data = Engfulllambda2c, aes(y= Engfulllambda2c[,1], x = Engfulllambda2c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda2c[,2],x=Engfulllambda2c[,3], colour = "Daily count")) + xlab("Date (2020)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 2 Cases")

p4
#################################### Wave 3 ##########################################
#######################################################################################



filter_data_cp_peak(Englandcases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")


p0 <- c(1000,40,1,1000,40,1)
lower = c(0,0,0,0,0,0)
upper = c(100000,100000,100000,100000,100000,100000)


Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,10,2,20,10,2)))



Englamdaup3c <- lambda_cond_int(842.100685  , 16.030243   , 2.833294  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3c <- lambda_cond_int(1014.475505,  114.623388,    4.797966, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England3") 

Engfulllambda3c <- as.data.frame(c(Englamdaup3c, Englamdadown3c))
Engfulllambda3c$Deaths <- modelling_data.df[c(1:121),3]
Engfulllambda3c$Date <- (modelling_data.df[c(1:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3c, aes(y= Engfulllambda3c[,1], x = Engfulllambda3c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3c[,2],x=Engfulllambda3c[,3]))




#########################Recovering 10% of data #######################################
#######################################################################################



Randomints <- function(n) {return(as.vector(sort(sample(1:n, ceiling(n/10)))))}

Confidencepct <- function(fulllambda, amount, pct) {
  count <- 0
  n <- nrow(fulllambda)
  Random <- as.vector(sort(sample(1:n, ceiling(n*amount/100))))
  realvalues <- subset(event_count[,2], event_count[,1] >0)
  for (i in 1:length(Random)){
    if (fulllambda[Random[i],1] < realvalues[Random[i]]*(1+pct/100) & fulllambda[Random[i],1] > realvalues[Random[i]]*(1-pct/100)){
      count <- count + 1
    }
    else {
      count <- count + 0
    }
  }
  return(count/length(Random))
}

Confint <- function(fulllambda, pct) {
  realvalues <- subset(event_count[,2], event_count[,1] >0)
  count <- 0
  for (i in 1:nrow(fulllambda)){
    if (fulllambda[i,1] < realvalues[i]*(1+pct/100) & fulllambda[i,1] > realvalues[i]*(1-pct/100)){
      count <- count + 1
    }
    else {
      count <- count + 0
    }
  }
  return(count/nrow(fulllambda))
}


asdf <- subset(event_count, event_count[,1] >0)
asdf <- as.data.frame(asdf)
asdf$pred <- Engfulllambda1

filter_data_cp_peak(Englanddeaths,"England1", 31)
Confidencepct(Engfulllambda1, 10, 10)
Confidencepct(Engfulllambda1, 10, 100)
Confidencepct(Engfulllambda1, 10, 2)
Confidencepct(Engfulllambda1, 100, 5)
Confidencepct(Engfulllambda1, 100, 10)
Confidencepct(Engfulllambda1, 100, 15)
Confint(Engfulllambda1, 15)

filter_data_cp_peak(Englanddeaths,"England2", 127)
Confidencepct(Engfulllambda2, 10, 10)
Confidencepct(Engfulllambda2, 10, 100)
Confidencepct(Engfulllambda2, 10, 2)
Confidencepct(Engfulllambda2, 100, 5)
Confidencepct(Engfulllambda2, 100, 10)
Confidencepct(Engfulllambda2, 100, 15)

filter_data_cp_peak(Englanddeaths,"England3", 92)
Confidencepct(Engfulllambda3, 10, 10)
Confidencepct(Engfulllambda3, 10, 100)
Confidencepct(Engfulllambda3, 10, 2)
Confidencepct(Engfulllambda3, 100, 5)
Confidencepct(Engfulllambda3, 100, 10)
Confidencepct(Engfulllambda3, 100, 15)

filter_data_cp_peak(Englandcases,"England2", 92)
Confidencepct(Engfulllambda2c, 100, 5)
Confidencepct(Engfulllambda2c, 100, 10)
Confidencepct(Engfulllambda2c, 100, 15)

filter_data_cp_peak(Englandcases,"England3", 52)
Confidencepct(Engfulllambda3c, 100, 5)
Confidencepct(Engfulllambda3c, 100, 10)
Confidencepct(Engfulllambda3c, 100, 15)


##########################VACINATIONS###################################################
#######################################################################################



Extractvacs <- function(data,wave, cps, remove){
  vacwave = data %>%
    filter(Country == wave)
  assign("vacwaveup", vacwave[remove:(cps+remove-1),4], envir = .GlobalEnv)
  assign("vacwavedown", vacwave[(cps+remove):(nrow(vacwave)),4], envir = .GlobalEnv)
}


VAClambda_cond_int = function(mu, gam,  alpha, beta, decay_times, decay_times_counts, vacwave){
  
  sum_decay = calc_decay_vals(beta, decay_times, decay_times_counts)[[1]]
  lambda = mu + gam * vacwave + alpha * sum_decay
  
  return(lambda)
}

Vacloglikelihood = function(p0){
  u1 <- p0[1]
  g1 <- p0[2]
  a1 <- p0[3]
  b1 <- p0[4]
  u2 <- p0[5]
  g2 <- p0[6]
  a2 <- p0[7]
  b2 <- p0[8]
  
  eval_lambda_up = VAClambda_cond_int(u1,g1, a1, b1, decay_times_upward, decay_times_counts_upward, vacwaveup)
  eval_lambda_down = VAClambda_cond_int(u2,g2, a2, b2, decay_times_downward, decay_times_counts_downward, vacwavedown)
  eval_lambda = c(eval_lambda_up, eval_lambda_down)
  event_counts = c(event_count_upward[,2], event_count_downward[,2])
  deviation = (sqrt(10)/10)*(sqrt(event_counts))
  
  log_lik = -sum(dnorm(ceil(eval_lambda),mean  = event_counts, sd = deviation,  log=TRUE))
  return(log_lik)
}


cumvacdata <- read.csv("Nationvacscum.csv")
cumvacdata <- cumvacdata[,c(4,5)]
cumvacdata <- cumvacdata %>% 
  rename(
    Date = date,
    vac = cumVaccinationSecondDoseUptakeByVaccinationDatePercentage
  )
cumvacdata <- cumvacdata %>% map_df(rev)
cumvacdata$vac <- as.numeric(cumvacdata$vac)
cumvacdata$Date <- as.Date(cumvacdata$Date)
cumvacdata <- as.data.frame(cumvacdata)

Vacdatadeaths <- merge(Englanddeaths,cumvacdata, by ="Date", all = T)
Vacdatadeaths[is.na(Vacdatadeaths)] <- 0
Vacdatadeaths[604,4] <- as.numeric(70.3)


#####################wave 3 deaths + cum vacs outside the hawkes sum ############################################
#######################################################################################


filter_data_cp_peak(Vacdatadeaths,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatadeaths,"England3",92,3)



p0 <- c(20,0.5,20,1,20,0.5,20,1)
lower = c(0,-100,0,0,0,-100,0,0)
upper = c(100000,100000,100000,100000,100000,100000,100000,100000)


Nelder_Mead(Vacloglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(5,1,5,1,5,1,5,1)))



Englamdaup3vac <- VAClambda_cond_int(21.3086653, -0.4660829 ,17.9560235,  2.8511856 , decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3vac, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac <- VAClambda_cond_int(20.6358959, -0.5664049, 12.5772924,  2.4720429, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3vac, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeaths %>%
  filter(Country == "England3") 

Engfulllambda3vac <- as.data.frame(c(Englamdaup3vac, Englamdadown3vac))
Engfulllambda3vac$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac, aes(y= Engfulllambda3vac[,1], x = Engfulllambda3vac[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac[,2],x=Engfulllambda3vac[,3]))



###############################wave 3 cases + vacs outside the hawkes sum #################################
#######################################################################################



Vacdatacases <- merge(Englandcases,cumvacdata, by ="Date", all = T)
Vacdatacases[is.na(Vacdatacases)] <- 0


filter_data_cp_peak(Vacdatacases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatacases,"England3",52,1)


p0 <- c(1000,0,0.5,0.5,1000,0,0.5,0.5)
lower = c(0,-100,0,0,0,-100,0,0)
upper = c(100000,100000,100000,100000,100000,100000,100000,100000)

Nelder_Mead(Vacloglikelihood, p0, lower, upper, control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,1,5,1,20,1,5,1)))



Englamdaup3cvac <- VAClambda_cond_int(924.5205113,-6.1735809,   0.8333325,   0.5676053 , decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3cvac, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3cvac <- VAClambda_cond_int(762.9145163,  -5.0056789,   4.9828166 ,  1.8141849, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3cvac, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3cvac <- as.data.frame(c(Englamdaup3cvac, Englamdadown3cvac))
Engfulllambda3cvac$Deaths <- modelling_data.df[c(1:121),3]
Engfulllambda3cvac$Date <- (modelling_data.df[c(1:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3cvac, aes(y= Engfulllambda3cvac[,1], x = Engfulllambda3cvac[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3cvac[,2],x=Engfulllambda3cvac[,3]))


##################################### vaccinations within the hawkes sum deaths ##############
#######################################################################################




calc_decay_valsvacs = function(beta, decay_timesvacs, decay_times_countsvacs, calc_grad=FALSE){
  
  # Calculations for decay function
  t = 1:max(unlist(decay_timesvacs))
  calc_decay = sapply(t, function(t) exp_trig(t, beta))
  decay_valsvacs = lapply(decay_timesvacs, function(x) calc_decay[x])
  decay_vals_counts = Map('*',decay_valsvacs,decay_times_countsvacs)
  sum_decayvacs = sapply(1:length(decay_vals_counts), function(x) sum(decay_vals_counts[[x]]))
  return(list(sum_decayvacs))
}


Extractvacs <- function(data,wave, cps, remove){
  vacwave = data %>%
    filter(Country == wave)
  assign("vacwaveup", vacwave[remove:(cps+remove-1),4], envir = .GlobalEnv)
  assign("vacwavedown", vacwave[(cps+remove):(nrow(vacwave)),4], envir = .GlobalEnv)
}


VAC2lambda_cond_int = function(mu, alpha,  beta, gamma,theta, decay_times, decay_times_counts, decay_timesvac, decay_times_countsvac){
  
  sum_decay = calc_decay_vals(beta, decay_times, decay_times_counts)[[1]]
  sum_decayvacs = calc_decay_valsvacs(theta, decay_timesvac, decay_times_countsvac)[[1]]
  lambda = mu + alpha * sum_decay + gamma * sum_decayvacs
  
  return(lambda)
}


Vacloglikelihood2 = function(p0){
  u1 <- p0[1]
  a1 <- p0[2]
  b1 <- p0[3]
  g1 <- p0[4]
  t1 <- p0[5]
  u2 <- p0[6]
  a2 <- p0[7]
  b2 <- p0[8]
  g2 <- p0[9]
  t2 <- p0[10]
  
  
  eval_lambda_up = VAC2lambda_cond_int(u1,a1,b1,g1,t1, decay_times_upward, decay_times_counts_upward,decay_timesvac_upward, decay_times_countsvac_upward)
  eval_lambda_down = VAC2lambda_cond_int(u2,a2,b2,g2,t2, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
  eval_lambda = c(eval_lambda_up, eval_lambda_down)
  event_counts = c(event_count_upward[,2], event_count_downward[,2])
  deviation = (sqrt(10)/10)*(sqrt(event_counts))
  
  log_lik = -sum(dnorm(ceil(eval_lambda),mean  = event_counts, sd = deviation,  log=TRUE))
  return(log_lik)
}



# Function to filter data 
filter_data_cp_peakvac = function(data, wave, cps, start){
  
  # Filter data to peak of interest
  modelling_data.df = data %>%
    filter(Country == wave) 
  
  # Start observation period once there is 10 events 
  t1_ind = start
  
  # Calculate number of days in sample
  start_date = modelling_data.df$Date[t1_ind]
  
  end_date = max(modelling_data.df$Date)
  
  max_T = as.numeric(difftime(end_date,  start_date, units="days")) + 1
  
  # Modelling data
  modelling_data.df %<>%
    filter(Date <= end_date) %>%
    mutate(time = as.numeric(difftime(Date, start_date, units="days") + 1))  %>%
    select(time, vac) 
  
  # Calculate individual events and their previous events
  event_countvac = as.matrix(sapply(modelling_data.df, as.numeric))
  decay_timesvac = lapply(1:max_T, function(t) as.numeric(t- event_countvac[event_count[,1] <t & event_countvac[,2] >0 ,1]))
  decay_times_countsvac = lapply(1:max_T, function(t) as.numeric(event_countvac[event_countvac[,1] <t & event_countvac[,2] >0 ,2]))
  
  # Create global variables
  assign("start_datevac", start_date, envir = .GlobalEnv)
  assign("event_countvac", event_count, envir = .GlobalEnv)
  
  
  # Upward trajectory    
  max_T_upward = cps[1]
  event_countvac_upward = event_countvac[t1_ind:(max_T_upward+t1_ind-1),]
  decay_timesvac_upward = decay_timesvac[1:max_T_upward]
  decay_times_countsvac_upward = decay_times_countsvac[1:max_T_upward]
  
  # Create global variables
  assign("event_countvac_upward", event_countvac_upward, envir = .GlobalEnv)
  assign("decay_timesvac_upward", decay_timesvac_upward, envir = .GlobalEnv)
  assign("decay_times_countsvac_upward", decay_times_countsvac_upward, envir = .GlobalEnv)
  
  # Downward trajectory   
  
  event_countvac_downward = event_countvac[(max_T_upward+t1_ind):(max_T+t1_ind-1),]
  decay_timesvac_downward = decay_timesvac[(max_T_upward+1):max_T]
  decay_times_countsvac_downward = decay_times_countsvac[(max_T_upward+1):max_T]
  assign("decay_timesvac_downward", decay_timesvac_downward, envir = .GlobalEnv)
  assign("decay_times_countsvac_downward", decay_times_countsvac_downward, envir = .GlobalEnv)
  assign("event_countvac_downward", event_countvac_downward, envir = .GlobalEnv)
  
  
  
  return()
}




filter_data_cp_peak(Vacdatadeaths,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatadeaths,"England3", 92,3)


p0 <- c(10,10,2,0,1,10,10,2,0,1)
lower = c(0,0,0,-100,0,0,0,0,-100,0)
upper = c(100000,100000,100000,100000,5,100000,100000,100000,100000,5)

Nelder_Mead(Vacloglikelihood2, p0, lower, upper,  control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,5,1,5,1,20,5,1,5,1)))



Englamdaup3vac2 <- VAC2lambda_cond_int(1.5557057, 10.2051705,  2.3993980, -0.1402749,  1.4745250, decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3vac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac2 <- VAC2lambda_cond_int( 18.0544642, 10.6495358,  2.4702442, -1.8555912,  2.0377260, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3vac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeaths %>%
  filter(Country == "England3") 

Engfulllambda3vac2 <- as.data.frame(c(Englamdaup3vac2, Englamdadown3vac2))
Engfulllambda3vac2$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac2$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3]))

########################Vacs within the hawkes sum wave 3 cases ###########################################
#######################################################################################


filter_data_cp_peak(Vacdatacases,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatacases,"England3", 92,3)


p0 <- c(10,10,2,0,1,10,10,2,0,1)
lower = c(0,0,0,-1000,0,0,0,0,-1000,0)
upper = c(100000,100000,100000,100000,5,100000,100000,100000,100000,5)

Nelder_Mead(Vacloglikelihood2, p0, lower, upper,  control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,5,1,5,1,20,5,1,5,1)))



Englamdaup3cvac2 <- VAC2lambda_cond_int( 11.594525,  8.918878,  2.283492, -3.093477  ,1.443064 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3cvac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3cvac2 <- VAC2lambda_cond_int(15.960188 , 8.914304  ,2.293654 , 1.050023  ,1.456511, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3cvac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3cvac2 <- as.data.frame(c(Englamdaup3cvac2, Englamdadown3cvac2))
Engfulllambda3cvac2$Deaths <- modelling_data.df[c(1:121),3]
Engfulllambda3cvac2$Date <- (modelling_data.df[c(1:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3cvac2, aes(y= Engfulllambda3cvac2[,1], x = Engfulllambda3cvac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3cvac2[,2],x=Engfulllambda3cvac2[,3]))



########################Daily vac figures within the sum deaths###########################################
#######################################################################################


Nationdailyvac <- read.csv("Nationsvacdaily.csv")
Nationdailyvac <- Nationdailyvac[,c(4,5)]
Nationdailyvac <- Nationdailyvac %>% 
  rename(
    Date = date,
    vac = newPeopleVaccinatedSecondDoseByVaccinationDate
  )


Nationdailyvac <- Nationdailyvac %>% map_df(rev)
Nationdailyvac$vac <- as.numeric(Nationdailyvac$vac)
Nationdailyvac$Date <- as.Date(Nationdailyvac$Date)
Nationdailyvac <- as.data.frame(Nationdailyvac)

Vacdailydeaths <- merge(Englanddeaths,Nationdailyvac, by ="Date", all = T)
Vacdailydeaths[is.na(Vacdailydeaths)] <- 0


filter_data_cp_peak(Vacdailydeaths,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdailydeaths,"England3", 92,3)


p0 <- c(10,10,1,0,1,10,10,1,0,1)
lower = c(-100,-100,-100,-100,-100,-100,-100,-100,-100,-100)
upper = c(100000,100000,100000,100000,100000,100000,100000,100000,100000,100000)

Nelder_Mead(Vacloglikelihood2, p0, lower, upper,  control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,5,1,5,1,20,5,1,5,1)))




Englamdaup3vac2 <- VAC2lambda_cond_int(-2.7710015470,  4.5607285166,  1.6983950991,  0.0001171509,  2.1618840612 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3vac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac2 <- VAC2lambda_cond_int( -5.1691887280,  7.3579329416,  2.0332573768, -0.0007103961 , 1.7656079497, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3vac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdailydeaths %>%
  filter(Country == "England3") 

Engfulllambda3vac2 <- as.data.frame(c(Englamdaup3vac2, Englamdadown3vac2))
Engfulllambda3vac2$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac2$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3]))


#########################daily pct change within the hwakes sum deaths##############################################
#######################################################################################


Vacdatadeathspct <- Vacdatadeaths
Vacdatadeathspct <- Vacdatadeathspct %>%
  mutate(dailyvac = c(vac[1],diff(vac)))
Vacdatadeathspct <- Vacdatadeathspct[,-4]
Vacdatadeathspct <- Vacdatadeathspct %>% 
  rename(
    vac = dailyvac
  )



filter_data_cp_peak(Vacdatadeathspct,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatadeathspct,"England3", 92,3)


p0 <- c(10,2,2,-2,2,10,2,2,-2,2)
lower = c(0,0,0,-100,-100,0,0,0,-100,-100)
upper = c(100000,100000,100000,100000,100000,100000,100000,100000,100000,100000)

Nelder_Mead(Vacloglikelihood2, p0, lower, upper,  control = list(maxfun = 100000, FtolAbs = 1e-15, FtolRel = 1e-20, xst = c(20,5,1,5,1,20,5,1,5,1)))




Englamdaup3vac2 <- VAC2lambda_cond_int( 0.5252423 , 5.9227713 , 1.9351691 ,-1.5801092 , 1.8717156 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3vac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac2 <- VAC2lambda_cond_int( 16.1875007,  6.4356654,  2.1998044  ,0.9053236,  1.4559980, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3vac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeathspct %>%
  filter(Country == "England3") 

Engfulllambda3vac2 <- as.data.frame(c(Englamdaup3vac2, Englamdadown3vac2))
Engfulllambda3vac2$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac2$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3]))




###################################### R value ########################################
#######################################################################################


weeklyR <- read.csv("NationRvalue.csv")
weeklyR$date <- as.Date(weeklyR$date, format = "%d/%m/%Y")
dailyR <- data.frame(date = seq(as.Date("2020-05-29"), as.Date("2021-11-12"),"days"))
R = merge(dailyR,weeklyR, by = "date", all =TRUE)
R = na.locf(R)
R <- R %>% 
  rename(
    Date = date
  )
R[,1] <- as.character(R[,1])
R <- R[,c(1,7)]
CasesR <- merge(Englandcases, R, by = "Date", all = TRUE)






###############################nelder mead multiple runs #####################
##############################################################################
library(svMisc)
Neldermeadmulti <- function(a,b,lower,upper,n){
  datamat <- matrix(0, nrow=n, ncol = length(a)+1)
  paramavg <- rep(0, length(a))
  for(i in 1:n){
    a1 <- runif(1, a[1], b[1])
    b1 <- runif(1, a[2], b[2])
    c1 <- runif(1, a[3], b[3])
    a2 <- runif(1, a[4], b[4])
    b2 <- runif(1, a[5], b[5])
    c2 <- runif(1, a[6], b[6])
    p0 <- c(a1,b1,c1,a2,b2,c2)
    params <- Nelder_Mead(loglikelihood, p0, lower, upper, control = list(maxfun = 100000))

    datamat[i,1:6] <- params$par
    datamat[i,7] <- params$fval
    progress(i,n)
    Sys.sleep(0.02)
    if (i == n) message("Done!")
  }
  
  
  return(datamat)
}

library(svMisc)
Neldermeadmultivacsoutside <- function(a,b,lower,upper,n){
  datamat <- matrix(0, nrow=n, ncol = length(a)+1)
  paramavg <- rep(0, length(a))
  for(i in 1:n){
    a1 <- runif(1, a[1], b[1])
    b1 <- runif(1, a[2], b[2])
    c1 <- runif(1, a[3], b[3])
    d1 <- runif(1, a[4], b[4])
    a2 <- runif(1, a[5], b[5])
    b2 <- runif(1, a[6], b[6])
    c2 <- runif(1, a[7], b[7])
    d2 <- runif(1, a[8], b[8])
    p0 <- c(a1,b1,c1,d1,a2,b2,c2,d2)
    params <- Nelder_Mead(Vacloglikelihood, p0, lower, upper, control = list(maxfun = 100000))
    
    datamat[i,1:8] <- params$par
    datamat[i,9] <- params$fval
    progress(i,n)
    Sys.sleep(0.02)
    if (i == n) message("Done!")
  }
  
  
  return(datamat)
}



library(svMisc)
Neldermeadmultivacsoutside2 <- function(a,b,lower,upper,n){
  datamat <- matrix(0, nrow=n, ncol = length(a)+1)
  paramavg <- rep(0, length(a))
  for(i in 1:n){
    a1 <- runif(1, a[1], b[1])
    b1 <- runif(1, a[2], b[2])
    c1 <- runif(1, a[3], b[3])
    d1 <- runif(1, a[4], b[4])
    e1 <- runif(1, a[5], b[5])
    a2 <- runif(1, a[6], b[6])
    b2 <- runif(1, a[7], b[7])
    c2 <- runif(1, a[8], b[8])
    d2 <- runif(1, a[9], b[9])
    e2 <- runif(1, a[10], b[10])
    p0 <- c(a1,b1,c1,d1,e1,a2,b2,c2,d2,e1)
    params <- Nelder_Mead(Vacloglikelihood2, p0, lower, upper, control = list(maxfun = 100000))
    
    datamat[i,1:10] <- params$par
    datamat[i,11] <- params$fval
    progress(i,n)
    Sys.sleep(0.02)
    if (i == n) message("Done!")
  }
  
  
  return(datamat)
}























set.seed(1)
testa <- c(600,0.1,0.1,600,0.1,0.1)
testb <- c(1300,80,10,1300,80,10)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)


goodtest <- subset(test, test[,7] <1900 & test[,5] >0.1)
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(40,80))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(40,60))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(4,6))
legend(7000,2, legend=c("a1", "b1", "c1", "a2", "b2", "c2"),
       col=c("red", "blue", "green", "purple", "black", "brown"), lty=c(1,1,1,1,1,1), cex = 0.3)


write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave2casesparam.csv", row.names = FALSE)


wave1deaths <- read.csv("wave1param.csv")
wave2cases <- read.csv("wave2casesparam.csv")





##### plots ###


par(mfrow = c(1,2))
x1 <- seq(1,6,1)
x2 <- seq(0,6,0.001)
y1 <- dgeom(x1, 0.5)
y2 <- dexp(x2,2)

geommat <- matrix(0,7,6)
geommat[1,] <- dgeom(x1-1, 0.5)
geommat[2,] <- dgeom(x1-1, 0.6)
geommat[3,] <- dgeom(x1-1, 0.7)
geommat[4,] <- dgeom(x1-1, 0.8)
geommat[5,] <- dgeom(x1-1, 0.9)
geommat[6,] <- dgeom(x1-1, 1)
geommat[7,] <- x1
geommat <- as.data.frame(geommat)

rownames(geommat) = c("0.5", "0.6", "0.7", "0.8" , "0.9","1", "x")
geommat <- t(geommat)
geommat <- round(geommat,4)
geommat <- as.data.frame(geommat)
p2 <- ggplot(geommat, aes(y = geommat[,1] , x = geommat[,7])) + geom_line(size = 1.2, aes(colour = "0.5")) +
  geom_line(aes(y = geommat[,2] , x = geommat[,7], colour = "0.6"), size =1.2) +
  geom_line(aes(y = geommat[,3] , x = geommat[,7], colour = "0.7"), size = 1.2) +
  geom_line(aes(y = geommat[,4] , x = geommat[,7], colour = "0.8"), size = 1.2) +
  geom_line(aes(y = geommat[,5] , x = geommat[,7],colour = "0.9"), size = 1.2) +
  geom_line(aes(y = geommat[,6] , x = geommat[,7],colour = "1"), size = 1.2) +
  labs(color='p') + scale_x_continuous(breaks = seq(1,10,1)) + scale_y_continuous(breaks = seq(0,1.2,0.2))+
  ylab("Density") + xlab("Observed x (integer values)") + ggtitle("Geometric Distribution Densities")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )

expmat <- matrix(0,6,6001)
expmat[1,] <- dexp(x2, 0.5)
expmat[2,] <- dexp(x2, 1)
expmat[3,] <- dexp(x2, 1.5)
expmat[4,] <- dexp(x2, 2)
expmat[5,] <- dexp(x2, 2.5)
expmat[6,] <- x2
expmat <- as.data.frame(expmat)

rownames(expmat) = c("0.5", "1.0", "1.5", "2.0" , "2.5", "y")
expmat <- t(expmat)
expmat <- round(expmat,4)
expmat <- as.data.frame(expmat)
p1 <-ggplot(expmat, aes(y = expmat[,1] , x = expmat[,6])) + geom_line(size = 1.2, aes(colour = "0.5")) +
  geom_line(aes(y = expmat[,2] , x = expmat[,6], colour = "1.0"), size =1.2) +
  geom_line(aes(y = expmat[,3] , x = expmat[,6], colour = "1.5"), size = 1.2) +
  geom_line(aes(y = expmat[,4] , x = expmat[,6], colour = "2.0"), size = 1.2) + ylim(0,2.5) +
  geom_line(aes(y = expmat[,5] , x = expmat[,6],colour = "2.5"), size = 1.2) + xlim(1,10) +
  labs(color='Lambda') + scale_x_continuous(breaks = seq(1,10,1)) + scale_y_continuous(breaks = seq(0,2.5,0.5))+
  ylab("Density") + xlab("Observed x") + ggtitle("Exponential Distribution Densities") + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )

library(cowplot)
plot_grid(p1,p2)
plot_grid(p3, p4, ncol= 1, nrow = 2)



##wave plot


Englandcasesmodel <- subset(Englandcases, Country != "England")
Englandcasesextra <- subset(Englandcases, Country == "England")
Englandcasesextra <- Englandcasesextra %>% rename( blank = n)
Englandcasesmodel <- merge(Englandcasesmodel, Englandcasesextra, by= 'Date', all = TRUE)

p5 <- ggplot(Englandcasesmodel, aes(x=Date, y = n, colour = Country.x)) + geom_line(size = 1, show.legend = FALSE)+ ylab("Cases")+xlab("Date (2020-2021)") + ggtitle("Daily Cases and Waves in England") + scale_x_date(date_labels = "%b", date_breaks = "2 month") + geom_point(aes(x = Date, y = blank, colour = "No wave"), show.legend = FALSE, size = 0.9)+ scale_color_manual(name="Wave",   labels=c("1","2","2", "No wave"), values = c("red", "purple", "blue", "black")) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13),plot.title = element_text(size = 14)) 
                                                                                                                                                                                                                                                                                                                    
Englanddeathsmodel <- subset(Englanddeaths, Country != "England")
Englanddeathsextra <- subset(Englanddeaths, Country == "England")
Englanddeathsextra <- Englanddeathsextra %>% rename( blank = n)
Englanddeathsmodel <- merge(Englanddeathsmodel, Englanddeathsextra, by= 'Date', all = TRUE)

p6 <- ggplot(Englanddeathsmodel, aes(x=Date, y = n, colour = Country.x)) + geom_line(size = 1) +xlab("Date (2020-2021)")+ ylab("Deaths") + ggtitle("Daily Deaths and Waves in England") + scale_x_date(date_labels = "%b", date_breaks = "2 month") + geom_point(aes(x = Date, y = blank, colour = "No wave"), size = 0.9) + scale_color_manual(name="Wave",   labels=c("1","2","2", "No wave"), values = c("red", "purple", "blue", "black")) + theme(legend.position = c(0.87,0.86)) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) 

ggarrange(p5,p6,nrow = 1, ncol = 2, widths = c(0.5, 0.5))

seq <- seq(1,1000,1)
func <- function(n){
  return(exp(-2.335 * n))
}
w1d <- sapply(seq, func)
exp(-2.335)/sum(w1d)

par(mfrow =c(1,2))

unifsam <- runif(200,0,1)
normsam <- rnorm(200,0,1)


ks.test(unifsam,normsam)

plot(ecdf(unifsam), col = "red", main = "200 samples from U[0,1] and N(0,1) dists")
plot(ecdf(normsam), col = "blue", add = TRUE)


unifsam1 <- runif(200,0,1)
unifsam2 <- runif(200,0,1)
ks.test(unifsam1,unifsam2)

plot(ecdf(unifsam1), col = "green", main = "200 samples from two U[0,1] dists")
plot(ecdf(unifsam2), col = "purple", add = TRUE)






































#############wave 1 deaths normal dist 10000 runs ###############################################
############################################################################


filter_data_cp_peak(Englanddeaths,"England1", 31)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(20,100,10,20,100,10)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)

test <- read.csv("wave1deathsnormal.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.2), mgp=c(2,1,0), cex.lab = 1.15)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(0,3), font = 1)
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(50,80))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(3.5,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(0,2))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(40,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(3.5,5))


par(mfrow = c(1,1))
Englamdaup1 <- lambda_cond_int(1.6442030, 57.3288598 , 3.9567726 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup1, type = "l")
plot(event_count_upward, type ="l")

Englamdadown1 <- lambda_cond_int( 0.5109245, 57.6471535  ,4.1421826, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown1, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England1") 

Engfulllambda1 <- as.data.frame(c(Englamdaup1, Englamdadown1))
Engfulllambda1$Deaths <- modelling_data.df[c(6:150),3]
Engfulllambda1$Date <- (modelling_data.df[c(6:150),1])

par(mfrow = c(1,1))
p3 <- ggplot(data = Engfulllambda1, aes(y= Engfulllambda1[,1], x = Engfulllambda1[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda1[,2],x=Engfulllambda1[,3], colour = "Daily counts"), size = 1.4) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 1 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     date_labels = "%b")

p3



#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave1deathsnormal.csv", row.names = FALSE)


#############wave 2 deaths normal dist 10000 runs ###############################################
############################################################################

filter_data_cp_peak(Englanddeaths,"England2", 127)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(40,100,10,40,100,10)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("wave2deathsnormal.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(40,80))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup2 <- lambda_cond_int( 3.1362661, 57.9681531,  4.0825585  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2 <- lambda_cond_int(0.6731667, 57.7496492,  4.1475132, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England2") 

Engfulllambda2 <- as.data.frame(c(Englamdaup2, Englamdadown2))
Engfulllambda2$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda2$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda2, aes(y= Engfulllambda2[,1], x = Engfulllambda2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda2[,2],x=Engfulllambda2[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave2deathsnormal.csv", row.names = FALSE)


###############wave 3 deaths normal dist #############################################
####################################################################################
filter_data_cp_peak(Englanddeaths, "England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(20,50,6,20,40,4)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("wave3deathsnormal.csv")

goodtest <- subset(test,  test[,7] < quantile(test[,7], probs = c(.4)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,3))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(15,35))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(5,15))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(15,35))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(1,5))


Englamdaup3 <- lambda_cond_int(0.4840261, 26.9390765,  3.3372209, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3 <- lambda_cond_int( 8.9259013, 21.6896782,  3.2549696, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3, type = "l")
plot(event_count_downward, type ="l") 

wave3params <- c(0.4840261, 26.9390765,  3.3372209,8.9259013, 21.6896782,  3.2549696)
loglikelihood(wave3params)

modelling_data.df = Englanddeaths %>%
  filter(Country == "England3") 

Engfulllambda3 <- as.data.frame(c(Englamdaup3, Englamdadown3))
Engfulllambda3$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3, aes(y= Engfulllambda3[,1], x = Engfulllambda3[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3[,2],x=Engfulllambda3[,3], colour = "Daily deaths")) + xlab("Date(2021)") + ylab("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths")

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3deathsnormal.csv", row.names = FALSE)

######################################Wave2 cases 10000 normal dist ###########################

filter_data_cp_peak(Englandcases,"England2", 142)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(200,10,1,200,10,1)
testb <- c(1000,50,6,1000,50,6)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("wave2casesparam.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] >0.1)
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of mu 1', col="red" , ylim = c(200,300))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of alpha 1', ylim = c(30,50),  main = "Empirical Mean of Wave 2 Case Parameters")
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of beta 1', ylim = c(2,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of mu 2', ylim = c(350,600))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of alpha 2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of beta 2', ylim = c(2,5))


Englamdaup2c <- lambda_cond_int(209.061334,  35.460665  , 3.668630, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2c <- lambda_cond_int(654.235055 , 41.209664  , 3.914012, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England2") 

Engfulllambda2c <- as.data.frame(c(Englamdaup2c, Englamdadown2c))
Engfulllambda2c$Deaths <- modelling_data.df[c(1:221),3]
Engfulllambda2c$Date <- (modelling_data.df[c(1:221),1])

par(mfrow = c(1,1))
p4 <- ggplot(data = Engfulllambda2c, aes(y= Engfulllambda2c[,1], x = Engfulllambda2c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda2c[,2],x=Engfulllambda2c[,3], colour = "Daily count")) + xlab("Date (2020)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 2 Cases")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               date_labels = "%b")

p4

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave2casesparam2nd.csv", row.names = FALSE)


#########################WAVE 3 NORMAL DIST 10000 ITERATIONS###################################


filter_data_cp_peak(Englandcases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")





set.seed(1)
testa <- c(600,0.1,0.1,600,0.1,0.1)
testb <- c(1200,70,8,1200,70,8)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("wave3casesparam.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] >0.1)
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(700,950))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(30,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(2,4))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(800,1100))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3casesparam.csv", row.names = FALSE)



Englamdaup3c <- lambda_cond_int(910.311171 , 30.969810  , 3.546793  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3c <- lambda_cond_int(913.841861,  34.924824 ,  3.670757, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England3") 

Engfulllambda3c <- as.data.frame(c(Englamdaup3c, Englamdadown3c))
Engfulllambda3c$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3c$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3c, aes(y= Engfulllambda3c[,1], x = Engfulllambda3c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3c[,2],x=Engfulllambda3c[,3]))













#################################### Vacs outside the hawkes sum deaths ###########

filter_data_cp_peak(Vacdatadeaths,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatadeaths,"England3",92,3)



set.seed(1)
testa <- c(0.1,-5,0.1,0.1,0.1,-5,0.1,0.1)
testb <- c(20,5,50,6,20,5,40,4)
testc <- c(0,-100,0,0,0,-100,0,0)
testd <- c(10000,10000,10000,10000,10000,10000,10000,10000)


test <- Neldermeadmultivacsoutside(testa,testb,testc,testd,10000)
test <- read.csv("wave3deathsvacsoutside.csv")

goodtest <- subset(test, test[,9] < quantile(test[,9], probs = c(.8)))
goodtest[,c(3,7)] <- log(goodtest[,c(3,7)])
paramavg <- c(0,0,0,0,0,0,0,0)
for(i in 1:8){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(3,7)] <- exp(paramavg[c(3,7)])


Englamdaup3vac <- VAClambda_cond_int(2.38960589, -0.03956341, 24.67580987,  3.24323407, decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3vac, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac <- VAClambda_cond_int(34.86122836, -0.41150298, 17.57351104  ,3.04356307, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3vac, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeaths %>%
  filter(Country == "England3") 

Engfulllambda3vac <- as.data.frame(c(Englamdaup3vac, Englamdadown3vac))
Engfulllambda3vac$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac, aes(y= Engfulllambda3vac[,1], x = Engfulllambda3vac[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac[,2],x=Engfulllambda3vac[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3deathsvacsoutside.csv", row.names = FALSE)

Vaceffect <- subset(Vacdatadeaths[,c(1,4)], Vacdatadeaths[,2] == "England3")
Vaceffect[c(1:92),2] <- Vaceffect[c(1:92),2] * -0.03956341
Vaceffect[c(93:121),2] <- Vaceffect[c(93:121),2] * -0.41150298
Vaceffect$pctvac <- subset(Vacdatadeaths[,c(4)], Vacdatadeaths[,2] == "England3")
ggplot(Vaceffect, aes(x= pctvac, y = vac)) + geom_line() +ylim(-30,10)
ggplot(Vaceffect, aes(x= Date, y = vac)) + geom_line() +ylim(-30,10)


#################################### Vacs outside the hawkes sum cases ###########
Vacdatacases <- merge(Englandcases,cumvacdata, by ="Date", all = T)
Vacdatacases[is.na(Vacdatacases)] <- 0


filter_data_cp_peak(Vacdatacases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatacases,"England3",52,1)



set.seed(1)
testa <- c(600,-30,0.1,0.1,600,-30,0.1,0.1)
testb <- c(1200,30,70,8,1200,30,70,8)
testc <- c(0,-1000,0,0,0,-1000,0,0)
testd <- c(10000,10000,10000,10000,10000,10000,10000,10000)

test <- Neldermeadmultivacsoutside(testa,testb,testc,testd,10000)
test <- read.csv("wave3casesvacsoutside.csv")

goodtest <- subset(test, test[,9] < quantile(test[,9], probs = c(.8)))
goodtest[,c(3,7)] <- log(goodtest[,c(3,7)])
paramavg <- c(0,0,0,0,0,0,0,0)
for(i in 1:8){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(3,7)] <- exp(paramavg[c(3,7)])


Englamdaup3vac <- VAClambda_cond_int(1114.514135 , -11.242597,    6.144356,    2.132725, decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3vac, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac <- VAClambda_cond_int(986.081719 ,  -1.303130  , 63.576978 ,   4.257654, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3vac, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3vac <- as.data.frame(c(Englamdaup3vac, Englamdadown3vac))
Engfulllambda3vac$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3vac$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac, aes(y= Engfulllambda3vac[,1], x = Engfulllambda3vac[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac[,2],x=Engfulllambda3vac[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3casesvacsoutside.csv", row.names = FALSE)

Vaceffect <- subset(Vacdatadeaths[,c(1,4)], Vacdatadeaths[,2] == "England3")
Vaceffect[c(1:92),2] <- Vaceffect[c(1:92),2] * -0.03956341
Vaceffect[c(93:121),2] <- Vaceffect[c(93:121),2] * -0.41150298
Vaceffect$pctvac <- subset(Vacdatadeaths[,c(4)], Vacdatadeaths[,2] == "England3")
ggplot(Vaceffect, aes(x= pctvac, y = vac)) + geom_line() +ylim(-30,10)
ggplot(Vaceffect, aes(x= Date, y = vac)) + geom_line() +ylim(-30,10)








#################DAILY VAC deaths WITHIN THE HAWKES SUM ##################


filter_data_cp_peak(Vacdatadeathspct,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatadeathspct,"England3", 92,3)


set.seed(1)
testa <- c(0.1,0.1,0.1,-50,0.1,0.1,0.1,0.1,-50,0.1)
testb <- c(20,50,6,20,6,20,40,4,20,6)
testc <- c(0,0,0,-1000,0,0,0,0,-1000,0)
testd <- c(10000,10000,10000,10000,10000,10000,10000,10000,10000,10000)

test <- Neldermeadmultivacsoutside2(testa,testb,testc,testd,10000)
#test <- read.csv("wave3casesvacsoutside.csv")


goodtest <- subset(test, test[,11] < quantile(test[,11], probs = c(.8)))
goodtest[,c(2,4,7,9)] <- log(goodtest[,c(2,4,7,9)])
paramavg <- c(0,0,0,0,0,0,0,0,0,0)
for(i in 1:10){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,4,7,9)] <- exp(paramavg[c(2,4,7,9)])
paramavg[c(4,9)] <- -paramavg[c(4,9)]
paramavg




Englamdaup3vac2 <- VAC2lambda_cond_int( 0.6754597,  28.7949228,   3.3931289, -13.3496547,   8.4476729 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3vac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac2 <- VAC2lambda_cond_int(  12.7107266 , 21.0196055  , 3.2613814, -12.8120997,  10.1187683, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3vac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeathspct %>%
  filter(Country == "England3") 

Engfulllambda3vac2 <- as.data.frame(c(Englamdaup3vac2, Englamdadown3vac2))
Engfulllambda3vac2$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac2$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3]))

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3deathsvacinside.csv", row.names = FALSE)




#####Daily vac cases within the hawkes sum ##############

filter_data_cp_peak(Vacdatacases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatacases,"England3", 52,1)


set.seed(1)
testa <- c(200,20,1,-10,0.1,200,20,1,-10,0.1)
testb <- c(2000,40,5,10,8,2000,40,5,10,8)
testc <- c(0,0,0,-10000,0,0,0,0,-10000,0)
testd <- c(10000,10000,10000,10000,10000,10000,10000,10000,10000,10000)

test <- Neldermeadmultivacsoutside2(testa,testb,testc,testd,10000)
#test <- read.csv("wave3casesvacsoutside.csv")

goodtest <- subset(test, test[,11] < quantile(test[,11], probs = c(.8)) & test[,4] < 0 & test[,9] < 0)
goodtest[,c(4,9)] <- -goodtest[,c(4,9)]
goodtest[,c(2,4,7,9)] <- log(goodtest[,c(2,4,7,9)])
paramavg <- c(0,0,0,0,0,0,0,0,0,0)
for(i in 1:10){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,4,7,9)] <- exp(paramavg[c(2,4,7,9)])
paramavg[c(4,9)] <- -paramavg[c(4,9)]
paramavg



Englamdaup3cvac2 <- VAC2lambda_cond_int(919.034245 ,  28.177983 ,   3.364140  , - 6.521545 ,   4.141068 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3cvac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3cvac2 <- VAC2lambda_cond_int(1112.172499,   31.333024  ,  3.559333 ,   -3.872889 ,   6.444376, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3cvac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3cvac2 <- as.data.frame(c(Englamdaup3cvac2, Englamdadown3cvac2))
Engfulllambda3cvac2$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3cvac2$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3cvac2, aes(y= Engfulllambda3cvac2[,1], x = Engfulllambda3cvac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3cvac2[,2],x=Engfulllambda3cvac2[,3]))

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/wave3casesvacinside.csv", row.names = FALSE)



######################################## BY age ###############################################################


Agedeaths <- read.csv("Deathsage.csv")
Agedeaths <- Agedeaths[,c(2,4,5,6)]
Agedeaths <- na.omit(Agedeaths)
Agedeaths <- Agedeaths %>% 
  rename(
    Date = date,
    Country = areaName,
    n = deaths,
    age = age
  )

Agedeaths <- Agedeaths[,c(2,1,4,3)]
Agedeaths <- Agedeaths %>% map_df(rev)
#Englanddeaths$n <- rollmean(Englanddeaths$n, 5, na.pad = TRUE , align = "center")
#Englanddeaths <- na.omit(Englanddeaths)
Agedeaths$n <- as.integer(Agedeaths$n)
Agedeaths$n <- as.numeric(Agedeaths$n)
Agedeaths$Date <- as.Date(Agedeaths$Date)
Agedeaths <- as.data.frame(Agedeaths)

Agedeaths <- subset(Agedeaths, !(Agedeaths$age %in% c("60+", "00_04", "00_59", "05_09", "10_14", "15_19","20_24")))

Country <- rep(0,605)
Country[1:150] <- rep("England1",150)
Country[195:400] <- rep("England2",206)
Country[460:580] <- rep("England3",121)
Country[151:194] <- rep("England",44)
Country[401:459] <- rep("England",59)
Country[580:605] <- rep("England",26)

Agedeaths20 <- subset(Agedeaths, (Agedeaths$age %in% c("25_29", "30_34", "35_39","40_44")))
Agedeaths20 <- aggregate(Agedeaths20[,c(3)], FUN="sum", by=list(Agedeaths20$Date))
Agedeaths20$Country <- Country
Agedeaths20 <- Agedeaths20 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths20 <- Agedeaths20[,c(1,3,2)]

Agedeaths45 <- subset(Agedeaths, (Agedeaths$age %in% c("45_49", "50_54", "55_59", "60_64", "65_69")))
Agedeaths45 <- aggregate(Agedeaths45[,c(3)], FUN="sum", by=list(Agedeaths45$Date))
Agedeaths45$Country <- Country
Agedeaths45 <- Agedeaths45 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths45 <- Agedeaths45[,c(1,3,2)]

Agedeaths70 <- subset(Agedeaths, (Agedeaths$age %in% c("70_74", "75_79", "80_84", "85_89", "89+")))
Agedeaths70 <- aggregate(Agedeaths70[,c(3)], FUN="sum", by=list(Agedeaths70$Date))
Agedeaths70$Country <- Country
Agedeaths70 <- Agedeaths70 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths70 <- Agedeaths70[,c(1,3,2)]

Agedeaths20$n <- rollmean(Agedeaths20$n, 5, na.pad = TRUE , align = "center")
Agedeaths20 <- na.omit(Agedeaths20)
Agedeaths20$n <- ceil(Agedeaths20$n)
Agedeaths45$n <- rollmean(Agedeaths45$n, 5, na.pad = TRUE , align = "center")
Agedeaths45 <- na.omit(Agedeaths45)
Agedeaths45$n <- ceil(Agedeaths45$n)
Agedeaths70$n <- rollmean(Agedeaths70$n, 5, na.pad = TRUE , align = "center")
Agedeaths70 <- na.omit(Agedeaths70)
Agedeaths70$n <- ceil(Agedeaths70$n)

###Visalisations###
Agedatadeaths <- merge(Agedeaths20, Agedeaths45 , by =c("Date", "Country"))
Agedatadeaths <- merge(Agedatadeaths, Agedeaths70 , by =c("Date", "Country"))
Agedatadeaths <- Agedatadeaths %>% 
  rename("20-44" = n.x,
         "45-69" = n.y,
         "70+" = n
  )
ggplot(Agedatadeaths, aes(x= Date, Agedatadeaths$`20-44`, colour = "20-44")) +geom_line() +geom_line(aes(x=Date,Agedatadeaths$`45-69`, colour = "45-69"))+geom_line(aes(x=Date,y=Agedatadeaths$`70+`, , colour = "70+"))


Agecases <- read.csv("Casesage.csv")
Agecases <- Agecases[,c(2,4,5,6)]
Agecases <- na.omit(Agecases)
Agecases <- Agecases %>% 
  rename(
    Date = date,
    Country = areaName,
    n = cases,
    age = age
  )

Agecases <- Agecases[,c(2,1,4,3)]
Agecases <- Agecases %>% map_df(rev)
#Englanddeaths$n <- rollmean(Englanddeaths$n, 5, na.pad = TRUE , align = "center")
#Englanddeaths <- na.omit(Englanddeaths)
Agecases$n <- as.integer(Agecases$n)
Agecases$n <- as.numeric(Agecases$n)
Agecases$Date <- as.Date(Agecases$Date)
Agecases <- as.data.frame(Agecases)

Agecases <- subset(Agecases, !(Agecases$age %in% c("60+", "00_04", "00_59", "05_09", "10_14", "15_19","20_24")))
Agecases <- subset(Agecases, Agecases$Date > "2020-03-01")

Country <- rep(0,605)
Country[1:150] <- rep("England1",150)
Country[195:400] <- rep("England2",206)
Country[460:580] <- rep("England3",121)
Country[151:194] <- rep("England",44)
Country[401:459] <- rep("England",59)
Country[580:605] <- rep("England",26)

Agecases20 <- subset(Agecases, (Agecases$age %in% c("25_29", "30_34", "35_39","40_44")))
Agecases20 <- aggregate(Agecases20[,c(3)], FUN="sum", by=list(Agecases20$Date))
Agecases20$Country <- Country
Agecases20 <- Agecases20 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases20 <- Agecases20[,c(1,3,2)]

Agecases45 <- subset(Agecases, (Agecases$age %in% c("45_49", "50_54", "55_59", "60_64", "65_69")))
Agecases45 <- aggregate(Agecases45[,c(3)], FUN="sum", by=list(Agecases45$Date))
Agecases45$Country <- Country
Agecases45 <- Agecases45 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases45 <- Agecases45[,c(1,3,2)]

Agecases70 <- subset(Agecases, (Agecases$age %in% c("70_74", "75_79", "80_84", "85_89", "89+")))
Agecases70 <- aggregate(Agecases70[,c(3)], FUN="sum", by=list(Agecases70$Date))
Agecases70$Country <- Country
Agecases70 <- Agecases70 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases70 <- Agecases70[,c(1,3,2)]

Agecases20$n <- rollmean(Agecases20$n, 5, na.pad = TRUE , align = "center")
Agecases20 <- na.omit(Agecases20)
Agecases45$n <- rollmean(Agecases45$n, 5, na.pad = TRUE , align = "center")
Agecases45 <- na.omit(Agecases45)
Agecases70$n <- rollmean(Agecases70$n, 5, na.pad = TRUE , align = "center")
Agecases70 <- na.omit(Agecases70)

###Visalisations###
Agedatacases <- merge(Agecases20, Agecases45 , by =c("Date", "Country"))
Agedatacases <- merge(Agedatacases, Agecases70 , by =c("Date", "Country"))
Agedatacases <- Agedatacases %>% 
  rename("20-44" = n.x,
         "45-69" = n.y,
         "70+" = n
  )
ggplot(Agedatacases, aes(x= Date, Agedatacases$`20-44`, colour = "20-44")) +geom_line() +geom_line(aes(x=Date,Agedatacases$`45-69`, colour = "45-69"))+geom_line(aes(x=Date,y=Agedatacases$`70+`, , colour = "70+"))


##########age 45 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths20,"England1", 29)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(2,60,5,2,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
#test <- read.csv("wave2deathsnormal.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup120 <- lambda_cond_int(0.2509354,
                                 30.4917179,
                                 4.1188409
                                  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup120, type = "l")
plot(event_count_upward, type ="l")

Englamdadown120 <- lambda_cond_int(
                                   1.2369513,
                                   25.5067962,
                                   2.8406203, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown120, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England1") 

Engfulllambda120 <- as.data.frame(c(Englamdaup120, Englamdadown120))
Engfulllambda120$Deaths <- modelling_data.df[c(1:136),3]
Engfulllambda120$Date <- (modelling_data.df[c(1:136),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda120, aes(y= Engfulllambda120[,1], x = Engfulllambda120[,3])) +geom_line() +geom_point(aes(y=Engfulllambda120[,2],x=Engfulllambda120[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave1death.csv", row.names = FALSE)




#####Age 20 wave 2 deaths #####


filter_data_cp_peak(Agedeaths20,"England2", 120)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(5,60,5,5,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age20wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220 <- lambda_cond_int(0.042021518 , 7.134938588,  2.158705238 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220 <- lambda_cond_int(0.008681788, 57.970346388,  4.182533744, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England2") 

Engfulllambda220 <- as.data.frame(c(Englamdaup220, Englamdadown220))
Engfulllambda220$Deaths <- modelling_data.df[c(11:206),3]
Engfulllambda220$Date <- (modelling_data.df[c(11:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda220, aes(y= Engfulllambda220[,1], x = Engfulllambda220[,3])) +geom_line() +geom_point(aes(y=Engfulllambda220[,2],x=Engfulllambda220[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave2death.csv", row.names = FALSE)



#####Age 20 wave 3 deaths #####


filter_data_cp_peak(Agedeaths20, "England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(5,60,5,5,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age20wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220 <- lambda_cond_int(0.042021518,  7.134938588,  2.158705238 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220 <- lambda_cond_int( 0.008681788, 57.970346388 , 4.182533744, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England2") 

Engfulllambda220 <- as.data.frame(c(Englamdaup220, Englamdadown220))
Engfulllambda220$Deaths <- modelling_data.df[c(10:120),3]
Engfulllambda220$Date <- (modelling_data.df[c(10:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda220, aes(y= Engfulllambda220[,1], x = Engfulllambda220[,3])) +geom_line() +geom_point(aes(y=Engfulllambda220[,2],x=Engfulllambda220[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave2death.csv", row.names = FALSE)


##########age 45 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths45,"England1", 29)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(10,60,5,20,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age45wave1death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] > 0 )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup145 <- lambda_cond_int(0.8906642, 32.7470112  ,3.4665409  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup145, type = "l")
plot(event_count_upward, type ="l")

Englamdadown145 <- lambda_cond_int(0.1667018, 36.3700139,  3.7269861, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown145, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England1") 

Engfulllambda145 <- as.data.frame(c(Englamdaup145, Englamdadown145))
Engfulllambda145$Deaths <- modelling_data.df[c(7:147),3]
Engfulllambda145$Date <- (modelling_data.df[c(7:147),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda145, aes(y= Engfulllambda145[,1], x = Engfulllambda145[,3])) +geom_line() +geom_point(aes(y=Engfulllambda145[,2],x=Engfulllambda145[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age45wave1death.csv", row.names = FALSE)



#####Age 45 wave 2 deaths #####


filter_data_cp_peak(Agedeaths45,"England2", 126)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(10,60,5,10,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age45wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup245 <- lambda_cond_int( 0.28908511, 35.35506749,  3.60984843  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup245, type = "l")
plot(event_count_upward, type ="l")

Englamdadown245 <- lambda_cond_int( 0.09481242 ,34.58332326 , 3.64890019, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown245, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England2") 

Engfulllambda245 <- as.data.frame(c(Englamdaup245, Englamdadown245))
Engfulllambda245$Deaths <- modelling_data.df[c(5:206),3]
Engfulllambda245$Date <- (modelling_data.df[c(5:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda245, aes(y= Engfulllambda245[,1], x = Engfulllambda245[,3])) +geom_line() +geom_point(aes(y=Engfulllambda245[,2],x=Engfulllambda245[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age45wave2death.csv", row.names = FALSE)


#####Age 45 wave 3 deaths #####


filter_data_cp_peak(Agedeaths45, "England3", 104)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,30,0.1)
testb <- c(10,60,5,20,40,3)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,4)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age45wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] > 0 )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup345 <- lambda_cond_int(0.03252644 ,32.45149606,  3.54732262 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup345, type = "l")
plot(event_count_upward, type ="l")

Englamdadown345 <- lambda_cond_int( 8.11204311, 27.36934475  ,3.84730038, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown345, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England3") 

Engfulllambda345 <- as.data.frame(c(Englamdaup345, Englamdadown345))
Engfulllambda345$Deaths <- modelling_data.df[c(4:119),3]
Engfulllambda345$Date <- (modelling_data.df[c(4:119),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda345, aes(y= Engfulllambda345[,1], x = Engfulllambda345[,3])) +geom_line() +geom_point(aes(y=Engfulllambda345[,2],x=Engfulllambda345[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age45wave3death.csv", row.names = FALSE)

##########age 70 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths70,"England1", 32)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(10,60,5,20,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age70wave1death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] > 0 )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup170 <- lambda_cond_int( 0.5713820 ,37.3105564,  3.5389946 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup170, type = "l")
plot(event_count_upward, type ="l")

Englamdadown170 <- lambda_cond_int(0.1270915, 35.3369525 , 3.6668892, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown170, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England1") 

Engfulllambda170 <- as.data.frame(c(Englamdaup170, Englamdadown170))
Engfulllambda170$Deaths <- modelling_data.df[c(4:146),3]
Engfulllambda170$Date <- (modelling_data.df[c(4:146),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda170, aes(y= Engfulllambda170[,1], x = Engfulllambda170[,3])) +geom_line() +geom_point(aes(y=Engfulllambda170[,2],x=Engfulllambda170[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age70wave1death.csv", row.names = FALSE)

##########age 70 + wave 2 death #############################################


filter_data_cp_peak(Agedeaths70,"England2", 129)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(10,60,5,20,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age70wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] > 0 )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup270 <- lambda_cond_int(1.6572342, 30.9496942  ,3.4734216 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup270, type = "l")
plot(event_count_upward, type ="l")

Englamdadown270 <- lambda_cond_int(0.3226746 ,30.4335420,  3.5307157, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown270, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England2") 

Engfulllambda270 <- as.data.frame(c(Englamdaup270, Englamdadown270))
Engfulllambda270$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda270$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda270, aes(y= Engfulllambda270[,1], x = Engfulllambda270[,3])) +geom_line() +geom_point(aes(y=Engfulllambda270[,2],x=Engfulllambda270[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age70wave2death.csv", row.names = FALSE)

###################################AGE 70 WAVE 3 DEATHS ####################

filter_data_cp_peak(Agedeaths70, "England3", 93)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,10,30,0.1)
testb <- c(10,60,5,60,40,3)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,4)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age70wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,5] > 0 )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup370 <- lambda_cond_int(0.1475067 ,31.4912751  ,3.5169747  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup370, type = "l")
plot(event_count_upward, type ="l")

Englamdadown370 <- lambda_cond_int(22.6099003, 26.4025945 , 3.8386668, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown370, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England3") 

Engfulllambda370 <- as.data.frame(c(Englamdaup370, Englamdadown370))
Engfulllambda370$Deaths <- modelling_data.df[c(3:119),3]
Engfulllambda370$Date <- (modelling_data.df[c(3:119),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda370, aes(y= Engfulllambda370[,1], x = Engfulllambda370[,3])) +geom_line() +geom_point(aes(y=Engfulllambda370[,2],x=Engfulllambda370[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age70wave3death.csv", row.names = FALSE)







#####################################CASES##########################################################



#####Age 20 wave 2 cases #####


filter_data_cp_peak(Agecases20,"England2", 115)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(500,60,5,500,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age20wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220c <- lambda_cond_int(264.14834 , 28.04372  , 3.43524 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220c <- lambda_cond_int(111.53177 , 27.47241,   3.46082, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases20 %>%
  filter(Country == "England2") 

Engfulllambda220c <- as.data.frame(c(Englamdaup220c, Englamdadown220c))
Engfulllambda220c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda220c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda220c, aes(y= Engfulllambda220c[,1], x = Engfulllambda220c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda220c[,2],x=Engfulllambda220c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave2case.csv", row.names = FALSE)


#####Age 20 wave 3 cases #####


filter_data_cp_peak(Agecases20,"England3", 41)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(1000,60,5,1000,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age20wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup320c <- lambda_cond_int(419.908888 , 25.401143  , 3.343029 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup320c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown320c <- lambda_cond_int(508.720087 , 26.535996  , 3.411279, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown320c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases20 %>%
  filter(Country == "England3") 

Engfulllambda320c <- as.data.frame(c(Englamdaup320c, Englamdadown320c))
Engfulllambda320c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda320c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda320c, aes(y= Engfulllambda320c[,1], x = Engfulllambda320c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda320c[,2],x=Engfulllambda320c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave3case.csv", row.names = FALSE)



#####Age 45 wave 2 cases #####


filter_data_cp_peak(Agecases45,"England2", 115)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(500,60,5,500,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age45wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup245c <- lambda_cond_int(243.436050,  31.085087  , 3.554825 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup245c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown245c <- lambda_cond_int(  60.567796 , 30.624778 ,  3.579514, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown245c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases45 %>%
  filter(Country == "England2") 

Engfulllambda245c <- as.data.frame(c(Englamdaup245c, Englamdadown245c))
Engfulllambda245c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda245c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda245c, aes(y= Engfulllambda245c[,1], x = Engfulllambda245c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda245c[,2],x=Engfulllambda245c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age45wave2case.csv", row.names = FALSE)


#####Age 45 wave 3 cases #####


filter_data_cp_peak(Agecases45,"England3", 44)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(1000,60,5,1000,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age45wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup345c <- lambda_cond_int(230.792066 , 18.058510 ,  3.270232 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup345c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown345c <- lambda_cond_int(550.572499 , 34.547712 ,  3.708280, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown345c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases45 %>%
  filter(Country == "England3") 

Engfulllambda345c <- as.data.frame(c(Englamdaup345c, Englamdadown345c))
Engfulllambda345c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda345c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda345c, aes(y= Engfulllambda345c[,1], x = Engfulllambda345c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda345c[,2],x=Engfulllambda345c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age45wave3case.csv", row.names = FALSE)



#####Age 70 wave 2 cases #####


filter_data_cp_peak(Agecases70,"England2", 117)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(100,60,5,100,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age70wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup270c <- lambda_cond_int(36.821542 ,30.494311  ,3.494842  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup270c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown270c <- lambda_cond_int(3.172496 ,29.723579,  3.514228, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown270c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases70 %>%
  filter(Country == "England2") 

Engfulllambda270c <- as.data.frame(c(Englamdaup270c, Englamdadown270c))
Engfulllambda270c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda270c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda270c, aes(y= Engfulllambda270c[,1], x = Engfulllambda270c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda270c[,2],x=Engfulllambda270c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age70wave2case.csv", row.names = FALSE)


#####Age 70 wave 3 cases #####


filter_data_cp_peak(Agecases70,"England3", 90)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.1,0.1,0.1,0.1,0.1)
testb <- c(200,60,5,200,60,5)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,10000,10000,10000,10000)

test <- Neldermeadmulti(testa,testb,testc,testd,10000)
test <- read.csv("age70wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) & test[,2] >0 & test[,5] >0  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup370c <- lambda_cond_int(  15.182660 , 35.145134 ,  3.614573, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup370c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown370c <- lambda_cond_int(108.177817  ,32.460516 ,  3.648677, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown370c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases70 %>%
  filter(Country == "England3") 

Engfulllambda370c <- as.data.frame(c(Englamdaup370c, Englamdadown370c))
Engfulllambda370c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda370c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda370c, aes(y= Engfulllambda370c[,1], x = Engfulllambda370c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda370c[,2],x=Engfulllambda370c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age70wave3case.csv", row.names = FALSE)

















#############wave 1 deaths normal dist 10000 runs ###############################################
##########################################################################################################################
############################################################################
##############################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
######################################################################################################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
######################################################################################################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
############################################################################


filter_data_cp_peak(Englanddeaths,"England1", 31)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.5,0,0.1,0.5,0)
testb <- c(20,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,1,10000,10000,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)

test <- read.csv("sd1deathsgeometric.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(5,10), font = 1)
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(3,8))
grid()
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(0.5,1))
grid()

par(mfrow = c(1,1))
Englamdaup1 <- lambda_cond_int(8.9535496, 1.0978246 ,0.7777346  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup1, type = "l")
plot(event_count_upward, type ="l")

Englamdadown1 <- lambda_cond_int(6.4711249 ,0.9289044, 0.8118869, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown1, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England1") 

Engfulllambda1 <- as.data.frame(c(Englamdaup1, Englamdadown1))
Engfulllambda1$Deaths <- modelling_data.df[c(6:150),3]
Engfulllambda1$Date <- (modelling_data.df[c(6:150),1])

par(mfrow = c(1,1))
p3 <- ggplot(data = Engfulllambda1, aes(y= Engfulllambda1[,1], x = Engfulllambda1[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda1[,2],x=Engfulllambda1[,3], colour = "Daily counts"), size = 1.1, alpha = 0.8) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 1 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     date_labels = "%b")

p3

par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda1[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 1 deaths")
plot(ecdf(Engfulllambda1[,1]), col = "blue", add = TRUE)

#par(mfrow = c(1,1))
#Englamdaup1 <- lambda_cond_int(8.1479814 ,1.0914868, 0.8221077, decay_times_upward, decay_times_counts_upward)
#par(mfrow=c(2,1))
#plot(Englamdaup1, type = "l")
#plot(event_count_upward, type ="l")

#Englamdadown1 <- lambda_cond_int(05.3334575, 0.9331249, 0.8045498, decay_times_downward, decay_times_counts_downward)
#plot(Englamdadown1, type = "l")
#plot(event_count_downward, type ="l") 
write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd1deathsgeometric3.csv", row.names = FALSE)


#############wave 2 deaths normal dist 10000 runs ###############################################
############################################################################

filter_data_cp_peak(Englanddeaths,"England2", 127)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,0.5,0,0.1,0.5,0)
testb <- c(40,1.5,1,40,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,10000,1,10000,10000,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("sd2deathsgeometric.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)
par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(10,20), font = 1)
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(10,20))
grid()
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(0.5,1))
grid()

Englamdaup2 <- lambda_cond_int(15.7236383 , 0.9994692,  0.6709230 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2 <- lambda_cond_int(13.5947802,  0.9122952,  0.8209458, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englanddeaths %>%
  filter(Country == "England2") 

Engfulllambda2 <- as.data.frame(c(Englamdaup2, Englamdadown2))
Engfulllambda2$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda2$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda2, aes(y= Engfulllambda2[,1], x = Engfulllambda2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda2[,2],x=Engfulllambda2[,3]))


ggplot(data = Engfulllambda2, aes(y= Engfulllambda2[,1], x = Engfulllambda2[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda2[,2],x=Engfulllambda2[,3], colour = "Daily counts"), size = 1, alpha = 0.8) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     date_labels = "%b")

ggplot(data = Engfulllambda2, aes(y= Engfulllambda2[,1], x = Engfulllambda2[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda2[,2],x=Engfulllambda2[,3], colour = "Daily counts"), size = 1.1, alpha = 0.8) + xlab("Date (2020-2021)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 2 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",date_labels = "%b")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd2deathsgeometric.csv", row.names = FALSE)
ks.test(Engfulllambda2[,1],c(event_count_upward[,2],event_count_downward[,2]))

plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 deaths")
plot(ecdf(Engfulllambda2[,1]), col = "blue", add = TRUE)
###############wave 3 deaths normal dist #############################################
####################################################################################
filter_data_cp_peak(Englanddeaths, "England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0.1,1,0,0.1,0.5,0)
testb <- c(10,1.5,1,10,1.5,1)
testc <- c(0,1,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("sd3deathsgeometric2.csv")

goodtest <- subset(test,  test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(10,20), font = 1)
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(10,20))
grid()
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(0.5,1))
grid()

Englamdaup3 <- lambda_cond_int(2.1438526 ,1.0160771, 0.4440982, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3 <- lambda_cond_int(5.1649106, 0.9167003, 0.7132997, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3, type = "l")
plot(event_count_downward, type ="l") 

modelling_data.df = Englanddeaths %>%
  filter(Country == "England3") 

Engfulllambda3 <- as.data.frame(c(Englamdaup3, Englamdadown3))
Engfulllambda3$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3, aes(y= Engfulllambda3[,1], x = Engfulllambda3[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3[,2],x=Engfulllambda3[,3], colour = "Daily deaths")) + xlab("Date(2021)") + ylab("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths")

write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd3deathsgeometric2.csv", row.names = FALSE)

p5 <- ggplot(data = Engfulllambda3, aes(y= Engfulllambda3[,1], x = Engfulllambda3[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3[,2],x=Engfulllambda3[,3], colour = "Daily counts"), size = 1.1, alpha=0.8) + xlab("Date (2021)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     date_labels = "%b")
plot_grid(p3,p4,p5,ncol = 1)

ks.test(Engfulllambda3[,1],c(event_count_upward[,2],event_count_downward[,2]))

plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 deaths")
plot(ecdf(Engfulllambda3[,1]), col = "blue", add = TRUE)


ggplot(data = Engfulllambda3, aes(y= Engfulllambda3[,1], x = Engfulllambda3[,3], colour = "Fitted line")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3[,2],x=Engfulllambda3[,3], colour = "Daily counts"), size = 1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 3 Deaths") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",date_labels = "%b")

######################################Wave2 cases 10000 normal dist ###########################

filter_data_cp_peak(Englandcases,"England2", 142)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(200,0.5,0,200,0.5,0)
testb <- c(1000,1.5,1,1000,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,0.95,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("sd2casesgeometric2.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)))
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(250,350), font = 1)
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(450,550))
grid()
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(0.5,1))
grid()

Englamdaup2c <- lambda_cond_int(296.4100206 ,  1.0038764  , 0.9367666 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup2c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown2c <- lambda_cond_int(477.2801638 ,  0.9310588  , 0.9606038, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown2c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England2") 

Engfulllambda2c <- as.data.frame(c(Englamdaup2c, Englamdadown2c))
Engfulllambda2c$Deaths <- modelling_data.df[c(1:221),3]
Engfulllambda2c$Date <- (modelling_data.df[c(1:221),1])

par(mfrow = c(1,1))
p4 <- ggplot(data = Engfulllambda2c, aes(y= Engfulllambda2c[,1], x = Engfulllambda2c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda2c[,2],x=Engfulllambda2c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020-2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 2 Cases")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 date_labels = "%b")

p4

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd2casesgeometric2.csv", row.names = FALSE)

ks.test(Engfulllambda2c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 cases")
plot(ecdf(Engfulllambda2c[,1]), col = "blue", add = TRUE)

#########################WAVE 3 NORMAL DIST 10000 ITERATIONS###################################


filter_data_cp_peak(Englandcases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")





set.seed(1)
testa <- c(600,0.5,0,600,0.5,0)
testb <- c(1200,1.5,1,1200,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,3,1,10000,3,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("sd3casesgeometric.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(3,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.5, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[1])), col="red" , ylim = c(750,850), font = 1)
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg of ", mu[2])), ylim = c(850,950))
grid()
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg of ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg of ", beta[2])), ylim = c(0.5,1))
grid()
#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd3casesgeometric.csv", row.names = FALSE)



Englamdaup3c <- lambda_cond_int(819.0365540 ,  1.0814259  , 0.7382588 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup3c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3c <- lambda_cond_int(897.3051729 ,  0.9264869,   0.9920328, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown3c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Englandcases %>%
  filter(Country == "England3") 

Engfulllambda3c <- as.data.frame(c(Englamdaup3c, Englamdadown3c))
Engfulllambda3c$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3c$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3c, aes(y= Engfulllambda3c[,1], x = Engfulllambda3c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3c[,2],x=Engfulllambda3c[,3]))



ggplot(data = Engfulllambda3c, aes(y= Engfulllambda3c[,1], x = Engfulllambda3c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3c[,2],x=Engfulllambda3c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 3 Cases")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        date_labels = "%b")




ks.test(Engfulllambda3c[,1],c(event_count_upward[,2],event_count_downward[,2]))

plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 cases")
plot(ecdf(Engfulllambda3c[,1]), col = "blue", add = TRUE)





#################################### Vacs outside the hawkes sum deaths ###########

filter_data_cp_peak(Vacdatadeaths,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatadeaths,"England3",92,3)



set.seed(1)
testa <- c(0.1,-1,0.5,0,0.1,-1,0.5,0)
testb <- c(10,1,1.5,1,10,1,1.5,1)
testc <- c(0,-100,0,0,0,-100,0,0)
testd <- c(10000,10000,2,1,10000,10000,2,1)


test <- Neldermeadmultivacsoutside(testa,testb,testc,testd,5000)
test <- read.csv("sd3deathsvacoutsidegeometric.csv")

goodtest <- subset(test, test[,9] < quantile(test[,9], probs = c(.8)))
goodtest[,c(3,7)] <- log(goodtest[,c(3,7)])
paramavg <- c(0,0,0,0,0,0,0,0)
for(i in 1:8){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(3,7)] <- exp(paramavg[c(3,7)])
Vacloglikelihood(paramavg)


par(mfrow = c(4,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.3, 0.3), mgp=c(2,1,0), cex.lab = 1.02)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg ", mu[1])), col="red" , ylim = c(3,6), font = 1)
grid()
plot(cumsum(goodtest[,5]) / seq_along(goodtest[,5]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg ", mu[2])), ylim = c(3,6))
grid()
plot(cumsum(goodtest[,2]) / seq_along(goodtest[,2]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg ", gamma[1])), ylim = c(-0.1,0))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg ", gamma[2])), ylim = c(-0.1,0))
grid()
plot((cumsum((goodtest[,3])) / seq_along(goodtest[,3])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,7])) / seq_along(goodtest[,7])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="dark green", xlab = 'Iteration', ylab = expression(paste("Avg ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,8]) / seq_along(goodtest[,8]), type = 'l', col="orange", xlab = 'Iteration', ylab = expression(paste("Avg ", beta[2])), ylim = c(0.5,1))
grid()

Englamdaup3vac <- VAClambda_cond_int(4.90186621, -0.09014538  ,1.03610913,  0.64271289 , decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3vac, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac <- VAClambda_cond_int(4.97353910, -0.03865846,  0.93654874,  0.58438042, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3vac, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeaths %>%
  filter(Country == "England3") 

Engfulllambda3vac <- as.data.frame(c(Englamdaup3vac, Englamdadown3vac))
Engfulllambda3vac$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac$Date <- (modelling_data.df[c(3:121),1])




par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac, aes(y= Engfulllambda3vac[,1], x = Engfulllambda3vac[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac[,2],x=Engfulllambda3vac[,3]))


ggplot(data = Engfulllambda3vac, aes(y= Engfulllambda3vac[,1], x = Engfulllambda3vac[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3vac[,2],x=Engfulllambda3vac[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Deaths") + ggtitle("Vac Model 1 Output for Wave 3 Deaths")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",date_labels = "%b")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd3deathsvacoutsidegeometric.csv", row.names = FALSE)

Vaceffect <- subset(Vacdatadeaths[,c(1,4)], Vacdatadeaths[,2] == "England3")
Vaceffect[c(1:92),2] <- Vaceffect[c(1:92),2] * -0.09014538
Vaceffect[c(93:121),2] <- Vaceffect[c(93:121),2] * -0.03865846
Vaceffect$pctvac <- subset(Vacdatadeaths[,c(4)], Vacdatadeaths[,2] == "England3")
ggplot(Vaceffect, aes(x= pctvac, y = vac)) + geom_line() +ylim(-10,0)
ggplot(Vaceffect, aes(x= Date, y = vac)) + geom_line() +ylim(-10,0) + ylab("Isolated effect of vaccination") + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),plot.title = element_text(size = 14),legend.text = element_text(size = 12) )



filter_data_cp_peak(Englanddeaths,"England3", 92)
Confidencepct(Engfulllambda3, 100, 5)
Confidencepct(Engfulllambda3, 100, 10)
Confidencepct(Engfulllambda3, 100, 15)
Confidencepct(Engfulllambda3vac, 100, 5)
Confidencepct(Engfulllambda3vac, 100, 10)
Confidencepct(Engfulllambda3vac, 100, 15)

par(mfrow = c(1,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda3vac[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Vac model 1 deaths")
plot(ecdf(Engfulllambda3vac[,1]), col = "blue", add = TRUE)


#################################### Vacs outside the hawkes sum cases ###########
Vacdatacases <- merge(Englandcases,cumvacdata, by ="Date", all = T)
Vacdatacases[is.na(Vacdatacases)] <- 0


filter_data_cp_peak(Vacdatacases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
Extractvacs(Vacdatacases,"England3",52,1)



set.seed(1)
testa <- c(600,-1,0.5,0,600,-1,0.5,0)
testb <- c(1200,1,1.5,1,1200,1,1.5,1)
testc <- c(0,-1000,0,0,0,-1000,0,0)
testd <- c(10000,10000,2,1,10000,10000,2,1)

test <- Neldermeadmultivacsoutside(testa,testb,testc,testd,5000)
test <- read.csv("sd3casesvacoutsidegeometric.csv")

goodtest <- subset(test, test[,9] < quantile(test[,9], probs = c(.8)))
goodtest[,c(3,7)] <- log(goodtest[,c(3,7)])
paramavg <- c(0,0,0,0,0,0,0,0)
for(i in 1:8){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(3,7)] <- exp(paramavg[c(3,7)])
Vacloglikelihood(paramavg)

Englamdaup3vacc <- VAClambda_cond_int(452.8224226 , -8.3552146,   1.0753597,   0.8442123, decay_times_upward, decay_times_counts_upward, vacwaveup)
par(mfrow=c(2,1))
plot(Englamdaup3vacc, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vacc <- VAClambda_cond_int(896.0602933,  -1.2454781,   0.9284824,   0.9871591, decay_times_downward, decay_times_counts_downward, vacwavedown)
plot(Englamdadown3vacc, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3vacc <- as.data.frame(c(Englamdaup3vacc, Englamdadown3vacc))
Engfulllambda3vacc$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3vacc$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vacc, aes(y= Engfulllambda3vacc[,1], x = Engfulllambda3vacc[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vacc[,2],x=Engfulllambda3vacc[,3]))


ggplot(data = Engfulllambda3vacc, aes(y= Engfulllambda3vacc[,1], x = Engfulllambda3vacc[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3vacc[,2],x=Engfulllambda3vacc[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Vac Model 1 Output for Wave 3 Cases")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",date_labels = "%b")



par(mfrow = c(4,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.3, 0.3), mgp=c(2,1,0), cex.lab = 1.02)
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg ", mu[1])), col="red" , ylim = c(800,900), font = 1)
grid()
plot(cumsum(goodtest[,5]) / seq_along(goodtest[,5]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg ", mu[2])), ylim = c(850,950))
grid()
plot(cumsum(goodtest[,2]) / seq_along(goodtest[,2]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg ", gamma[1])), ylim = c(-10,0))
grid()
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg ", gamma[2])), ylim = c(-5,0))
grid()
plot((cumsum((goodtest[,3])) / seq_along(goodtest[,3])), type = 'l', col="blue", xlab = 'Iteration', ylab = expression(paste("Avg ", alpha[1])), ylim = c(0.5,1.5))
grid()
plot((cumsum((goodtest[,7])) / seq_along(goodtest[,7])), type = 'l', col="black", xlab = 'Iteration', ylab = expression(paste("Avg ", alpha[2])), ylim = c(0.5,1.5))
grid()
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="dark green", xlab = 'Iteration', ylab = expression(paste("Avg ", beta[1])), ylim = c(0.5,1))
grid()
plot(cumsum(goodtest[,8]) / seq_along(goodtest[,8]), type = 'l', col="orange", xlab = 'Iteration', ylab = expression(paste("Avg ", beta[2])), ylim = c(0.5,1))
grid()

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sd3casesvacoutsidegeometric.csv", row.names = FALSE)

Vaceffect <- subset(Vacdatadeaths[,c(1,4)], Vacdatadeaths[,2] == "England3")
Vaceffect[c(1:92),2] <- Vaceffect[c(1:92),2] * -8.3552146
Vaceffect[c(93:121),2] <- Vaceffect[c(93:121),2] * -1.2454781
Vaceffect$pctvac <- subset(Vacdatadeaths[,c(4)], Vacdatadeaths[,2] == "England3")
ggplot(Vaceffect, aes(x= pctvac, y = vac)) + geom_line() +ylim(-700,10)
ggplot(Vaceffect, aes(x= Date, y = vac)) + geom_line() +ylim(-700,10) + ylab("Vaccination effect") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) )





filter_data_cp_peak(Englandcases,"England2", 92)
Confidencepct(Engfulllambda2c, 100, 5)
Confidencepct(Engfulllambda2c, 100, 10)
Confidencepct(Engfulllambda2c, 100, 15)

filter_data_cp_peak(Englandcases,"England3", 52)
Confidencepct(Engfulllambda3c, 100, 5)
Confidencepct(Engfulllambda3c, 100, 10)
Confidencepct(Engfulllambda3c, 100, 15)
filter_data_cp_peak(Englandcases,"England3", 52)
Confidencepct(Engfulllambda3vac, 100, 5)
Confidencepct(Engfulllambda3vac, 100, 10)
Confidencepct(Engfulllambda3vac, 100, 15)


ks.test(Engfulllambda3vacc[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Vacs model 1 cases")
plot(ecdf(Engfulllambda3vacc[,1]), col = "blue", add = TRUE)



#################DAILY VAC deaths WITHIN THE HAWKES SUM ##################


filter_data_cp_peak(Vacdatadeathspct,"England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatadeathspct,"England3", 92,3)


set.seed(1)
testa <- c(0.1,1,0.5,-0.5,0,0.1,0.5,0.5,-0.5,0)
testb <- c(20,1.5,1,0.5,1,20,1,1,0.5,1)
testc <- c(0,0,0,-1000,0,0,0,0,-1000,0)
testd <- c(10000,2,1,10000,1,10000,2,1,10000,1)

test <- Neldermeadmultivacsoutside2(testa,testb,testc,testd,5000)
test <- read.csv("sddeathvacinside.csv")


goodtest <- subset(test, test[,11] < quantile(test[,11], probs = c(.8)))
goodtest[,c(2,4,7,9)] <- log(goodtest[,c(2,4,7,9)])
paramavg <- c(0,0,0,0,0,0,0,0,0,0)
for(i in 1:10){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,4,7,9)] <- exp(paramavg[c(2,4,7,9)])
paramavg[c(4,9)] <- -paramavg[c(4,9)]
paramavg
Vacloglikelihood2(paramavg)



par(mfrow = c(2,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 0.6, 0.3), mgp=c(2,1,0), cex.lab = 1.02)

plot(cumsum(goodtest[,4]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = expression(paste("Avg (deaths) ", gamma[1])), col="red" , ylim = c(-0.1,0.1), font = 1)
grid()
plot(cumsum(goodtest[,9]) / seq_along(goodtest[,5]), type = 'l', col="purple", xlab = 'Iteration', ylab = expression(paste("Avg (deaths) ", gamma[2])), ylim = c(-0.1,0.1))
grid()
plot(cumsum(goodtest[,4]+6) / seq_along(goodtest[,2]), type = 'l', col="green", xlab = 'Iteration', ylab = expression(paste("Avg (cases) ", gamma[1])), ylim = c(-1,1))
grid()
plot(cumsum(goodtest[,9]-2.9) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = expression(paste("Avg (cases) ", gamma[2])), ylim = c(-4,4))
grid()

Englamdaup3vac2 <- VAC2lambda_cond_int(8.05245770 , 0.88319566,  0.70683819,  0.04771116  ,0.53317275, decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3vac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3vac2 <- VAC2lambda_cond_int( 10.08406253 , 0.86877870,  0.77051962,  0.06323450, 0.54285583, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3vac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatadeathspct %>%
  filter(Country == "England3") 

Engfulllambda3vac2 <- as.data.frame(c(Englamdaup3vac2, Englamdadown3vac2))
Engfulllambda3vac2$Deaths <- modelling_data.df[c(3:121),3]
Engfulllambda3vac2$Date <- (modelling_data.df[c(3:121),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3]))



ggplot(data = Engfulllambda3vac2, aes(y= Engfulllambda3vac2[,1], x = Engfulllambda3vac2[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3vac2[,2],x=Engfulllambda3vac2[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Deaths") + ggtitle("Vac Model 2 Output for Wave 3 Deaths")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",date_labels = "%b")


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sddeathvacinside.csv", row.names = FALSE)

filter_data_cp_peak(Englanddeaths,"England3", 92)
Confidencepct(Engfulllambda3, 100, 5)
Confidencepct(Engfulllambda3, 100, 10)
Confidencepct(Engfulllambda3, 100, 15)
Confidencepct(Engfulllambda3vac2, 100, 5)
Confidencepct(Engfulllambda3vac2, 100, 10)
Confidencepct(Engfulllambda3vac2, 100, 15)


par(mfrow = c(1,2), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)

ks.test(Engfulllambda3vac2[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Vacs model 2 deaths")
plot(ecdf(Engfulllambda3vac2[,1]), col = "blue", add = TRUE)

#####Daily vac cases within the hawkes sum ##############
Vacdatacases <- merge(Englandcases,cumvacdata, by ="Date", all = T)
Vacdatacases[is.na(Vacdatacases)] <- 0



filter_data_cp_peak(Vacdatacases,"England3", 52)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")
filter_data_cp_peakvac(Vacdatacases,"England3", 52,1)


set.seed(1)
testa <- c(200,0,0,-0.5,0,200,0,0,-0.5,0)
testb <- c(500,2,1,0.5,1,500,2,1,0.5,1)
testc <- c(0,0,0,-1000,0,0,0,0,-1000,0)
testd <- c(10000,2,1,10000,1,10000,2,1,10000,1)

test <- Neldermeadmultivacsoutside2(testa,testb,testc,testd,5000)
test <- read.csv("sdcasesvacinside.csv")

goodtest <- subset(test, test[,11] < quantile(test[,11], probs = c(.8)) )
goodtest[,c(4,9)] <- -goodtest[,c(4,9)]
goodtest[,c(2,4,7,9)] <- log(goodtest[,c(2,4,7,9)])
paramavg <- c(0,0,0,0,0,0,0,0,0,0)
for(i in 1:10){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,4,7,9)] <- exp(paramavg[c(2,4,7,9)])
paramavg[c(4,9)] <- -paramavg[c(4,9)]
paramavg
Vacloglikelihood2(paramavg
                  )


Englamdaup3cvac2 <- VAC2lambda_cond_int(348.2559800 ,  1.0549663  , 0.9927746 , -6.0984405   ,0.5991013 , decay_times_upward, decay_times_counts_upward, decay_timesvac_upward, decay_times_countsvac_upward)
par(mfrow=c(2,1))
plot(Englamdaup3cvac2, type = "l")
plot(event_count_upward, type ="l")

Englamdadown3cvac2 <- VAC2lambda_cond_int(355.7877289  , 0.9427739,   0.9960142,   2.8911780,   0.4326736, decay_times_downward, decay_times_counts_downward, decay_timesvac_downward, decay_times_countsvac_downward)
plot(Englamdadown3cvac2, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Vacdatacases %>%
  filter(Country == "England3") 

Engfulllambda3cvac2 <- as.data.frame(c(Englamdaup3cvac2, Englamdadown3cvac2))
Engfulllambda3cvac2$Deaths <- modelling_data.df[c(1:68),3]
Engfulllambda3cvac2$Date <- (modelling_data.df[c(1:68),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda3cvac2, aes(y= Engfulllambda3cvac2[,1], x = Engfulllambda3cvac2[,3])) +geom_line() +geom_point(aes(y=Engfulllambda3cvac2[,2],x=Engfulllambda3cvac2[,3]))



ggplot(data = Engfulllambda3cvac2, aes(y= Engfulllambda3cvac2[,1], x = Engfulllambda3cvac2[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda3cvac2[,2],x=Engfulllambda3cvac2[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Vac Model 2 Output for Wave 3 Cases")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",date_labels = "%b")

#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdcasesvacinside.csv", row.names = FALSE)

filter_data_cp_peak(Englandcases,"England3", 52)
Confidencepct(Engfulllambda3c, 100, 5)
Confidencepct(Engfulllambda3c, 100, 10)
Confidencepct(Engfulllambda3c, 100, 15)
filter_data_cp_peak(Englandcases,"England3", 52)
Confidencepct(Englamdadown3cvac2, 100, 5)
Confidencepct(Englamdadown3cvac2, 100, 10)
Confidencepct(Englamdadown3cvac2, 100, 15)


ks.test(Engfulllambda3cvac2[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Vacs model 2 cases")
plot(ecdf(Engfulllambda3cvac2[,1]), col = "blue", add = TRUE)


######################################## BY age ###############################################################


Agedeaths <- read.csv("Deathsage.csv")
Agedeaths <- Agedeaths[,c(2,4,5,6)]
Agedeaths <- na.omit(Agedeaths)
Agedeaths <- Agedeaths %>% 
  rename(
    Date = date,
    Country = areaName,
    n = deaths,
    age = age
  )

Agedeaths <- Agedeaths[,c(2,1,4,3)]
Agedeaths <- Agedeaths %>% map_df(rev)
#Englanddeaths$n <- rollmean(Englanddeaths$n, 5, na.pad = TRUE , align = "center")
#Englanddeaths <- na.omit(Englanddeaths)
Agedeaths$n <- as.integer(Agedeaths$n)
Agedeaths$n <- as.numeric(Agedeaths$n)
Agedeaths$Date <- as.Date(Agedeaths$Date)
Agedeaths <- as.data.frame(Agedeaths)

Agedeaths <- subset(Agedeaths, !(Agedeaths$age %in% c("60+", "00_04", "00_59", "05_09", "10_14", "15_19","20_24")))

Country <- rep(0,605)
Country[1:150] <- rep("England1",150)
Country[195:400] <- rep("England2",206)
Country[460:580] <- rep("England3",121)
Country[151:194] <- rep("England",44)
Country[401:459] <- rep("England",59)
Country[580:605] <- rep("England",26)

Agedeaths20 <- subset(Agedeaths, (Agedeaths$age %in% c("25_29", "30_34", "35_39","40_44")))
Agedeaths20 <- aggregate(Agedeaths20[,c(3)], FUN="sum", by=list(Agedeaths20$Date))
Agedeaths20$Country <- Country
Agedeaths20 <- Agedeaths20 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths20 <- Agedeaths20[,c(1,3,2)]

Agedeaths45 <- subset(Agedeaths, (Agedeaths$age %in% c("45_49", "50_54", "55_59", "60_64", "65_69")))
Agedeaths45 <- aggregate(Agedeaths45[,c(3)], FUN="sum", by=list(Agedeaths45$Date))
Agedeaths45$Country <- Country
Agedeaths45 <- Agedeaths45 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths45 <- Agedeaths45[,c(1,3,2)]

Agedeaths70 <- subset(Agedeaths, (Agedeaths$age %in% c("70_74", "75_79", "80_84", "85_89", "89+")))
Agedeaths70 <- aggregate(Agedeaths70[,c(3)], FUN="sum", by=list(Agedeaths70$Date))
Agedeaths70$Country <- Country
Agedeaths70 <- Agedeaths70 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agedeaths70 <- Agedeaths70[,c(1,3,2)]

Agedeaths20$n <- rollmean(Agedeaths20$n, 5, na.pad = TRUE , align = "center")
Agedeaths20 <- na.omit(Agedeaths20)
Agedeaths20$n <- ceil(Agedeaths20$n)
Agedeaths45$n <- rollmean(Agedeaths45$n, 5, na.pad = TRUE , align = "center")
Agedeaths45 <- na.omit(Agedeaths45)
Agedeaths45$n <- ceil(Agedeaths45$n)
Agedeaths70$n <- rollmean(Agedeaths70$n, 5, na.pad = TRUE , align = "center")
Agedeaths70 <- na.omit(Agedeaths70)
Agedeaths70$n <- ceil(Agedeaths70$n)

###Visalisations###
Agedatadeaths <- merge(Agedeaths20, Agedeaths45 , by =c("Date", "Country"))
Agedatadeaths <- merge(Agedatadeaths, Agedeaths70 , by =c("Date", "Country"))
Agedatadeaths <- Agedatadeaths %>% 
  rename("20-44" = n.x,
         "45-69" = n.y,
         "70+" = n
  )
plot1 <- ggplot(Agedatadeaths, aes(x= Date, Agedatadeaths$`20-44`, colour = "20-44")) +geom_line() +geom_line(aes(x=Date,Agedatadeaths$`45-69`, colour = "45-69"))+geom_line(aes(x=Date,y=Agedatadeaths$`70+`, , colour = "70+")) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Daily death figures by Age group") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "2 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      date_labels = "%b")



ggplot(Agedatadeaths, aes(x= Date, Agedatadeaths$`20-44`, colour = "20-44")) +geom_line()  + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Daily death figures by in the 20-44 age group") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "2 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 date_labels = "%b")


Agecases <- read.csv("Casesage.csv")
Agecases <- Agecases[,c(2,4,5,6)]
Agecases <- na.omit(Agecases)
Agecases <- Agecases %>% 
  rename(
    Date = date,
    Country = areaName,
    n = cases,
    age = age
  )

Agecases <- Agecases[,c(2,1,4,3)]
Agecases <- Agecases %>% map_df(rev)
#Englanddeaths$n <- rollmean(Englanddeaths$n, 5, na.pad = TRUE , align = "center")
#Englanddeaths <- na.omit(Englanddeaths)
Agecases$n <- as.integer(Agecases$n)
Agecases$n <- as.numeric(Agecases$n)
Agecases$Date <- as.Date(Agecases$Date)
Agecases <- as.data.frame(Agecases)

Agecases <- subset(Agecases, !(Agecases$age %in% c("60+", "00_04", "00_59", "05_09", "10_14", "15_19","20_24")))
Agecases <- subset(Agecases, Agecases$Date > "2020-03-01")

Country <- rep(0,605)
Country[1:150] <- rep("England1",150)
Country[195:400] <- rep("England2",206)
Country[460:580] <- rep("England3",121)
Country[151:194] <- rep("England",44)
Country[401:459] <- rep("England",59)
Country[580:605] <- rep("England",26)

Agecases20 <- subset(Agecases, (Agecases$age %in% c("25_29", "30_34", "35_39","40_44")))
Agecases20 <- aggregate(Agecases20[,c(3)], FUN="sum", by=list(Agecases20$Date))
Agecases20$Country <- Country
Agecases20 <- Agecases20 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases20 <- Agecases20[,c(1,3,2)]

Agecases45 <- subset(Agecases, (Agecases$age %in% c("45_49", "50_54", "55_59", "60_64", "65_69")))
Agecases45 <- aggregate(Agecases45[,c(3)], FUN="sum", by=list(Agecases45$Date))
Agecases45$Country <- Country
Agecases45 <- Agecases45 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases45 <- Agecases45[,c(1,3,2)]

Agecases70 <- subset(Agecases, (Agecases$age %in% c("70_74", "75_79", "80_84", "85_89", "89+")))
Agecases70 <- aggregate(Agecases70[,c(3)], FUN="sum", by=list(Agecases70$Date))
Agecases70$Country <- Country
Agecases70 <- Agecases70 %>% 
  rename(
    Date = Group.1,
    Country = Country,
    n = x
  )
Agecases70 <- Agecases70[,c(1,3,2)]

Agecases20$n <- rollmean(Agecases20$n, 5, na.pad = TRUE , align = "center")
Agecases20 <- na.omit(Agecases20)
Agecases45$n <- rollmean(Agecases45$n, 5, na.pad = TRUE , align = "center")
Agecases45 <- na.omit(Agecases45)
Agecases70$n <- rollmean(Agecases70$n, 5, na.pad = TRUE , align = "center")
Agecases70 <- na.omit(Agecases70)

###Visalisations###
Agedatacases <- merge(Agecases20, Agecases45 , by =c("Date", "Country"))
Agedatacases <- merge(Agedatacases, Agecases70 , by =c("Date", "Country"))
Agedatacases <- Agedatacases %>% 
  rename("20-44" = n.x,
         "45-69" = n.y,
         "70+" = n
  )
plot2 <- ggplot(Agedatacases, aes(x= Date, Agedatacases$`20-44`, colour = "20-44")) +geom_line() +geom_line(aes(x=Date,Agedatacases$`45-69`, colour = "45-69"))+geom_line(aes(x=Date,y=Agedatacases$`70+`, , colour = "70+")) + xlab("Date (2020)") + ylab ("Cases") +ggtitle("Daily case figures by Age group") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "2 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     date_labels = "%b")


plot2

plot_grid(plot1,plot2, ncol = 1)

##########age 45 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths20,"England1", 29)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(2,1.5,1,2,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,1.5,1,10000,1.5,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
#test <- read.csv("wave2deathsnormal.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup120 <- lambda_cond_int(0.2509354,
                                 30.4917179,
                                 4.1188409
                                 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup120, type = "l")
plot(event_count_upward, type ="l")

Englamdadown120 <- lambda_cond_int(
  1.2369513,
  25.5067962,
  2.8406203, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown120, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England1") 

Engfulllambda120 <- as.data.frame(c(Englamdaup120, Englamdadown120))
Engfulllambda120$Deaths <- modelling_data.df[c(1:136),3]
Engfulllambda120$Date <- (modelling_data.df[c(1:136),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda120, aes(y= Engfulllambda120[,1], x = Engfulllambda120[,3])) +geom_line() +geom_point(aes(y=Engfulllambda120[,2],x=Engfulllambda120[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/age20wave1death.csv", row.names = FALSE)




#####Age 20 wave 2 deaths #####


filter_data_cp_peak(Agedeaths20,"England2", 120)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(2,1.5,1,2,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,1.5,1,10000,1.5,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age20wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220 <- lambda_cond_int(0.8003272 ,0.8120289 ,0.6236840, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220 <- lambda_cond_int(0.8686452, 0.7081667 ,0.7564169, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England2") 

Engfulllambda220 <- as.data.frame(c(Englamdaup220, Englamdadown220))
Engfulllambda220$Deaths <- modelling_data.df[c(11:206),3]
Engfulllambda220$Date <- (modelling_data.df[c(11:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda220, aes(y= Engfulllambda220[,1], x = Engfulllambda220[,3])) +geom_line() +geom_point(aes(y=Engfulllambda220[,2],x=Engfulllambda220[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage20wave2death.csv", row.names = FALSE)



#####Age 20 wave 3 deaths #####


filter_data_cp_peak(Agedeaths20, "England3", 92)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")


set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(2,1.5,1,2,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,1.5,1,10000,1.5,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age20wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220 <- lambda_cond_int(0.9339432, 0.6264573, 0.5460338, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220 <- lambda_cond_int(1.0122897, 0.8655758, 0.6539218, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths20 %>%
  filter(Country == "England2") 

Engfulllambda220 <- as.data.frame(c(Englamdaup220, Englamdadown220))
Engfulllambda220$Deaths <- modelling_data.df[c(10:120),3]
Engfulllambda220$Date <- (modelling_data.df[c(10:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda220, aes(y= Engfulllambda220[,1], x = Engfulllambda220[,3])) +geom_line() +geom_point(aes(y=Engfulllambda220[,2],x=Engfulllambda220[,3]))


write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage20wave2death.csv", row.names = FALSE)


##########age 45 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths45,"England1", 29)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(10,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age45wave1death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup145 <- lambda_cond_int(4.7897359, 1.1416158 ,0.6286505, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup145, type = "l")
plot(event_count_upward, type ="l")

Englamdadown145 <- lambda_cond_int(6.6631338, 0.7626478, 0.7556282, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown145, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England1") 

Engfulllambda145 <- as.data.frame(c(Englamdaup145, Englamdadown145))
Engfulllambda145$Deaths <- modelling_data.df[c(7:147),3]
Engfulllambda145$Date <- (modelling_data.df[c(7:147),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda145, aes(y= Engfulllambda145[,1], x = Engfulllambda145[,3])) +geom_line() +geom_point(aes(y=Engfulllambda145[,2],x=Engfulllambda145[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage45wave1death.csv", row.names = FALSE)


ggplot(data = Engfulllambda145, aes(y= Engfulllambda145[,1], x = Engfulllambda145[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda145[,2],x=Engfulllambda145[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020)") + ylab ("Deaths") + ggtitle("Hawkes Model Output for Wave 1 Deaths for ages 45-69")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   date_labels = "%b")


par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda145[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 1 deaths 45-69")
plot(ecdf(Engfulllambda145[,1]), col = "blue", add = TRUE)

#####Age 45 wave 2 deaths #####


filter_data_cp_peak(Agedeaths45,"England2", 126)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(10,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age45wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup245 <- lambda_cond_int(4.1315191 ,0.9746539 ,0.6352272, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup245, type = "l")
plot(event_count_upward, type ="l")

Englamdadown245 <- lambda_cond_int( 7.1516521 ,0.8387810 ,0.7444404, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown245, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England2") 

Engfulllambda245 <- as.data.frame(c(Englamdaup245, Englamdadown245))
Engfulllambda245$Deaths <- modelling_data.df[c(5:206),3]
Engfulllambda245$Date <- (modelling_data.df[c(5:206),1])

par(mfrow = c(1,1))
plot3 <- ggplot(data = Engfulllambda245, aes(y= Engfulllambda245[,1], x = Engfulllambda245[,3], colour = "Fitted line")) +geom_line(size = 1.4) +geom_point(aes(y=Engfulllambda245[,2],x=Engfulllambda245[,3], colour = "Daily counts"), size = 1.4) + xlab("Date (2020)") + ylab ("Deaths") +ggtitle("Hawkes Model Output for Wave 2 deaths in age 45-69") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           date_labels = "%b")

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage45wave2death.csv", row.names = FALSE)

par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda245[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 deaths 45-69")
plot(ecdf(Engfulllambda245[,1]), col = "blue", add = TRUE)

#####Age 45 wave 3 deaths #####


filter_data_cp_peak(Agedeaths45, "England3", 104)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(10,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age45wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)))
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup345 <- lambda_cond_int( 2.7890434,  0.8409931,  0.5624359, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup345, type = "l")
plot(event_count_upward, type ="l")

Englamdadown345 <- lambda_cond_int(10.1286066,  0.5853363,  0.6597069, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown345, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths45 %>%
  filter(Country == "England3") 

Engfulllambda345 <- as.data.frame(c(Englamdaup345, Englamdadown345))
Engfulllambda345$Deaths <- modelling_data.df[c(4:119),3]
Engfulllambda345$Date <- (modelling_data.df[c(4:119),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda345, aes(y= Engfulllambda345[,1], x = Engfulllambda345[,3])) +geom_line() +geom_point(aes(y=Engfulllambda345[,2],x=Engfulllambda345[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage45wave3death.csv", row.names = FALSE)

ggplot(data = Engfulllambda345, aes(y= Engfulllambda345[,1], x = Engfulllambda345[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda345[,2],x=Engfulllambda345[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Deaths") + ggtitle("Hawkes Model Output for Wave 3 Deaths for ages 45-69")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         date_labels = "%b")
par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda345[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 deaths 45-69")
plot(ecdf(Engfulllambda345[,1]), col = "blue", add = TRUE)

##########age 70 + wave 1 death #############################################


filter_data_cp_peak(Agedeaths70,"England1", 32)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(10,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)
test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age70wave1death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)) )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup170 <- lambda_cond_int( 4.3082906, 1.0841971, 0.8236651 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup170, type = "l")
plot(event_count_upward, type ="l")

Englamdadown170 <- lambda_cond_int(5.4367533 ,0.9126579 ,0.7714841, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown170, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England1") 

Engfulllambda170 <- as.data.frame(c(Englamdaup170, Englamdadown170))
Engfulllambda170$Deaths <- modelling_data.df[c(4:146),3]
Engfulllambda170$Date <- (modelling_data.df[c(4:146),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda170, aes(y= Engfulllambda170[,1], x = Engfulllambda170[,3])) +geom_line() +geom_point(aes(y=Engfulllambda170[,2],x=Engfulllambda170[,3]))


ggplot(data = Engfulllambda170, aes(y= Engfulllambda170[,1], x = Engfulllambda170[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda170[,2],x=Engfulllambda170[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020)") + ylab ("Deaths") + ggtitle("Hawkes Model Output for Wave 1 Deaths for ages 70+")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         date_labels = "%b")


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage70wave1death.csv", row.names = FALSE)


par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda170[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 1 deaths 70+")
plot(ecdf(Engfulllambda170[,1]), col = "blue", add = TRUE)

##########age 70 + wave 2 death #############################################


filter_data_cp_peak(Agedeaths70,"England2", 129)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(10,1.5,1,20,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)
test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age70wave2death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)


par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup270 <- lambda_cond_int(4.0984866, 1.0145848 ,0.7609752 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup270, type = "l")
plot(event_count_upward, type ="l")

Englamdadown270 <- lambda_cond_int(5.9688581, 0.9187230, 0.8193549, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown270, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England2") 

Engfulllambda270 <- as.data.frame(c(Englamdaup270, Englamdadown270))
Engfulllambda270$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda270$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda270, aes(y= Engfulllambda270[,1], x = Engfulllambda270[,3])) +geom_line() +geom_point(aes(y=Engfulllambda270[,2],x=Engfulllambda270[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage70wave2death.csv", row.names = FALSE)


ggplot(data = Engfulllambda270, aes(y= Engfulllambda270[,1], x = Engfulllambda270[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda270[,2],x=Engfulllambda270[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020)") + ylab ("Deaths") + ggtitle("Hawkes Model Output for Wave 2 Deaths for ages 70+")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           date_labels = "%b")
par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda270[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 deaths 70+")
plot(ecdf(Engfulllambda270[,1]), col = "blue", add = TRUE)

###################################AGE 70 WAVE 3 DEATHS ####################

filter_data_cp_peak(Agedeaths70, "England3", 93)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,10,0.5,0)
testb <- c(4,1.5,1,50,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)
test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age70wave3death.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8)))
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(10,50))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(10,20))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(20,50))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup370 <- lambda_cond_int( 0.7919019 , 1.0150411,  0.5635343, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup370, type = "l")
plot(event_count_upward, type ="l")

Englamdadown370 <- lambda_cond_int(29.2438452,  0.4933094,  0.7071519, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown370, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agedeaths70 %>%
  filter(Country == "England3") 

Engfulllambda370 <- as.data.frame(c(Englamdaup370, Englamdadown370))
Engfulllambda370$Deaths <- modelling_data.df[c(3:119),3]
Engfulllambda370$Date <- (modelling_data.df[c(3:119),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda370, aes(y= Engfulllambda370[,1], x = Engfulllambda370[,3])) +geom_line() +geom_point(aes(y=Engfulllambda370[,2],x=Engfulllambda370[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage70wave3death.csv", row.names = FALSE)




ggplot(data = Engfulllambda370, aes(y= Engfulllambda370[,1], x = Engfulllambda370[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda370[,2],x=Engfulllambda370[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Deaths") + ggtitle("Hawkes Model Output for Wave 3 Deaths for ages 70+")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       date_labels = "%b")


par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda370[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 deaths 70+")
plot(ecdf(Engfulllambda370[,1]), col = "blue", add = TRUE)

#####################################CASES##########################################################



#####Age 20 wave 2 cases #####


filter_data_cp_peak(Agecases20,"England2", 115)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(500,1.5,1,500,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age20wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup220c <- lambda_cond_int(201.4202636 ,  0.9903044,   0.9907526  , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup220c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown220c <- lambda_cond_int(132.6631525,   0.9326351,   0.9705047, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown220c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases20 %>%
  filter(Country == "England2") 

Engfulllambda220c <- as.data.frame(c(Englamdaup220c, Englamdadown220c))
Engfulllambda220c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda220c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
plot4 <- ggplot(data = Engfulllambda220c, aes(y= Engfulllambda220c[,1], x = Engfulllambda220c[,3], colour = "Fitted line")) +geom_line(size = 1.4) +geom_point(aes(y=Engfulllambda220c[,2],x=Engfulllambda220c[,3], colour = "Daily counts"), size = 1.4) + xlab("Date (2020)") + ylab ("Cases") +ggtitle("Hawkes Model Output for Wave 2 cases in ages 20-44") + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ) + scale_colour_discrete(name = "Key")+ scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               date_labels = "%b")

plot_grid(plot3, plot4 , nrow =2)
#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage20wave2case.csv", row.names = FALSE)


par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda220c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 cases 20-44")
plot(ecdf(Engfulllambda220c[,1]), col = "blue", add = TRUE)

#####Age 20 wave 3 cases #####


filter_data_cp_peak(Agecases20,"England3", 41)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")

set.seed(1)
testa <- c(0,0.5,0,4000,0.5,0)
testb <- c(1000,1.5,1,6000,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age20wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup320c <- lambda_cond_int( 485.9455380 ,   1.1073095,    0.7789149, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup320c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown320c <- lambda_cond_int(4243.0172252 ,   0.4772270,    0.9441937, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown320c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases20 %>%
  filter(Country == "England3") 

Engfulllambda320c <- as.data.frame(c(Englamdaup320c, Englamdadown320c))
Engfulllambda320c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda320c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda320c, aes(y= Engfulllambda320c[,1], x = Engfulllambda320c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda320c[,2],x=Engfulllambda320c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage20wave3case.csv", row.names = FALSE)


ggplot(data = Engfulllambda320c, aes(y= Engfulllambda320c[,1], x = Engfulllambda320c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda320c[,2],x=Engfulllambda320c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 3 Cases for ages 20-44")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       date_labels = "%b")

par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda320c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 cases 20-44")
plot(ecdf(Engfulllambda320c[,1]), col = "blue", add = TRUE)

#####Age 45 wave 2 cases #####


filter_data_cp_peak(Agecases45,"England2", 115)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,1,0.5,0,0.5,0)
testb <- c(1000,1.5,1,1000,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age45wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup245c <- lambda_cond_int(349.9262865 , 0.9591448,   0.9567228, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup245c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown245c <- lambda_cond_int(  214.3431693  , 0.9077502  , 0.9222305, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown245c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases45 %>%
  filter(Country == "England2") 

Engfulllambda245c <- as.data.frame(c(Englamdaup245c, Englamdadown245c))
Engfulllambda245c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda245c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda245c, aes(y= Engfulllambda245c[,1], x = Engfulllambda245c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda245c[,2],x=Engfulllambda245c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage45wave2case.csv", row.names = FALSE)

ggplot(data = Engfulllambda245c, aes(y= Engfulllambda245c[,1], x = Engfulllambda245c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda245c[,2],x=Engfulllambda245c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020-2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 2 Cases for ages 45-69")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            date_labels = "%b")
par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda245c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 cases 45-64")
plot(ecdf(Engfulllambda245c[,1]), col = "blue", add = TRUE)

#####Age 45 wave 3 cases #####


filter_data_cp_peak(Agecases45,"England3", 44)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(200,1.5,1,1000,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,1.5,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age45wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup345c <- lambda_cond_int(94.4840170,   1.0200919   ,0.9906567, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup345c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown345c <- lambda_cond_int(436.1110142  , 0.9167318,   0.9983176, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown345c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases45 %>%
  filter(Country == "England3") 

Engfulllambda345c <- as.data.frame(c(Englamdaup345c, Englamdadown345c))
Engfulllambda345c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda345c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda345c, aes(y= Engfulllambda345c[,1], x = Engfulllambda345c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda345c[,2],x=Engfulllambda345c[,3]))


#write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage45wave3case.csv", row.names = FALSE)


ggplot(data = Engfulllambda345c, aes(y= Engfulllambda345c[,1], x = Engfulllambda345c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda345c[,2],x=Engfulllambda345c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 3 Cases for ages 45-69")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 date_labels = "%b")
par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda345c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 cases 45-64")
plot(ecdf(Engfulllambda345c[,1]), col = "blue", add = TRUE)

#####Age 70 wave 2 cases #####


filter_data_cp_peak(Agecases70,"England2", 117)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(200,1.5,1,100,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age70wave2case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot((cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(0,2))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(0,1))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot((cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(0,2))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(0,1))


Englamdaup270c <- lambda_cond_int(65.7206163 , 0.9743309,  0.9378179, decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup270c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown270c <- lambda_cond_int(28.9209288,  0.9150555,  0.8480038, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown270c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases70 %>%
  filter(Country == "England2") 

Engfulllambda270c <- as.data.frame(c(Englamdaup270c, Englamdadown270c))
Engfulllambda270c$Deaths <- modelling_data.df[c(1:206),3]
Engfulllambda270c$Date <- (modelling_data.df[c(1:206),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda270c, aes(y= Engfulllambda270c[,1], x = Engfulllambda270c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda270c[,2],x=Engfulllambda270c[,3]))


write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage70wave2case.csv", row.names = FALSE)


ggplot(data = Engfulllambda270c, aes(y= Engfulllambda270c[,1], x = Engfulllambda270c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda270c[,2],x=Engfulllambda270c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2020-2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 2 Cases for ages 70+")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            date_labels = "%b")

par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda270c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 2 cases 70+")
plot(ecdf(Engfulllambda270c[,1]), col = "blue", add = TRUE)

#####Age 70 wave 3 cases #####


filter_data_cp_peak(Agecases70,"England3", 90)
plot(event_count_upward[,2], type = "l")
plot(event_count_downward[,2], type = "l")



set.seed(1)
testa <- c(0,0.5,0,0,0.5,0)
testb <- c(200,1.5,1,200,1.5,1)
testc <- c(0,0,0,0,0,0)
testd <- c(10000,2,1,10000,2,1)

test <- Neldermeadmulti(testa,testb,testc,testd,5000)
test <- read.csv("age70wave3case.csv")

goodtest <- subset(test, test[,7] < quantile(test[,7], probs = c(.8))  )
goodtest[,c(2,5)] <- log(goodtest[,c(2,5)])
paramavg <- c(0,0,0,0,0,0)
for(i in 1:6){
  paramavg[i] <- mean(goodtest[,i])
}
paramavg[c(2,5)] <- exp(paramavg[c(2,5)])

loglikelihood(paramavg)

par(mfrow = c(2,3))
plot(cumsum(goodtest[,1]) / seq_along(goodtest[,1]), type = 'l', xlab = 'Iteration', ylab = 'Running Average of a1', col="red" , ylim = c(0,5))
plot(exp(cumsum((goodtest[,2])) / seq_along(goodtest[,2])), type = 'l', col="blue", xlab = 'Iteration', ylab = 'Running Average of b1', ylim = c(50,10))
plot(cumsum(goodtest[,3]) / seq_along(goodtest[,3]), type = 'l', col="green", xlab = 'Iteration', ylab = 'Running Average of c1', ylim = c(3,5))
plot(cumsum(goodtest[,4]) / seq_along(goodtest[,4]), type = 'l', col="purple", xlab = 'Iteration', ylab = 'Running Average of a2', ylim = c(0,5))
plot(exp(cumsum((goodtest[,5])) / seq_along(goodtest[,5])), type = 'l', col="black", xlab = 'Iteration', ylab = 'Running Average of b2', ylim = c(30,70))
plot(cumsum(goodtest[,6]) / seq_along(goodtest[,6]), type = 'l', col="brown", xlab = 'Iteration', ylab = 'Running Average of c2', ylim = c(3,5))


Englamdaup370c <- lambda_cond_int(76.7856562,  0.9738012,  0.6160530 , decay_times_upward, decay_times_counts_upward)
par(mfrow=c(2,1))
plot(Englamdaup370c, type = "l")
plot(event_count_upward, type ="l")

Englamdadown370c <- lambda_cond_int(99.1654799 , 0.9017397 , 0.7524824, decay_times_downward, decay_times_counts_downward)
plot(Englamdadown370c, type = "l")
plot(event_count_downward, type ="l")  

modelling_data.df = Agecases70 %>%
  filter(Country == "England3") 

Engfulllambda370c <- as.data.frame(c(Englamdaup370c, Englamdadown370c))
Engfulllambda370c$Deaths <- modelling_data.df[c(1:120),3]
Engfulllambda370c$Date <- (modelling_data.df[c(1:120),1])

par(mfrow = c(1,1))
ggplot(data = Engfulllambda370c, aes(y= Engfulllambda370c[,1], x = Engfulllambda370c[,3])) +geom_line() +geom_point(aes(y=Engfulllambda370c[,2],x=Engfulllambda370c[,3]))


write.csv(test,"C:/Users/afbda/Desktop/Hawkes processes clean/sdage70wave3case.csv", row.names = FALSE)


ggplot(data = Engfulllambda370c, aes(y= Engfulllambda370c[,1], x = Engfulllambda370c[,3], colour = "Fitted values")) +geom_line(size = 1.2) +geom_point(aes(y=Engfulllambda370c[,2],x=Engfulllambda370c[,3], colour = "Daily count"), size =1.1, alpha = 0.8) + xlab("Date (2021)") + ylab ("Cases") + ggtitle("Hawkes Model Output for Wave 3 Cases for ages 70+")+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 14),legend.text = element_text(size = 12) ,legend.title=element_text(size=13)) + scale_colour_discrete(name = "Key") + scale_x_date(date_breaks = "1 month",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               date_labels = "%b")
par(mfrow = c(2,3), mai = c(0.8,0.5,0.45,0.3), cex = 1, mar = c(4, 4, 1.1, 0.5), mgp=c(2,1,0), cex.lab = 1.1)
ks.test(Engfulllambda370c[,1],c(event_count_upward[,2],event_count_downward[,2]))
plot(ecdf(c(event_count_upward[,2], event_count_downward[,2])), col = "red", main = "Wave 3 cases 70+")
plot(ecdf(Engfulllambda370c[,1]), col = "blue", add = TRUE)





