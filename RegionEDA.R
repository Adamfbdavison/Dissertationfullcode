setwd("C:/Users/afbda/Desktop/New folder/RegionData")

library(ggplot2)
library(scales)
library(zoo)
library(grid)
library(viridis)
library(gridExtra)
library(cowplot)
library(egg)
library(lattice)
library(hawkes)
#remove.packages(hawkesbow)
#devtools::install_github("fcheysson/hawkesbow", force = TRUE)
library(hawkesbow)
install.packages("ptproc")

#Reading data#
GenderAgeCases <- read.csv("RegionGenderCasesAge.csv")
DailyCasesDeaths <- read.csv("RegionNewCasesDeaths.csv")
AgeVacCasesDeaths <- read.csv("RegionVacCaseDeath.csv")
GenderAgeCases$date <- as.Date(GenderAgeCases$date, format = "%d/%m/%Y")
DailyCasesDeaths$date <- as.Date(DailyCasesDeaths$date, format = "%d/%m/%Y")
AgeVacCasesDeaths$date <- as.Date(AgeVacCasesDeaths$date, format = "%d/%m/%Y")


Nationcases <- read.csv("NationCase.csv")
NationCasesage <- read.csv("NationCasesbyAge.csv")
NationDeathsage <- read.csv("NationDeathsbyAge.csv")
NationVacage <- read.csv("NationVacbyAge.csv")
Nationcases$date <- as.Date(Nationcases$date)
NationCasesage$date <- as.Date(NationCasesage$date)
NationDeathsage$date <- as.Date(NationDeathsage$date)
NationVacage$date <- as.Date(NationVacage$date)

weeklyR <- read.csv("NationRvalue.csv")
weeklyR$date <- as.Date(weeklyR$date, format = "%d/%m/%Y")
dailyR <- data.frame(date = seq(as.Date("2020-05-29"), as.Date("2021-11-12"),"days"))
R = merge(dailyR,weeklyR, by = "date", all =TRUE)
R = na.locf(R)

Pop <- read.csv("Englandpop.csv")
Pop$Population <- as.numeric(gsub(",","",Pop$Population,fixed=TRUE))

#####NATION EDA###########

#Combination of the age and vacination data

Nationdata <- merge(NationVacage,NationCasesage, by = c("date", "age", "areaCode", "areaName", "areaType"))
Nationdata <- merge(Nationdata,NationDeathsage, by = c("date", "age", "areaCode", "areaName", "areaType"))

#For NA values i chose to repeat the previous days figure as due to the nature of infections I believe this should be a sufficient estimate
summary(Nationdata)
Nationdata$age <- factor(Nationdata$age)
Nationdata <- merge(Nationdata, Pop, by = "age")
Nationdata <- merge(Nationdata, R[,c(1,5,6,7)], by = "date")
Nationdata$Population <- sapply(Nationdata$Population, as.numeric)


summary(Nationcases)
Nationcases <- na.locf(Nationcases)

summary(DailyCasesDeaths)
DailyCasesDeaths <- na.locf(DailyCasesDeaths)

#Looking at the national distribution of cases and deaths

rollingnationcases <- rollmean(Nationcases$newCasesByPublishDate,7, fill = 0, align = "left")
Nationcases$rollingmean <- rollingnationcases

colfunc <- colorRampPalette(c("white", "red"))

grad_by_val <- function(x, y, cols = blues9) {
  
  y <- y[order(x)]
  ys <- (y - min(y)) / diff(range(y))
  cols <- colorRamp(cols)(ys) / 256
  colnames(cols) <- c("red", "green", "blue")
  cols <- apply(cols, 1, function(z) do.call(rgb, as.list(z)))
  im <- matrix(cols, ncol = length(x))
  rasterGrob(
    image = im,
    width = unit(0.91, "npc"),
    height = unit(1, "npc"),
    y = 0.55,
    interpolate = TRUE
  )
}


ggplot(Nationcases, aes(x = date, y = Nationcases$newCasesByPublishDate)) + annotation_custom(
  grob = grad_by_val(Nationcases$date, Nationcases$newDeaths28DaysByDeathDate, cols = colfunc(80))) +
  geom_point(shape = 4, alpha = 0.5) +geom_area(aes( y=Nationcases$rollingmean), fill ="light blue", alpha = 0.75, outline.type = "upper") +geom_line(aes(y=Nationcases$rollingmean)) + xlab("Date (Jan 2020 to Nov 2021)") + ylab("New cases") + ggtitle("Daily Covid-19 cases in England") + scale_x_date(breaks = date_breaks("months"),labels = date_format("%b"))

ggplot(Nationcases, aes(x = date, y=Nationcases$rollingmean)) + annotation_custom(
  grob = grad_by_val(Nationcases$date, Nationcases$newDeaths28DaysByDeathDate, cols = colfunc(80))) +
  geom_area(fill ="light blue", alpha = 0.75, outline.type = "full")+geom_line(aes(y=Nationcases$rollingmean)) +xlab("Date (Jan 2020 to Nov 2021)") + ylab("New cases") + ggtitle("Daily Covid-19 cases in England") + scale_x_date(breaks = date_breaks("months"),labels = date_format("%b"))

ggplot(Nationcases, aes(x = date, y=Nationcases$newDeaths28DaysByDeathDate)) + annotation_custom(
  grob = grad_by_val(Nationcases$date, Nationcases$newCasesByPublishDate, cols = colfunc(80))) +
  geom_area(fill ="light blue", alpha = 0.75, outline.type = "full")+geom_line(aes(y=Nationcases$newDeaths28DaysByDeathDate)) +xlab("Date (Jan 2020 to Nov 2021)") + ylab("Number of deaths") + ggtitle("Daily Covid-19 deaths in England") + scale_x_date(breaks = date_breaks("months"),labels = date_format("%b")) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )



p3 <- ggplot(Nationdata, aes(x=age, y = cases)) +geom_boxplot(color="darkblue", fill="lightblue") + ylab("Number of cases") + xlab("Age range") + ggtitle(("Distribution of daily cases by age group"))  + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )

p4 <- ggplot(Nationdata, aes(x=age, y = deaths)) +geom_boxplot(color="darkblue", fill="lightblue")+ ylab("Number of deaths") + xlab("Age range") + ggtitle(("Distribution of daily deaths by age group")) + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )
ggarrange(p3, p4, heights = c(0.5, 0.5))

ggplot(Nationdata, aes(x=factor(age), y = VaccineRegisterPopulationByVaccinationDate, group = 1)) +geom_point(color="darkblue", fill="lightblue") + geom_line() +xlab("Age range") + ylab("Population on vaccine register") +ggtitle(("Population on the vaccine register by age group"))# + geom_point(aes(y=Population))

ggplot(Nationdata, aes(x=date, y = cumVaccinationSecondDoseUptakeByVaccinationDatePercentage, color = age)) + geom_line() + xlab("Date") + ylab("Percentage fully vaccinated") + scale_x_date(breaks = date_breaks("months"),labels = date_format("%b"))
                                                                                                                                                                                              

caserate <- ggplot(Nationdata, aes(x=date, y = rollingRate.x, color = age)) + geom_line( show.legend = TRUE)+ theme(legend.position = c(1.07,0.75)) + xlab("Date (Dec 2020 - Nov 2021)") + ylab("Rolling case rate")+ scale_x_date(breaks = date_breaks("months"),labels = date_format("%b"))+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) ) + theme(plot.margin=unit(c(0.3,2.6,0.3,0.3),"cm"))
deathrate <- ggplot(Nationdata, aes(x=date, y = rollingRate.y, color = age)) + geom_line() +xlab("Date") + ylab("Rolling death rate")+ scale_x_date(breaks = date_breaks("months"),labels = date_format("%b"))+ theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )
ggarrange(caserate, deathrate, widths  = c(0.5, 0.5))

p1 <- ggplot(Nationdata, aes(x = date, y = newPeopleVaccinatedFirstDoseByVaccinationDate ,color= age)) +geom_point(show.legend = FALSE) +  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) + ylab("No. of first doses per day")+ scale_x_date(breaks = date_breaks("months"),labels = date_format("%b")) + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) ) + scale_y_continuous(labels = comma)
p2 <- ggplot(Nationdata, aes(x = date, y = newPeopleVaccinatedSecondDoseByVaccinationDate ,color= age)) +geom_point() + theme(legend.position = c(0.95,1.2)) + ylab("No. of second doses per day")+ xlab("Date")+ scale_x_date(breaks = date_breaks("months"),labels = date_format("%b")) + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) )
ggarrange(p1, p2, heights = c(0.5, 0.5))

Rvalue <- ggplot(Nationdata, aes(x=date, y = AvgR)) + geom_line() + xlab("Date (Dec 2020 - Nov 2021)") + ylab("R value")+  theme(axis.title.x = element_blank(),axis.text.x = element_blank()) + scale_x_date(breaks = date_breaks("months"),labels = date_format("%b")) + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),plot.title = element_text(size = 16),legend.text = element_text(size = 12) ) + theme(plot.margin=unit(c(0.3,2.6,0.3,0.3),"cm"))
ggarrange(Rvalue, caserate, heights  = c(0.5, 0.5))

ggplot(Nationdata, aes(x=cumVaccinationCompleteCoverageByVaccinationDatePercentage, y = cumVaccinationSecondDoseUptakeByVaccinationDatePercentage, color = age)) + geom_line() + xlab("Pct 'complete' vacccinated") + ylab("Percentage fully vaccinated")



###Basic models###

par(mfrow = c(2,1))

lattice::xyplot(Nationdata$rollingRate.y ~ Nationdata$cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, groups=Nationdata$age,
       auto.key=list(title="Age", corner=c(0.1, 1)))
lattice::xyplot(Nationdata$rollingRate.y ~ Nationdata$cumVaccinationSecondDoseUptakeByVaccinationDatePercentage, groups=Nationdata$age,
                auto.key=list(title="Age", corner=c(1, 1)))




lm <- lm(rollingRate.y ~ 0 + age + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, Nationdata)
summary(lm)
plot(lm)
ggplot(aes(x = fitted(lm), y = rstandard(lm)), data = Nationdata )+ xlab("Fitted values")+ ylab(" standardised residuals") + geom_point(aes(colour = age), shape = 18)
lm2 <- lm(rollingRate.y ~ 0 + age + cumVaccinationSecondDoseUptakeByVaccinationDatePercentage, Nationdata)
summary(lm2)
plot(lm2)
ggplot(aes(x = fitted(lm2), y = rstandard(lm2)), data = Nationdata )+ xlab("Fitted values")+ ylab(" standardised residuals") + geom_point(aes(colour = age), shape = 18)

lm3 <- lm(rollingRate.y ~  0 + cases + age + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, Nationdata)
summary(lm3)
par(mfrow = c(2,2))
plot(lm3)
lm4 <- lm(rollingRate.y ~  0 + cases + age + cumVaccinationSecondDoseUptakeByVaccinationDatePercentage, Nationdata)
summary(lm4)
par(mfrow = c(2,2))
plot(lm4)

agelm1 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "25_29"))
agelm2 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "30_34"))
agelm3 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "35_39"))
agelm4 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "40_44"))
agelm5 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "45_49"))
agelm6 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "50_54"))
agelm7 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "55_59"))
agelm8 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "60_64"))
agelm9 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "65_69"))
agelm10 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "70_74"))
agelm11 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "75_79"))
agelm12 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "80_84"))
agelm13 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "85_89"))
agelm14 <- lm(rollingRate.y ~ 0  + cases + cumVaccinationFirstDoseUptakeByVaccinationDatePercentage, subset(Nationdata, age == "90+"))
summary(agelm1)
summary(agelm14)
plot(agelm14)

#Amalgamate by region so we can work with the national data when required
#CaseDeaths <- aggregate(DailyCasesDeaths[ ,5:7], FUN="sum", by=list(DailyCasesDeaths$date))
#CaseDeaths[,1] <- as.Date(CaseDeaths[,1])

#ggplot(DailyCasesDeaths, aes(x = date, y = DailyCasesDeaths$newCasesByPublishDate, colour = DailyCasesDeaths$areaName))+geom_line()  + facet_wrap(~DailyCasesDeaths$areaName) 
#ggplot(DailyCasesDeaths, aes(x = date, y = DailyCasesDeaths$newDeaths, colour = DailyCasesDeaths$areaName))+geom_line()  + facet_wrap(~DailyCasesDeaths$areaName)


set.seed(102030)




x <- hawkes(100, fun = 0.7, repr = 0.7, family = "exp", rate = 2)
par(mfrow = c(2,1))
plot(x, intensity = TRUE, precision = 1e3)
plot(x, intensity = FALSE, precision = 1e3)
xdiscrete <- discrete(x, binsize = 0.1)
xparams <- mle(x$p, kern = "exp", end = 100)
xparams$par
xdisparams <- whittle(xdiscrete, kern = "exp")
xdisparams$par
plot(intensity(x,0:100))
plot(xdiscrete, type = "l")
par(mfrow = c(1,1))
y <- hawkes(20, fun = 0.7, repr = 0.7, family = "exp", rate = 0.800)
plot(y, intensity = TRUE, precision = 1e3, cex.axis = 2)
plot(y, intensity = FALSE, precision = 1e3, cex.axis = 2)

yparams <- mle(y$p, kern = "exp", end = 20)
yparams$par

hawkesdata <- Nationcases[,c(4,5)]
time <- seq(1,639,1)
hawkesdata$time <- time
hawkesdata <- hawkesdata[,c(2,3)]


arrivals <- function(x){
  arrivaldata <- c()
  for(i in 1:length(x)){
    times <- seq(i-1,i,1/x[i])
    arrivaldata <- c(arrivaldata,times)
  }
  return((arrivaldata))
}
##cases 2nd peak
par(mfrow = c(1,1))
tester <- Nationcases[c(295:355),5]
tester <- rev(tester)
tester <- as.numeric(tester)
plot(tester)

arrivaldata <- arrivals(tester)

params <- mle(arrivaldata, kern = "exp", end = 60)
params$par
params$model$ddloglik(arrivaldata, 60)
MLE <- hawkes(60, fun = 11.85811        , repr = 0.9999000     , family = "exp", rate = 97595.83065)
discretize <- discrete(MLE, binsize=1)
par(mfrow = c(2,1))
plot(MLE, intensity = TRUE, precision = 1e3)
plot(tester, type ="l")
plot(discretize, type ="l")


Discrete <- whittle(tester, kern= "exp", binsize = 1)
Discrete$par
DiscreteMLE <- hawkes(60, fun = 6477.1197         , repr = 0.9900    , family = "exp", rate = 10.4093)
discretize1 <- discrete(DiscreteMLE, binsize=1)
par(mfrow = c(1,1))
plot(DiscreteMLE, intensity = TRUE, precision = 1e2)
abline(tester1, type ="l")

##deaths 2nd peak

testerd <- Nationcases[c(565:600),6]
testerd <- rev(testerd)
testerd <- as.numeric(testerd)
plot(testerd)
arrivaldata <- arrivals(testerd)

params <- mle(arrivaldata, kern = "exp", end = 35)
params$par
params$model$ddloglik(arrivaldata, 35)
MLEd <- hawkes(35, fun = 6.436058     , repr = 0.999900     , family = "exp", rate = 1915.543637)
par(mfrow = c(1,1))
plot(MLEd, intensity = TRUE, precision = 1e3)
plot(testerd, type ="l")


Discrete <- whittle(testerd, kern= "exp", binsize = 1)
Discrete$par
DiscreteMLE <- hawkes(35, fun = 107.6940616       , repr = 0.9684457      , family = "exp", rate = 2.4665819)
par(mfrow = c(2,1))
deathsdiscrete <- discrete(DiscreteMLE, binsize =1)
plot(DiscreteMLE, intensity = TRUE, precision = 1e2)

plot(testerd, type ="l")
plot(deathsdiscrete, type = "l")


averagehawkes <- function(data, n){
  datamat <- matrix(0,n, length(data)-1)
  for(i in 1:n) {
    DiscreteMLE <- hawkes(35, fun = 107.6940616       , repr = 0.9684457      , family = "exp", rate = 2.4665819)
    deathsdiscrete <- discrete(DiscreteMLE, binsize =1)
    datamat[i,] <- deathsdiscrete
  }
  return(datamat)
}

sample1000 <- averagehawkes(testerd,1000)

mataverage <- function(data){
  x <- matrix(0,1,ncol(data))
  for(i in 1:ncol(data)){
    x[,i] <- mean(data[,i])
  }
  return(x)
}
avgsample <- mataverage(sample1000)

testerd2 <- Nationcases[,6]
testerd2 <- rev(testerd2)
testerd2 <- as.numeric(testerd2)

Discrete2 <- whittle(testerd2, kern= "exp", binsize = 1)
Discrete2$par
Discrete2MLE <- hawkes(639, fun = 1.198099        , repr = 0.990000       , family = "exp", rate = 2.207828)
par(mfrow = c(1,1))
plot(Discrete2MLE, intensity = TRUE, precision = 1e2)
abline(testerd2, type ="l")










library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#devtools::install_github("ImperialCollegeLondon/epidemia")

library(epidemia)
library(lubridate)
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(EpiEstim)
data("Flu1918")
print(Flu1918)
date <- as.Date("1918-01-01") + seq(0, along.with = c(NA, Flu1918$incidence))
data <- data.frame(city = "Baltimore", cases = c(NA, Flu1918$incidence),
                   date = date)
rt <- epirt(formula = R(city, date) ~ 1 + rw(prior_scale = 0.01),
            prior_intercept = normal(log(2), 0.2), link = 'log')
obs <-  epiobs(formula = cases ~ 0 + offset(rep(1,93)), link = "identity",
               i2o = rep(.25,4))
inf <- epiinf(gen = Flu1918$si_distr)
inf_ext <-  epiinf(gen = Flu1918$si_distr, latent = TRUE,
                   prior_aux = normal(10,2))
args <- list(rt = rt, obs = obs, inf = inf, data = data, iter = 2e3,
             seed = 12345)
args_ext <- args; args_ext$inf <- inf_ext
system.time(fm1 <- do.call(epim, args))
system.time(fm2 <- do.call(epim, args_ext))
spaghetti_rt(fm1)
spaghetti_rt(fm2)


dataC <- subset(Nationdata, Nationdata$age == "55_59")
dataC <- dataC[,c(1,16,19)]

date <- dataC[,1]
data <- data.frame(region = "England", cases = dataC$cases,
                date = date)
rt <- epirt(formula = R(region, date) ~ 1 + rw(prior_scale = 0.01),
            prior_intercept = normal(log(2), 0.2), link = 'log')
obs <-  epiobs(formula = cases ~ 0 + offset(rep(1,323)), link = "identity",
               i2o = rep(1/7,7))
#inf <- epiinf(gen = Flu1918$si_distr)
inf_ext <-  epiinf(gen = rep(1/323,323), latent = TRUE,
                   prior_aux = normal(10,2))
args <- list(rt = rt, obs = obs, inf = inf_ext, data = data, iter = 100,
             seed = 12345)
#args_ext <- args; args_ext$inf <- inf_ext
system.time(fm1 <- do.call(epim, args))


dataA <- subset(Nationdata, Nationdata$age == "55_59")
dataA <- dataA[,c(1,16,19)]
date <- dataA[,1]
#dataB <- filter(dataA, date > date[which(cumsum(deaths) > 10)[3] - 30])
rt <- epirt(formula = R(region, date) ~ 0 ,
            prior = shifted_gamma(shape = 1/6, scale = 1, shift = log(1.05)/6),
            prior_covariance = decov(shape = c(2, rep(0.5, 5)), scale = 0.25),
            link = scaled_logit(6.5))
inf <- epiinf(gen = rep(1/323,323), seed_days = 6)
deaths <- epiobs(formula = cases ~ 1, i2o = rep(6,323),
                 prior_intercept = normal(0,0.2), link = scaled_logit(0.02))
args <- list(rt = rt, inf = inf, obs = deaths, data = data, seed = 12345,
             refresh = 0)
pr_args <- c(args, list(algorithm = "sampling", iter = 1000, prior_PD = TRUE))
fm_prior <- do.call(epim, pr_args)
plot_rt(fm_prior)
