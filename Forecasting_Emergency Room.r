# Case Study for Admissions in an Emergency Department

# Ken Zabiy Muhammad ARIEF
# Fahreza Baskara HEDIANDRA

# Input Library
library(zoo)
library(xts)
library(tidyverse)
library(tseries)
library(forecast)

setwd("PATH_TO_THE_CSV_FILE")
# Read the input file using read.zoo
wr_zoo <- read.zoo("data_emergency.csv", header=TRUE, sep = ",", FUN=as.POSIXct, format="%Y-%m-%d %H:%M:%S", tz='UTC')
# Transform it into an xts object
wr_xts <- xts(wr_zoo)

# Create a regular time index to comply with the time series framework
step = 60 #in minutes
seqmins <- seq(ISOdatetime(2020,12,7,0,0,0), ISOdatetime(2020,12,13,23,59,59), as.difftime(step, units="mins"))
# Create a new time series using this new regular index 
wr_regxts <- xts(na.omit(na.locf(merge(xts(,seqmins), wr_xts))[seqmins]), order.by=seqmins, frequency=24*60/step)
head(wr_regxts, n=30)

# Plotting to see the data distribution
wr_baseplot <- ggplot(data=wr_regxts, aes(x=Index, y=x)) +
  geom_line(aes(x=Index, y=x, color="original data")) +
  scale_color_manual(name="Patients", values=c("original data"="darkblue")) +
  xlab("Month") + ylab("Nb of Patients") +
  ggtitle("Evolution of the Number of Patients (7-13 Dec 2020)")
print(wr_baseplot)

ggplot(data = NULL, aes(x = index(wr_regxts) , y = wr_regxts$x)) +
  geom_line(color="darkblue")

# Data Analysis 

#Split test-train
train.period <- '/2020-12-12'
trainset <- wr_regxts[train.period] 
trainset
test.period <- '2020-12-13/2020-12-14'
testset <- wr_regxts[test.period]
testset

#Additive data breakdown
trainset <- ts(trainset, frequency=24)
plot(decompose(as.ts(trainset)))

#Multiplicative data breakdown
plot(decompose(as.ts(trainset), type='multiplicative'))

#Test of Stationary
library(urca)
summary(ur.df(log(trainset)))
trainset.d1 <- diff.xts(log(trainset))
summary(ur.df(na.omit(trainset.d1)))
trainset.d2 <- diff.xts(log(trainset), differences = 2)
summary(ur.df(na.omit(trainset.d2)))

#Validate the Seasonality
library(TSA)
pg <-periodogram(na.omit(trainset.d1))

#Check Period of Season
pdat <- data.frame(freq=pg$freq, amp=pg$spec)
pord <- pdat[order(pdat$amp, decreasing=TRUE),]
freq5 <- head(pord, 5)
print(1/freq5$freq)

#The frequency is then added to the beginning of the coding


#Forecasting - Linear Regression
#Forecasting, creating the equation
wr_regxts$log.dat <- log(wr_regxts$x) #creating new column based on the logarithmic value of x
train.id <- seq(from=1, to=nrow(wr_regxts[train.period]), by=1)
test.id <- seq(nrow(wr_regxts[train.period])+1, nrow(wr_regxts[train.period])+(1*24))
allperiods.id <- c(train.id, test.id) #
length(allperiods.id)
nrow(wr_regxts)

#Define functions
wr_regxts$f1 <- allperiods.id
wr_regxts$f2 <- c(sin(2*pi*allperiods.id/24))
wr_regxts$f3 <- c(cos(2*pi*allperiods.id/24))
print(tail(wr_regxts, 1*24))

#Use the log data
linreg.model <- lm(log.dat ~ f1 + f2 + f3,
                   data = wr_regxts,
                   subset = train.id)
summary(linreg.model)

#Use the original data
gen.linreg.model <- glm(x ~ f1 + f2 + f3,
                        family = gaussian(link = "log"),
                        data = wr_regxts,
                        subset = index(wr_regxts) < '2020-12-13')
summary(gen.linreg.model)

#Create new explanatory variable
new.explvar <- wr_regxts[test.period, grepl('f1|f2|f3', names(wr_regxts))]
forecast.lm <- predict.lm(linreg.model, new.explvar, h=(24),
                          interval = "prediction", se.fit = TRUE, level = 0.95)
forecast.glm <- predict.glm(gen.linreg.model, new.explvar, h=24,
                            interval = "prediction", se.fit = TRUE)
wr_regxts$lm.fit <- exp(c(linreg.model$fitted.values, forecast.lm$fit[,'fit']))
wr_regxts$lm.lwr <- c(rep(NA, nrow(wr_regxts[train.period])),
                     exp(forecast.lm$fit[,'lwr']))
wr_regxts$lm.upr <- c(rep(NA, nrow(wr_regxts[train.period])),
                     exp(forecast.lm$fit[,'upr']))

#Link Function
ilink <- family(gen.linreg.model)$linkinv
wr_regxts$glm.fit <- c(gen.linreg.model$fitted.values, ilink(forecast.glm$fit))
wr_regxts$glm.lwr <- c(rep(NA, nrow(wr_regxts[train.period])),
                      ilink(qnorm(0.025,
                                  mean=forecast.glm$fit,
                                  sd = forecast.glm$se.fit)))
wr_regxts$glm.upr <- c(rep(NA, nrow(wr_regxts[train.period])),
                      ilink(qnorm(0.975,
                                  mean=forecast.glm$fit,
                                  sd = forecast.glm$se.fit)))
#Show Data
tail(wr_regxts, 16)

#Plot the Linear Model
plot.lm <- wr_baseplot +
  geom_ribbon(data = wr_regxts, #Interval
              aes(ymin = lm.lwr, ymax = lm.upr, fill="LM", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = lm.fit, color="LM")) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "LM"="red")) +
  scale_fill_manual(name = "Prediction interval",
                    values=c("LM"="red"))
print(plot.lm)

#Plot the General Linear Model
plot.glm <- plot.lm +
  geom_ribbon(data = wr_regxts, #Interval
              aes(ymin = glm.lwr, ymax = glm.upr, fill="GLM", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = glm.fit, color="GLM")) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "LM"="red",
                              "GLM"="green"))+
  scale_fill_manual(name ="Prediction intervals",
                    values=c("LM"="red",
                             "GLM"="green"))
print(plot.glm)


# Forecast - Autoregressive
ar.model <- ar(wr_regxts$log.dat[train.period], order.max=24) #Function AR
forecast.ar <- forecast(ar.model, h=(24), level = 0.95)
wr_regxts$ar.fit <- c(rep(NA, nrow(wr_regxts[train.period])),
                     exp(forecast.ar$mean))
wr_regxts$ar.lwr <- c(rep(NA, nrow(wr_regxts[train.period])),
                     exp(forecast.ar$lower))
wr_regxts$ar.upr <- c(rep(NA, nrow(wr_regxts[train.period])),
                     exp(forecast.ar$upper))

#Plot AR
plot.ar <- wr_baseplot +
  geom_ribbon(data = wr_regxts,
              aes(ymin = ar.lwr, ymax = ar.upr, fill="AR", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = ar.fit, color="AR")) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "AR"="orange")) +
  scale_fill_manual(name = "Prediction interval",
                    values=c("AR"="orange"))
print(plot.ar)

#Accuracy LM-GLM-AR
forecast.res <- as.data.frame(wr_regxts[test.period])
accuracy(forecast.res$lm.fit, forecast.res$x)
accuracy(forecast.res$glm.fit, forecast.res$x)
accuracy(forecast.res$ar.fit, forecast.res$x)


#Forecast - Exponential Smoothing

#Single ES
es.model <- ets(wr_regxts$x[train.period], model="MNN")
forecast.es <- forecast.ets(es.model, h=(24), level = 0.95)
wr_regxts$es.fit <- c(es.model$fitted, forecast.es$mean)
wr_regxts$es.lwr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.es$lower)
wr_regxts$es.upr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.es$upper)

#Double ES
hw.model <- ets(wr_regxts$x[train.period], model="MMN")
forecast.hw <- forecast.ets(hw.model, h=(24), level = 0.95)
wr_regxts$hw.fit <- c(hw.model$fitted, forecast.hw$mean)
wr_regxts$hw.lwr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.hw$lower)
wr_regxts$hw.upr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.hw$upper)

#Triple ES
hws.model <- ets(wr_regxts$x[train.period], model="MMM")
forecast.hws <- forecast.ets(hws.model, h=(24), level = 0.95)
wr_regxts$hws.fit <- c(hws.model$fitted, forecast.hws$mean)
wr_regxts$hws.lwr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.hws$lower)
wr_regxts$hws.upr <- c(rep(NA, nrow(wr_regxts[train.period])), forecast.hws$upper)

#Plotting ES
plot.esmodels <- wr_baseplot +
  geom_ribbon(data = wr_regxts,
              aes(ymin = es.lwr, ymax = es.upr, fill="ES", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = es.fit, color="ES")) +
  geom_ribbon(data = wr_regxts,
              aes(ymin = hw.lwr, ymax = hw.upr, fill="HW", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = hw.fit, color="HW")) +
  geom_ribbon(data = wr_regxts,
              aes(ymin = hws.lwr, ymax = hws.upr, fill="HWS", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = hws.fit, color="HWS")) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "ES"="yellow",
                              "HW"="green",
                              "HWS"="purple")) +
  scale_fill_manual(name = "Prediction interval",
                    values=c("ES"="yellow",
                             "HW"="green",
                             "HWS"="purple"))
print(plot.esmodels)

#Acuracy for Exponential Smothing
forecast.res <- as.data.frame(wr_regxts[test.period])
accuracy(forecast.res$es.fit, forecast.res$x)
accuracy(forecast.res$hw.fit, forecast.res$x)
accuracy(forecast.res$hws.fit, forecast.res$x)


#Forecast - ARIMA

#To get q = 1,2,4
acf(na.omit(diff.xts(wr_regxts$x, log=TRUE)))

#To get p = 1,2,6
pacf(na.omit(diff.xts(wr_regxts$x, log=TRUE)))

#To get Q
acf(na.omit(diff.xts(wr_regxts$x, lag=24, log=TRUE))) #lag=seasonality?

#To get P
pacf(na.omit(diff.xts(wr_regxts$x, lag=24, log=TRUE)))

#Accuracy Calculation for ARIMA
model_smry.arima <- data.frame()

#Accuracy Loop for ARIMA
p_list <- list(1,2,6)
q_list <- list(1,2,4)

for(p in p_list){
  for(q in q_list){
    for(d in 0:0){
      arima.model <- Arima(wr_regxts$log.dat[train.period], order=c(p,d,q), 
                          seasonal=list(order=c(2,0,0), period=24), method="CSS")
      forecast <- forecast(arima.model, h=24, level=0.95)
      wr_regxts$fit <- exp(c(arima.model$fitted, forecast$mean))
      forecast.res <- as.data.frame(wr_regxts[test.period])
      acc <- accuracy(forecast.res$fit, forecast.res$x)
      # gather everything into a single data frame 
      acc_ext <- data.frame(# information from accuracy function
        acc,
        # arima order
        p,
        d,
        q)
      # add arima summary
      model_smry.arima <- rbind(model_smry.arima, acc_ext)
    }
  }
}  

#Show summary
model_smry.arima

#Calculate Choosen ARIMA
arima202.model <- Arima(wr_regxts$log.dat[train.period], order=c(2,0,2),
                        seasonal=list(order=c(2,0,0), period=(24), method="CSS"))
forecast.a202 <- forecast(arima202.model, h=(24), level = 0.95)
wr_regxts$a202.fit <- exp(c(arima202.model$fitted, forecast.a202$mean))
wr_regxts$a202.lwr <- c(rep(NA, nrow(wr_regxts[train.period])),
                        exp(forecast.a202$lower))
wr_regxts$a202.upr <- c(rep(NA, nrow(wr_regxts[train.period])),
                        exp(forecast.a202$upper))
wr_regxts[test.period]

#Plotting Choosen ARIMA
plot.arimamodels <- wr_baseplot +
  geom_ribbon(data = wr_regxts,
              aes(ymin = a202.lwr, ymax = a202.upr,
                  fill="ARIMA(2,0,2)", color=NULL),
              alpha=0.1) +
  geom_line(data = wr_regxts, aes(y = a202.fit, color="ARIMA(2,0,2)")) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "ARIMA(2,0,2)"="red")) +
  scale_fill_manual(name = "Prediction interval",
                    values=c("ARIMA(2,0,2)"="red"))
print(plot.arimamodels)

#Testing the Model into another period of data set
arima202.model <- Arima(wr_regxts$log.dat[train.period], order=c(2,0,2),
                        seasonal=list(order=c(2,0,0), period=(24), method="CSS"))

# Create a regular time index to comply with the time series framework
step = 60 #in minutes
seqtrials <- seq(ISOdatetime(2021,1,28,0,0,0), ISOdatetime(2021,1,29,23,59,59), as.difftime(step, units="mins"))
# Create a new time series using this new regular index 
trialset <- xts(na.omit(na.locf(merge(xts(,seqtrials), wr_xts))[seqtrials]), order.by=seqtrials, frequency=24*60/step)
head(trialset, n=30)

trialtrain.period <- '/2021-01-28'
trialtrainset <- trialset[trialtrain.period] 
trialtrainset
trialtest.period <- '2021-01-29/2021-01-30'
trialtestset <- trialset[trialtest.period]
trialtestset

arima202.model <- Arima(wr_regxts$log.dat[train.period], order=c(2,0,2),
                        seasonal=list(order=c(2,0,0), period=(24), method="CSS"))
forecast.trial <- forecast(arima202.model, h=(24), level = 0.95)
trialset$trial.lwr <- c(rep(NA, nrow(trialset[trialtrain.period])),
                        exp(forecast.trial$lower))
trialset$trial.upr <- c(rep(NA, nrow(trialset[trialtrain.period])),
                        exp(forecast.trial$upper))

# Plotting the new trial data distribution
trial_baseplot <- ggplot(data=trialset, aes(x=Index, y=x)) +
  geom_line(aes(x=Index, y=x, color="original data")) +
  scale_color_manual(name="Patients", values=c("original data"="darkblue")) +
  xlab("Month") + ylab("Nb of Patients") +
  ggtitle("Evolution of the Number of Patients (28-30 January 2021)")
print(trial_baseplot)

#Plotting The Trial
plot.trialmodels <- trial_baseplot +
  geom_ribbon(data = trialset,
              aes(ymin = trial.lwr, ymax = trial.upr,
                  fill="ARIMA(2,0,2)", color=NULL),
              alpha=0.1) +
  scale_color_manual(name="Fitted predictions",
                     values=c("original data"="darkblue",
                              "ARIMA(2,0,2)"="red")) +
  scale_fill_manual(name = "Prediction interval",
                    values=c("ARIMA(2,0,2)"="red"))
print(plot.trialmodels)

# Merci Beaucoup !