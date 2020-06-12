############################################################################################
######                              Term Paper                                        ######
######  Forecast and Forecast Evaluation of the Australian inflation rate             ######
######  - Does adding more variables improve the forecast's accuracy?                 ######
######                                                                                ######
###### Name: Niklas Landsberg                                                         ######
###### T-Number: 940421-T857                                                          ######
######                                                                                ######
############################################################################################

############################################################################################
########### 1.0 Libraries & Functions ######################################################
############################################################################################
library(readr)
library(dplyr)
library(forecast)
library(urca)
library(lmtest)
library(vars)
library(quantmod)
library(portes)
library(dynlm)
library(tseries)
library(readxl)

ADF.test <- function(data) {
  ADF <- function(type, data) {
    require(urca)
    result1 <- ur.df(data, type = type, lags = 3 * frequency(data), 
                     selectlags = "AIC")
    DETERM <- ifelse(type == "trend", 2, ifelse(type == "drift", 
                                                1, 0))
    LAGS <- length(coefficients(result1@testreg)[, 1]) - 
      DETERM - 1
    result2 <- cbind(t(result1@teststat), result1@cval, coefficients(result1@testreg)["z.lag.1", 
                                                                                      1], LAGS)
    round(result2, 2)
  }
  types <- c("trend", "drift", "none")
  result3 <- apply(t(types), 2, ADF, data)
  cat(rep("#", 20), "\n")
  cat(rep("#", 20), "\n")
  cat("Augmented Dickey--Fuller test\n")
  cat(rep("#", 20), "\n")
  cat("type:", "  trend ", "drift ", "none\n")
  cat("AR1:   ", result3[[1]][1, 5] + 1, " ", result3[[2]][1, 
                                                           5] + 1, " ", result3[[3]][5] + 1, "\n")
  cat("lags:  ", result3[[1]][1, 6], "   ", result3[[2]][1, 
                                                         6], "   ", result3[[3]][6], "\n")
  cat(rep("#", 20), "\n")
  result5 <- rbind(result3[[1]][c(1, 3), 1:4], result3[[2]][1:2, 
                                                            1:4], result3[[3]][1:4])
  rownames(result5)[5] <- "tau1"
  result5
}

bias.test <- function(h,FE){
  require(dynlm); require(lmtest);require(sandwich)
  model <- dynlm(FE[,h] ~ 1)
  matrix <- NeweyWest(model,
                      lag = h-1)
  round(coeftest(model, 
                 vcov. = matrix)[[4]],
        2)
}

DM.TEST <- function(ld, h) {
  require(dynlm); require(lmtest); require(sandwich) 
  fe1 <- get(ld[1]); fe2 <- get(ld[2])
  res <- dynlm(I(fe1[, h]^2 - fe2[, h]^2) ~ 1)
  mat <- NeweyWest(res, lag = h - 1) 
  round(coeftest(res, vcov. = mat)[4], 2)
}

JJ.test <- function (data, test = c("trace", "eigen"), extdata = NULL) {
  if (is.ts(data) == FALSE) 
    stop("Data is not a ts-object!")
  else VS <- function(type, data) {
    VARselect(y = data, type = type)
  }
  types <- t(c("none", "const", "trend"))
  result1 <- apply(types, 2, VS, data)
  result2 <- sapply(result1, function(x) x$selection[[4]])
  result3 <- list(list(types[1], result2[1]), list(types[2], 
                                                   result2[2]), list(types[3], result2[3]))
  TT <- function(Type, data) {
    if (Type[[2]] < 2) 
      Type[[2]] <- 2
    TSTAT <- rev(ca.jo(data, type = test, ecdet = Type[[1]], 
                       K = Type[[2]], dumvar = extdata)@teststat)
    CVAL <- ca.jo(data, type = test, ecdet = Type[[1]], K = Type[[2]], 
                  dumvar = extdata)@cval
    CVAL <- CVAL[nrow(CVAL):1, 2]
    list(TSTAT, CVAL)
  }
  result5 <- lapply(result3, TT, data)
  result6 <- cbind(result5[[1]][[1]], result5[[1]][[2]], result5[[2]][[1]], 
                   result5[[2]][[2]], result5[[3]][[1]], result5[[3]][[2]])
  result7 <- round(result6, 2)
  if (test == "trace") 
    TESTTYPE <- "Johansen's tracetest."
  else TESTTYPE <- "Maximal eigenvalue test."
  colnames(result7) <- c("tstat none", "cval none", "tstat const", 
                         "cval const", "tstat trend", "cval trend")
  cat(rep("#", 40), "\n")
  cat(rep("#", 40), "\n")
  cat("The JOHANSEN PROCEDURE\n")
  cat(rep("#", 40), "\n")
  cat("Test: ", paste(TESTTYPE), "\n")
  cat("Lags selected by AIC. Level of significance: 5 %.\n")
  cat("lags: none = ", result2[1], " const = ", result2[2], 
      " trend = ", result2[3], "\n")
  cat(rep("#", 40), "\n")
  result7
}




############################################################################################
########### 2 Introduction #################################################################
############################################################################################
########### 2.1 Dataset  ###################################################################
############################################################################################ 

raw_pi <- read_excel('RBA_pi.xls') # Loading the inflation rate data

pi <- raw_pi[12:402, # Quarterly inflation rate
             12]
pi <- ts(pi,
         start = c(1982,3),
         end = c(2020,1),
         frequency = 4)

raw_pi_exp <- read_excel('RBA_pi_exp.xls') # Loading the inflation rate expectation data

pi_exp <- raw_pi_exp[26:148, # 3-month Business inflation expectation
             3]

pi_exp <- ts(pi_exp,
             start = c(1989,3),
             end = c(2020,3),
             frequency = 4)

raw_rgdp <- read_excel('RBA_GDP.xls') # Loading the GDP data

rgdp <- raw_rgdp[15:252, #Year-ended real GDP growth
                 3]

rgdp <- ts(rgdp,
           start = c(1960,3),
           end = c(2019,4),
           frequency = 4)

raw_u <- read_excel('RBA_u.xls') # Loading the unemployment data

u <- raw_u[13:517, # Unemployed persons as percentage of labour force	
           11]     # Dropping the first two observations because we will aggregate the unemployment
                   # rate to quarters and the two obervations in the first quarter 1973 will lead to 
                   # inadequate quarters in the aggregation function. Thus, the two observation will not be 
                   # considered. Otherwise it would have been loaded and the by [-c(11:12), ] 'deloaded.

u <- ts(u,
        start = c(1978,4),
        end = c(2020,4),
        frequency = 12)

u <- aggregate(u, 
               FUN = mean, 
               nfrequency=4)

var.data <- ts.intersect(pi, pi_exp, rgdp, u) 

colnames(var.data) <- c('Pi','Pi_exp','RGDP', 'U') 


(c(sum(is.na(var.data[,'Pi'])),
  sum(is.na(var.data[,'Pi_exp'])),
  sum(is.na(var.data[,'RGDP'])),
  sum(is.na(var.data[,'U']))))

# The dataset is complete, there is no missing value in any of them. 
# If a missing value would be there it would result in a TRUE which is a 1 
# and hence the sum would not be 0

var.data <- window(var.data,
                   start = c(1989,3),
                   frequency = 4)

colnames(var.data) <- c('Pi','Pi_exp','RGDP', 'U') 


# Slicing the data to estimate models for the pseudo-out-of-sample forecast
# The data period to fit the model is 91 quarters long and starts in 1989 Q3 and
# ends in 2012 Q1. Thus covering around 75% of the observations.

var.data.est <- window(var.data,
                       start = c(1989,3),
                       end = c(2012,1),
                       frequency = 4)

colnames(var.data.est) <- c('Pi',
                            'Pi_exp',
                            'RGDP',
                            'U') 

cor(var.data.est) # Checking if the additional predictors correlate with the inflation rate in order to be able as predictors

#       Pi      Pi_exp        RGDP          U
# Pi      1.0000000 -0.30729497  0.19782139 -0.4795721
# Pi_exp -0.3072950  1.00000000 -0.01134665  0.1961288
# RGDP    0.1978214 -0.01134665  1.00000000 -0.0228805
# U      -0.4795721  0.19612878 -0.02288050  1.0000000


# Conclusion: Each variable is well correlated with the inflation rate and thus can be used
# to predict future variability in the inflation rate

coeftest(lm(Pi~RGDP, var.data.est))

#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept) 0.026852   0.443368  0.0606  0.95184  
# RGDP        0.238708   0.125380  1.9039  0.06016 .

# The p-value is 0.06 and the real gdp growth rate correlates with the inflation rate, the effect is unequal zero for confidence
# level of 10%. The model suffers from omitted variable bias, but it is still enough to ensure that real gdp growth rate 
# correlates with inflation rate.

coeftest(lm(Pi~U, var.data.est))

#             Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  3.826422   0.612816  6.2440 1.422e-08 ***
# U           -0.439753   0.085292 -5.1559 1.512e-06 ***

# The unemployment rate has a significant effect on inflation rate. However, the model suffers from OVB and hence, it can
# only be concluded that the unemployment rate correlates with inflation rate

coeftest(lm(Pi~Pi_exp, var.data.est))

#               Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  2.08462    0.46299  4.5025 2.027e-05 ***
#  Pi_exp      -0.64232    0.21084 -3.0464  0.003047 ** 


# The inflation expectation have a significant effect on inflation rate. However, the model suffers from OVB and hence, it can
# only be concluded that the inflation rate expectation correlates with inflation rate

## Conclusion: The chosen variables are suitable predictors for the inflation rate.

############################################################################################
########### 2.1 Data Visualisation #########################################################
############################################################################################
head(var.data.est)

tail(var.data.est)

# The data is going to be visualized to get a first indication of a stationary process and
# to see if the data suffers from extreme values

par(mfrow = c(3,2),
    cex = 0.8,
    mex = 0.8,
    cex.lab= 0.8,
    cex.main =0.8,
    cex.axis = 0.8)

ts.plot(var.data.est[,'Pi'], 
        main = 'Australian Inflation Rate',
        col = 'blue',
        ylab = 'Pi')

abline(h = mean(var.data.est[,'Pi']), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The inflation rate might not be stationary, the variance does not seem to be time
# invariant and it has a trend

# plotting the first difference
ts.plot(diff(var.data.est[,'Pi']), 
        main = '1. Difference Australian Inflation Rate',
        col = 'blue',
        ylab = 'dPi')

abline(h = mean(diff(var.data.est[,'Pi'])), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The first difference does look stationary

ts.plot(var.data.est[,'Pi_exp'], 
        main = 'Australian 3-Months Inflation Expectation',
        col = 'blue',
        ylab = 'Pi_exp')

abline(h = mean(var.data.est[,'Pi_exp']), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The inflation expectation does not look stationary, the variance is 
# clearly time variant and it has a trend downwards.

# plotting the first difference
ts.plot(diff(var.data.est[,'Pi_exp']), 
        main = '1. Difference Australian 3-Months Inflation Expectation',
        col = 'blue',
        ylab = 'dPi_exp')

abline(h = mean(diff(var.data.est[,'Pi_exp'])), 
       col = 'red',
       lty = 'dashed')


# Conclusion: The difference of the inflation expectation does look stationary, with 
# time invariant mean and variance. It will be tested by ADF and KPSS test

ts.plot(var.data.est[,'RGDP'], 
        main = 'Australian Real GDP Growth Rate',
        col = 'blue',
        ylab = 'RGDP')

abline(h = mean(var.data.est[,'RGDP']), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The real gdp growth rate has some higher volatility in 1993/1995.
# Expect the outlier in the beginning of the time series, the variance seems to be 
# more or less time invariant with the exception of the start of the sample. 
# Hence, the process could be stationary and the ADF and KPSS test will be done
# on the real gdp growth rate.

ts.plot(diff(var.data.est[,'RGDP']), 
        main = '1. Difference Real GDP Growth Rate',
        col = 'blue',
        ylab = 'dRGDP')

abline(h = mean(diff(var.data.est[,'RGDP'])), 
       col = 'red',
       lty = 'dashed')

ts.plot(var.data.est[,'U'], 
        main = 'Australian Unemployment Rate',
        col = 'blue',
        ylab = 'U')

abline(h = mean(var.data.est[,'U']), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The unemployment rate clearly does not look stationary and the variance and
# the mean are clearly not time invariant.

ts.plot(diff(var.data.est[,'U']), 
        main = '1. Difference Australian Unemployment Rate',
        col = 'blue',
        ylab = 'dU')

abline(h = mean(diff(var.data.est[,'U'])), 
       col = 'red',
       lty = 'dashed')

# Conclusion: The difference of the unemployment rate does look stationary 
# with a more or less constant variance. The variance has a higher value in 1990 and 2008 
# The differenced unemployment rate will be tested in the ADF & KPSS test.


## Overall conclusion:

# From the plot it seems that:
# Inflation rate has a deterministic trend, because of the positive slope and the increase in levels
# over time but no intercept.
# Inflation rate expectation has a deterministic trend, because of the negative slope and the decrease
# in levels over time and an intercept.
# Real GDP growth rate seems to have just an intercept because it oscillates at a higher level.
# The unemployment rate an intercept because it starts higher than zero and has a deterministic
# trend which leads to a decrease in levels over time.

## These commentaries are important to fit the VAR type later in the paper

############################################################################################
########### 3 Estimation ###################################################################
############################################################################################
########### 3.1 Stationarity ###############################################################
############################################################################################
########### 3.1.1 Augmented Dicky-Fuller Test ##############################################
############################################################################################
# Testing inflation in levels

ADF.test(var.data.est[,'Pi'])


# Results for inflation rate:
#       statistic  1pct  5pct 10pct
# tau3     -2.07 -4.04 -3.45 -3.15   , we CANNOT reject the H0 of a unit root. 
# phi3      2.28  8.73  6.49  5.47   , we CANNOT reject the H0 of a unit root with trend. 
# tau2     -1.23 -3.51 -2.89 -2.58   , we CANNOT reject the H0 of a unit root.
# phi1      1.46  6.70  4.71  3.86   , we CANNOT reject the H0 of a unit root with drift. 
# tau1     -0.22 -2.60 -1.95 -1.61   , we CANNOT reject the H0 of a unit root.

# Conclusion: The differenced inflation will be tested for stationary because the inflation in levels
# is not stationary.

ADF.test(diff(var.data.est[,'Pi']))

# Results for differenced inflation rate:

#       statistic  1pct  5pct 10pct
# tau3     -8.19 -4.04 -3.45 -3.15   , we reject the H0 of a unit root. 
# phi3     33.61  8.73  6.49  5.47   , we reject the H0 of a unit root with trend.
# tau2     -8.24 -3.51 -2.89 -2.58   , we reject the H0 of a unit root. 
# phi1     33.91  6.70  4.71  3.86   , we reject the H0 of a unit root with drift.
# tau1     -8.09 -2.60 -1.95 -1.61   , we reject the H0 of a unit root. 

# Conclusion: The difference inflation rate is an I(0) process and is stationary and inflation
# is thus differenced stationary.


# Testing inflation expectation in levels

ADF.test(var.data.est[,'Pi_exp'])

# Results for  3-month inflation expectation:

#       statistic  1pct  5pct 10pct
# tau3     -3.70 -4.04 -3.45 -3.15   , we reject the H0 of a unit root.
# phi3      6.90  8.73  6.49  5.47   , we reject the H0 of a unit root with trend. 
# tau2     -2.99 -3.51 -2.89 -2.58   , we reject the H0 of a unit root.
# phi1      4.51  6.70  4.71  3.86   , we CANNOT reject the H0 of a unit root with drift
# tau1     -0.83 -2.60 -1.95 -1.61   , we CANNOT reject the H0 of a unit root 

# Conclusion: The 3-month inflation rate expectation is stationary with trend and drift.
# These results will be tested by the KPSS test.


# Testing the real GDP growth rate in levels

ADF.test(var.data.est[,'RGDP'])

# Results for real gdp growth rate:
#       statistic  1pct  5pct 10pct
# tau3     -3.40 -4.04 -3.45 -3.15   , we CANNOT reject the H0 of a unit root.
# phi3      6.89  8.73  6.49  5.47   , we reject the H0 of a unit root with trend
# tau2     -2.66 -3.51 -2.89 -2.58   , we CANNOT reject the H0 of a unit root
# phi1      3.63  6.70  4.71  3.86   , we CANNOT reject the H0 of a unit root with drift
# tau1     -0.03 -2.60 -1.95 -1.61   , we CANNOT reject the H0 of a unit root 


# Conclusion: The real gdp growth rate does not reject the unit root in levels except at phi3. However,
#the differenced real GDP growth rate will be tested becasue phi3 is an F-test and if the real gdp has
# unit root the test theta = b = 0 might not be rejected if real gdp has an trend b <>0. The results
# of tau 2 counter the rejcetion becasue setting the trend to zero, the equaltion rejects the H0 of 
# a unit root in an equation with drift.

ADF.test(diff(var.data.est[,'RGDP']))

# Results for the differenced real gdp growth rate:

#       statistic  1pct  5pct 10pct
# tau3     -4.19 -4.04 -3.45 -3.15   , we reject the H0 of a unit root.
# phi3      9.94  8.73  6.49  5.47   , we reject the H0 of a unit root with trend
# tau2     -4.42 -3.51 -2.89 -2.58   , we reject the H0 of a unit root
# phi1      9.83  6.70  4.71  3.86   , we reject the H0 of a unit root with drift
# tau1     -4.46 -2.60 -1.95 -1.61   , we CANNOT reject the H0 of a unit root 

# Conclusion: The differenced real gdp growth rate is differened stationary. Hence, the process
# is an I(1) process. However, we note that the lags are outside the unit root. The KPSS test will be 
# performed to check these results.


# Testing the unemployment rate in levels

ADF.test(var.data.est[,'U'])

# Results for the unemployment rate:

#       statistic  1pct  5pct 10pct
# tau3     -1.80 -4.04 -3.45 -3.15   , we CANNOT reject the H0 of a unit root
# phi3      2.95  8.73  6.49  5.47   , we CANNOT reject the H0 of a unit root with trend
# tau2     -2.21 -3.51 -2.89 -2.58   , we CANNOT reject the H0 of a unit root
# phi1      4.11  6.70  4.71  3.86   , we CANNOT reject the H0 of a unit root with drift
# tau1     -2.37 -2.60 -1.95 -1.61   , we reject the H0 of a unit root 


## Conclusion: The unemployment rate can only reject the H0 of a unit root without trend
# and drift. The KPSS test will be performed to check the results of the ADF test.


############################################################################################
########### 3.1.2 KPSS Test ################################################################
############################################################################################
summary(ur.kpss(diff(var.data.est[,'Pi']))) 

# Results for inflation rate:
# The test statistic is 0.0476  < 0.463 The H0 of stationary is not rejected.
# Conclusion: The results are in line with ADF test

summary(ur.kpss(var.data.est[,'Pi_exp']))

# Results for inflation expectation:
# The test statistic is 0.8909 > 0.463. The H0 of stationary is rejected.
# Conclusion: The results are not in line with ADF test

ADF.test(diff(var.data.est[,'Pi_exp'])) # Perform ADF on the differenced inflation exp. rate

#       statistic  1pct  5pct 10pct
# tau3     -5.36 -4.04 -3.45 -3.15   , we CANNOT reject the H0 of a unit root
# phi3     14.38  8.73  6.49  5.47   , we CANNOT reject the H0 of a unit root with trend
# tau2     -5.36 -3.51 -2.89 -2.58   , we CANNOT reject the H0 of a unit root
# phi1     14.37  6.70  4.71  3.86   , we CANNOT reject the H0 of a unit root with drift
# tau1     -5.39 -2.60 -1.95 -1.61   , we reject the H0 of a unit root 

summary(ur.kpss(diff(var.data.est[,'Pi_exp'])))
# The test statistic is 0.1283 < 0.463. The H0 of stationary is not rejected.
# Conclusion: The results are in line with ADF test.


summary(ur.kpss(diff(var.data.est[,'RGDP'])))

# Results for differenced real gdp growth rate:
# The test statistic is 0.0673 < 0.463. The H0 of stationary is not rejected.
# Conclusion: The results are in line with ADF test

summary(ur.kpss(var.data.est[,'U']))

# Results for unemployment rate:
# The test statistic is 1.7642 > 0.463. The H0 of stationary is rejected.
# Conclusion: The results are not in line with ADF test. 

ADF.test(diff(var.data.est[,'U'])) # Perform ADF on the differenced unemployment rate

# Results for differenced unemployment rate:

#       statistic  1pct  5pct 10pct
# tau3     -5.15 -4.04 -3.45 -3.15   , we reject the H0 of a unit root
# phi3     13.40  8.73  6.49  5.47   , we reject the H0 of a unit root with trend
# tau2     -4.69 -3.51 -2.89 -2.58   , we reject the H0 of a unit root
# phi1     11.02  6.70  4.71  3.86   , we reject the H0 of a unit root with drift
# tau1     -4.23 -2.60 -1.95 -1.61   , we reject the H0 of a unit root 


# Conclusion: The differenced unemployment rate is a stationary process of I(0) according to 
# the ADF test. The results will be challenged by the KPSS test.

summary(ur.kpss(diff(var.data.est[,'U']))) # Perform KPSS test on differenced unemployment rate

# Results for unemployment rate:
# The test statistic is 0.2652  < 0.463. The H0 of stationary not is rejected.
# Conclusion: The results are in line with ADF test.

## Overall conclusion:

# The inflation rate is I(1)
# The inflation expectation is I(1) 
# The real gdp growth rate is I(1)
# The unemployment rate is I(1)

############################################################################################
########### 3.2 ARIMA(p,d,q) model #########################################################
############################################################################################
########### 3.2.1 Identfication ############################################################
############################################################################################
acf(diff(var.data.est[,'Pi']), 
    col = 'blue',
    main = 'Autocorrelogram')

# Result ACF:
# General: The ACF fades away

pacf(diff(var.data.est[,'Pi']), 
     col = 'blue',
     main = 'Partial Autocorrelogram')

# PACF
# General: The PACF cuts off at lag 1. It has an insignificant spike between lag 1 and 2.

# Conclusion: The differenced inflation rate seems to be an autoregressive model of lag order 1
# because of the fading away of the acf and the cut off in the pacf. Two models will be fitted
# and tested

ar1.result <- arima(diff(var.data.est[,'Pi']), 
                    order = c(1, 0, 0))

ar2.result <- arima(diff(var.data.est[,'Pi']), 
                    order = c(2, 0, 0))

############################################################################################
########### 3.2.2 Information Criteria #####################################################
############################################################################################

AIC(ar1.result, 
    ar2.result)

#             df      AIC
# ar1.result  3 268.2143
# ar2.result  4 269.3297

# The AIC favours the AR(1,0,0) model with the lowest value among the tested models.

BIC(ar1.result, 
    ar2.result)

#            df      BIC
# ar1.result  3 275.7137
# ar2.result  4 279.3290

# The BIC favours the AR(1,0,0) model with the lowest value among the tested models.

# Conclusion: Both information criteria favour AR(1) as a model. This means inflation is 
# an ARIMA(1,1,0) model.


############################################################################################
########### 3.2.3 Box-Jenkins-Procedure ####################################################
############################################################################################

LjungBox(ar1.result, 
         lags = 1:10)

# lags statistic df   p-value
#  1 0.2889103  0 0.0000000
#  2 0.4543636  1 0.5002700
#  3 0.7354917  2 0.6922931
#  4 1.2993370  3 0.7292905
#  5 3.7365284  4 0.4428354
#  6 3.9105934  5 0.5623589
#  7 4.0709470  6 0.6670758
#  8 5.1127047  7 0.6462115
#  9 7.7473389  8 0.4585320
# 10 7.9004848  9 0.5442053

# Conclusion: For each lag the p-value is above 5%.  Hence, the H0 of no autocorrelation 
# cannot be rejected.The ARIMA(1,1,0) model does not suffer from autocorrelation.

jarque.bera.test(residuals(ar1.result))

# Conclusion: The p-value is  0.3876 significantly above 5%.
# Hence, the H0 of normally distributed residuals is NOT rejected. 
# Thus, the residuals are normally distributed.

############################################################################################
########### 3.4 VAR 1 Inflation rate and real gdp growth rate ##############################
############################################################################################
########### 3.4.1 Identification ##########################################################
############################################################################################

# Pi is I(1) an RGDP is I(1).
# Thus, these variables have to be tested for cointegration meaning there could be a
# lineare combination of those that is stationary. The test used for it is the Engle-Granger
# test.

var1.data.est <-var.data.est[, c('Pi', 'RGDP')]  # data to estimate the model
                             
coint <- lm(RGDP ~ Pi,            
            data = var1.data.est)

# In the process of writing this script, doing Pi ~ RGDP has led to no cointegration 
# thereby leading to biased estimation of the forecasts. The forecasts have been in a 
# range of 30 basis points which means they were very low in magnitude.
# The forecasts were almost zero. Even after changing the type and the lag order 
# did not alleviate this problem. However, I noticed that changing the regression 
# order led to cointegration. Performing the forecast in that matter led to 
# more reasonable forecasts. Hence, I continue with the assumption of cointegrated variables. 

ADF.test(residuals(coint))

#      statistic  1pct  5pct 10pct
# tau2     -4.06 -3.51 -2.89 -2.58

# Cointegration
# Results of the Engle Granger test:
# The test statistic is -4.06 < - 3.41 (critival value for one regressor according to S&W table)
# The null of a unit root in the residuals is rejected and thus the variables are cointegrated.
# A cointegrated relationship might also be economically meaningful due to the business cycle.
# If the economy is in a boom, the prices usually tends to increase and if the economy is contracting
# the prices tends to decrease.

# Conclusion: The inflation rate and the real GDP growth rate are cointegrated and a VAR model 
# in levels is fitted
  
# From the visualization of the processes in levels, the inflation rate has a trend.
# The real gdp growth rate has an intercept. 
# Hence, type 'const' is specified because differencing takes away the trend an the 
# cointegration has the same effect.

VARselect(var1.data.est, # Choosing the lag order for the VAR model
          type = 'const')$selection 


# Results of the selection:
# AIC(n)  HQ(n)  SC(n) FPE(n) 
#  6      3      2      6 

# Conclusion: 
# While AIC tends to take a lag order which is often the true lag order of the 
# data generating process, the BIC chooses a lag order that does better short-term
# forecasts. Furthermore, the BIC punishes more complex models. The purpose of the paper is 
# to assess forecasting precision. Hence, the lag order of 2 is chosen according to BIC.

var1.result <- VAR(var1.data.est, #Fitting the model VAR(2)
                   p = 2,
                   type ='const')

coeftest(var1.result$varresult$Pi)

#           Estimate Std. Error t value  Pr(>|t|)    
# Pi.l1    0.412445   0.099727  4.1357 8.350e-05 ***
# RGDP.l1  0.129807   0.122755  1.0574    0.2933    
# Pi.l2    0.509300   0.109113  4.6677 1.142e-05 ***
# RGDP.l2 -0.117680   0.120185 -0.9792    0.3303    
# const    0.108164   0.273714  0.3952    0.6937 

# Conclusion: The effect of the first lage and second lag in the model are clearly different from
# zero while for the other coefficients the effect is questionable 

############################################################################################
########### 3.4.2 Box-Jenkins-procedure ####################################################
############################################################################################
serial.test(var1.result) # Perform multivariate Ljung Box test

# data:  Residuals of VAR object var1.result
# Chi-squared = 63.052, df = 56, p-value = 0.2411

# Conclusion: the p-value is 0.2411 and thus we cannot reject the null of no autocorrelation.
# Hencce, the model is well specified. 

normality.test(var1.result)$jb.mul$JB # Perform multivariate Jarque-Bera Test

# data:  Residuals of VAR object var1.result
# Chi-squared = 0.92905, df = 4, p-value = 0.9204

# Conclusion: The p value is 0.9204 and thus the null of normally distributed residuals
# is NOT rejected. The residuals are normally distributed.


roots(var1.result) # Testing the stability of the VAR model
# [1] 0.9504595 0.7326898 0.5606058 0.1320653

# Conclusion: The eigenvalues are real and distinct. The model is stable because the eigenvalues are 
# inside the unit interval


############################################################################################
########### 3.5 VAR 2 Inflation rate, real gdp growth rate and unemployment rate ###########
############################################################################################
########### 3.5.1 Identification ###########################################################
############################################################################################

# Pi is I(1), unemployment rate is I(1) and RGDP is I(1)
# The three variables could share a stochastic trend and thus a linear combination of these
# might exist that represents a long term relationship and is stationary.
# In order to assess wether there exists a cointegration, the Johansen Procedure is performed because three
# variable are tested for cointegration and it permits more than one cointegrated relationship,
# which the Engle Granger test does not.

var2.data.est <- var.data.est[,c('Pi','RGDP','U')] # data to estimate the model

JJ.test(var2.data.est, test = "trace") # Perform the Johansen Trace test

#              tstat none cval none tstat const cval const tstat trend cval trend
# r = 0  |      37.52     31.52       38.21      34.91       72.74      42.44
# r <= 1 |       9.43     17.95       10.11      19.96       35.30      25.32
# r <= 2 |       0.64      8.18        1.28       9.24        7.47      12.25

# Results of the Johansen trace test
# Lags are selected by AIC, here the lags is 2
# For r = 0 the H0 of no cointegration is rejected for none, const and trend
# For r <= 1 the H0 is not rejected for none, const and trend. Hence, r = 1

JJ.test(var2.data.est, test = "eigen") # Perform the Johansen–Juselius’ maximum eigenvalue test

#             tstat none cval none tstat const cval const tstat trend cval trend
# r = 0  |      28.09     21.07       28.10      22.00       37.44      25.54
# r <= 1 |       8.79     14.90        8.83      15.67       27.83      18.96
# r <= 2 |       0.64      8.18        1.28       9.24        7.47      12.25


# Results of the Johansen–Juselius’ maximum eigenvalue test trace test
# Lags are selected by AIC, here the lags is 2
# For r = 0 the H0 of no cointegration is rejected for none, const and trend
# For r <= 1 the H0 is not rejected for none, const and trend. Hence, r = 1

## Conclusion: Both test one cointegrated relationship

# From the visualisation, the concluded the following specifications of the processes:
# The inflation rate has a trend. The trend was removed due to differencing
# The real gdp growth rate has an intercept. 
# The unemployment rate has an intercept and trend. The trend was removed due to differencing
# The ecdet 'const' will be specified because the intercept could remain and give the model
# stability.

var2.data.test.t <- ca.jo(var2.data.est, 
                          type ="trace",
                          ecdet = 'const',
                          K = 2)         # Number of lags. Here, K = 2

vecm2 <- cajorls(var2.data.test.t, 
                 r = 1)                  # Number of cointegrated relationship. Here, r = 1

vecm2$rlm # VECM model has two lags and one cointegrated relationship
          # Looking at the model equation.

# Coefficients:
#               Pi.d       RGDP.d     U.d      
# ect1       0.002726   0.019686   0.002982
# Pi.dl1    -0.562691   0.122934  -0.009125
# RGDP.dl1   0.080768  -0.245268  -0.092597
# U.dl1     -0.313904  -0.800923   0.418216


var2.result <- vec2var(var2.data.test.t, # Transform the VECM back to VAR
                        r = 1)

############################################################################################
########### 3.5.2 Box-Jenkins-procedure ####################################################
############################################################################################

serial.test(var2.result) # Perform multivariate Ljung Box test

# data:  Residuals of VAR object var2.result
# Chi-squared = 116.93, df = 129, p-value = 0.7686

# Conclusion: the p-value is 0.7686 and thus we cannot reject the null of no autocorrelation.
# This is good because the model does not suffer from autocorrelation

normality.test(var2.result)$jb.mul$JB # Perform multivariate Jarque-Bera Test

# data:  Residuals of VAR object var2.result
# Chi-squared = 1.9158, df = 6, p-value = 0.9273

# Conclusion: The p value is 0.9273 and thus the null of normally distributed residuals
# is NOT rejected. Hence, the residuals are normally distributed.

############################################################################################
########### 3.6 VAR 3 Infl. rate, real gdp growth rate, unemployment rate and infl. expectation
############################################################################################
########### 3.6.1 Identification ###########################################################
############################################################################################

# Pi is I(1), Unemployment rate is I(1), RGDP is I(1), and Pi_exp is I(1) are integrated of
# the same order. The four variables could share a stochastic trend and thus a linear 
# combination of these might exist that represents a long term relationship and is stationary.
# In order to assess wether there exists a cointegration, the Johansen Procedure is performed because three
# variable are tested for cointegration and it permits more than one cointegrated relationship,
# which the Engle Granger test does not.

var3.data.est <- var.data.est[, c('Pi', 'RGDP', 'U', 'Pi_exp')]


JJ.test(var3.data.est, test = "trace") # Perform the Johansen Procedure

#             tstat none cval none tstat const cval const tstat trend cval trend
# r = 0  |      65.49     48.28       67.53      53.12       91.41      62.99
# r <= 1 |      32.22     31.52       34.24      34.91       54.39      42.44
# r <= 2 |       7.52     17.95        9.55      19.96       22.76      25.32
# r <= 3 |       0.01      8.18        1.88       9.24        7.16      12.25


# Results for the Johansen Test:
# For r = 0 the H0 of no cointegration is rejected for none, const and trend
# For r <= 1 the H0 is rejected for none, but not for const and trend.  Hence r = 1 with lag 2
# For r <= 2 the H0 is NOT rejected for none, const and trend.
# For r <= 3 the H0 is NOT rejected, but r = 1 because it is an algorithm and the estimation
# stops after r <= 1 does not reject the H0


JJ.test(var3.data.est, test = "eigen") # Perform the Johansen–Juselius’ maximum eigenvalue test

#               tstat none cval none tstat const cval const tstat trend cval trend
# r = 0  |      33.27     27.14       33.28      28.14       37.02      31.46
# r <= 1 |      24.70     21.07       24.70      22.00       31.62      25.54
# r <= 2 |       7.51     14.90        7.67      15.67       15.60      18.96
# r <= 3 |       0.01      8.18        1.88       9.24        7.16      12.25

# Results for the Johansen-Juselius  maximum eigenvalue test:
# For r = 0 the H0 of no cointegration is rejected for none, const and trend
# For r <= 1 the H0 is rejected for none, const and trend.
# For r <= 2 the H0 is NOT rejected for none, const and trend. Hence r = 2 with lag 2
# For r <= 3 the H0 is NOT rejected, but r = 2 because it is an algorithm and the estimation

## Conclusion: The Johansen trace test detects one cointegrated relationship while the Johansen-Juselus'
# maximum eigenvalue two cointegrated relationships detected. I will continue with the results of the 
# trace test because Johansen and Juselius (1990) recommended to use the trace test results if the results
# between the two statistics differ. Thus r = 1 and we use a const and two lags.


var.data3.test.t <- ca.jo(var3.data.est, 
                          type ="trace",   
                          ecdet = 'const', # including a const in the model
                          K = 2)           # two lags

vecm3 <- cajorls(var.data3.test.t,         # estimating the VECM
                 r = 1)

vecm3$rlm # The model is a VECM  with 2 lags and 2 cointegrated relationships

# Coefficients:
#             Pi.d        RGDP.d      U.d         Pi_exp.d  
# ect1         0.0038600   0.0147233  -0.0006999   0.0032094
# Pi.dl1      -0.5683306   0.0771445   0.0049735  -0.0379962
# RGDP.dl1     0.0271420  -0.3559901  -0.0370523  -0.0873184
# U.dl1       -0.4254167  -0.5156898   0.5815385  -0.8066495
# Pi_exp.dl1  -0.1763757  -0.2441853  -0.1174512  -0.3999980

var3.result <- vec2var(var.data3.test.t, # Transform the VECM back to VAR
                       r = 1)

var3.result

############################################################################################
########### 3.6.2 Box-Jenkins-procedure ####################################################
############################################################################################
serial.test(var3.result) # Perform multivariate Ljung Box test

# data:  Residuals of VAR object var3.result
# Chi-squared = 204.11, df = 228, p-value = 0.8706

# Conclusion: The p-value is 0.8706 and thus the null of no autocorrelation is NOT rejected.
# Hence, the model does not suffer from autocorrelation.

normality.test(var3.result)$jb.mul$JB # Perform multivariate Jarque-Bera Test

# data:  Residuals of VAR object var3.result
# Chi-squared = 9.5175, df = 8, p-value = 0.3005

# Conclusion: The p value is 0.3005 and thus the null of normally distributed residuals
# is not rejected.  Hence, the residuals are normally distributed

############################################################################################
########### 3.7 Pseudo-out-of-sample forecasts #############################################
############################################################################################
########### 3.7.1 Pseudo-out-of-sample rekursive forecast with horizon 2, 4 and 8 ##########
############################################################################################
########### 3.7.1.1 AR1 ####################################################################
############################################################################################

# The model for the Australian inflation rate is ARIMA(1,1,0). In order to adjust the model 
# at each origin the loop will make use of the auto.arima function. However, the paper underscores 
# that the model is an ARIMA(1,1,0) and inflation behaves like that.


AR.RE <- ts(matrix(NA,
                  nrow = 30,
                  ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(AR.RE) <- paste("Horizon", 1:8, sep = " ")

for(i in 1:30){                                               
  EndDate       <- 2012 + (i - 1) / 4                   
  
  data          <- window(var.data[,'Pi'],
                          end = EndDate)                    
  
  model        <- auto.arima(data, ic = 'bic') # BIC does asymptotically better short term forecasts                       
   
  AR.RE[i, ] <- forecast(model,
                         h = 8)$mean            
} # Forecasting the inflation rate with the ARIMA model

head(AR.RE, n = 30)
# Conclusion: The results look fine and the loop seems to have worked correctly

############################################################################################
########### 3.7.1.2 VAR 1 Inflation rate and real gdp growth rate  #########################
############################################################################################

# From the identification of the VAR1,we fitted the a VAR(2) in levels

# data to evaluate and fit the forecasting model
var1.data <- var.data[,c('Pi', 'RGDP')]
colnames(var1.data) <- c('Pi','RGDP')


VAR1.RE <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                     nrow = 30,
                     ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(VAR1.RE) <- paste("Horizon", 
                          1:8, 
                          sep = " ")

for(i in 1:30){                                               
  EndDate       <- 2012 + (i - 1) / 4                   
  
  data          <- window(var1.data,
                          end = EndDate)                    
  
  model        <- VAR(data,
                      p = 2,
                      type = 'const')
  
  VAR1.RE[i, ] <- predict(model,
                          n.ahead = 8)$fcst$Pi[,"fcst"]       
} # Forecasting the inflation rate with the VAR1 model

head(VAR1.RE, n = 30)
# Conclusion: The loop seems to have worked correctly. 


############################################################################################
########### 3.7.1.3 VAR 2 Inflation rate, real gdp growth rate and unemployment rate #######
############################################################################################

# For VAR(2) we fittted an VECM with lag order 2 and 1 cointegrated relationship

var2.data <-var.data[,c('Pi', 'RGDP', 'U')]

VAR2.RE <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                     nrow = 30,
                     ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(VAR2.RE) <- paste("Horizon", 1:8, sep = " ")

for(i in 1:30){                                               
  
  EndDate       <- 2012 + (i - 1) / 4                   
  
  data          <- window(var2.data,
                          end = EndDate)                    
  
  var.data2.test.t <- ca.jo(data, 
                            type ="trace", 
                            ecdet = 'const',
                            K = 2)
  
  vecm2 <- cajorls(var.data2.test.t, 
                   r = 1)
  
  var2.result <- vec2var(var.data2.test.t, # Transform the VECM back to VAR
                         r = 1)
  
  VAR2.RE[i, ] <- predict(var2.result,
                          n.ahead = 8)$fcst$Pi[,"fcst"]          
} # Forecasting the inflation rate with the VAR2 model

head(VAR2.RE, n = 30)
# Conclusion: The loop seems to have worked correctly

############################################################################################
########### 3.7.1.4 VAR 3 Infl rate, real gdp growth rate, unemployment rate and infl expectation 
############################################################################################

# For VAR 3 an VECM of lag order 2 with 1 cointegrated relationships.

var3.data <- var.data[, c('Pi','RGDP','U','Pi_exp')]


VAR3.RE <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                     nrow = 30,
                     ncol = 8),
                start = c(2012,2),
                frequency = 4)

colnames(VAR3.RE) <- paste("Horizon", 1:8, sep = " ")


for(i in 1:30){                                               
  EndDate       <- 2012  + (i - 1) / 4                   
  
  data          <- window(var3.data, end = EndDate)                  
  
  var.data3.test.t <- ca.jo(data, 
                            type ="trace",
                            ecdet = "const", 
                            K = 2)
  
  vecm3 <- cajorls(var.data3.test.t, r = 1)
  
  var3.result <- vec2var(var.data3.test.t, r = 1) # Transform the VECM back to VAR
                         
  VAR3.RE[i, ] <- predict(var3.result, n.ahead = 8)$fcst$Pi[,"fcst"]          
} # Forecasting the inflation rate with the VAR3 model

head(VAR3.RE, n = 30)
# Conclusion: The loop seems to have worked correctly

############################################################################################
########### 3.7.2 Pseudo-out-of-sample rolling forecast with horizon 2, 4 and 8 ############
############################################################################################
############################################################################################
########### 3.7.2.1 AR1  ###################################################################
############################################################################################

# The model for the Australian inflation rate is ARIMA(1,1,0). In order to adjust the model at each origin
# the loop will make use of the auto.arima function. However, the paper underscores that the 
# model is an AR(1) and inflation behaves like that.

AR.RO <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                  nrow = 30,
                  ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(AR.RO) <- paste("Horizon", 
                        1:8, 
                        sep = " ")

for(i in 1:30){
  
  StartDate <- 1989 + 2/4 + (i-1) / 4
  EndDate <- 2012 + (i-1) / 4
  
  data <- window(var.data[,'Pi'],
                 start = StartDate,
                 end = EndDate)
  
  model <- auto.arima(data, ic ='bic')
  
  AR.RO[i,] <- forecast(model,
                        h = 8)$mean
} # Forecasting the inflation rate with the ARIMA model

head(AR.RO, n =30)
# Conclusion: The loop seems to have worked correctly

############################################################################################
########### 3.7.2.2 VAR 1 Inflation rate and real gdp growth rate ##########################
############################################################################################

# The model fitted is a VAR(2) in levels according to the identification

VAR1.RO <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                     nrow = 30,
                     ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(VAR1.RO) <- paste("Horizon", 1:8, sep = " ")


for(i in 1:30){
  
  StartDate <- 1989 + 2/4 + (i - 1) / 4 
  EndDate <- 2012 + (i - 1) / 4          
  
  data <- window(var1.data,
                 start = StartDate,
                 end = EndDate)
  
  model <- VAR(data,
               p = 2,
               type = 'const')
  
  VAR1.RO[i,] <- predict(model,
                        n.ahead = 8)$fcst$Pi[,"fcst"]
} # Forecasting the inflation rate with the VAR1 model

head(VAR1.RO, n = 30)
# Conclusion: The loop seems to have worked

############################################################################################
########### 3.7.2.3 VAR 2 Inflation rate, real gdp growth rate and unemployment rate #######
############################################################################################

# For VAR(2) we fittted a VECM with lag order 2 and 1 cointegrated relationship

VAR2.RO <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                    nrow = 30,
                    ncol = 8),
                start = c(2012,2),
                frequency = 4)

colnames(VAR2.RO) <- paste("Horizon", 
                          1:8, 
                          sep = " ")

for(i in 1:30){
  
  StartDate <- 1989 + 2/4 + (i-1) / 4
  EndDate <- 2012 + (i-1) / 4
  
  data <- window(var2.data,
                 start = StartDate,
                 end = EndDate)
  
  var.data2.test.t <- ca.jo(data, 
                            type ="trace",
                            ecdet = 'const',
                            K = 2)
  
  vecm2 <- cajorls(var.data2.test.t, 
                   r = 1)
  
  var2.result <- vec2var(var.data2.test.t, 
                         r = 1)
  
  VAR2.RO[i,] <- predict(var2.result,
                         n.ahead = 8)$fcst$Pi[,"fcst"]
} # Forecasting the inflation rate with the VAR2 model

head(VAR2.RO, n = 30)
# Conclusion: The loop seems to have worked

############################################################################################
########### 3.7.2.4 VAR 3 Infl rate, real gdp growth rate, unemployment rate and infl expectation
############################################################################################

# The model is a VECM  with 2 lags and 1 cointegrated relationships

VAR3.RO <- ts(matrix(NA, # Creating an empty shell to fill with the forecast
                    nrow = 30,
                    ncol = 8),
              start = c(2012,2),
              frequency = 4)

colnames(VAR3.RO) <- paste("Horizon", 1:8, sep = " ")


for(i in 1:30){
  
  StartDate <- 1989 + 2/4 + (i-1) / 4
  EndDate <- 2012 + (i-1) / 4
  
  data <- window(var3.data,
                 start = StartDate,
                 end = EndDate)
  
  var.data3.test.t <- ca.jo(data,
                            type = 'trace',
                            ecdet = 'const',
                            K = 2)
  
  vecm3 <- cajorls(var.data3.test.t, 
                   r = 1)
  
  var3.result <- vec2var(var.data3.test.t, # Transform the VECM back to VAR
                         r = 1)
  
  VAR3.RO[i,] <- predict(var3.result,
                        n.ahead = 8)$fcst$Pi[,"fcst"]
} # Forecasting the inflation rate with the VAR3 model

head(VAR3.RO, n = 30)
# Conclusion: The loop seems to have worked.

############################################################################################
########### 3.7.3 Forecast Errors ##########################################################
############################################################################################
########### 3.7.3.1 Outcome ################################################################
############################################################################################
OUTCOME.2 <- ts(matrix(NA, # Creating an empty shell to fill the outcome with for horizon 2
                       nrow = 30,
                       ncol = 2),
                start =c(2012,2),
                frequency = 4)

colnames(OUTCOME.2) <- paste("Horizon", 1:2, sep = " ")

for(i in 1:30){                                               
  StartDate <- 2012 + 1/4 + (i - 1) / 4           
  EndDate  <- StartDate + 1/4
  
  
  OUTCOME.2[i, ] <- window(var.data[,'Pi'],
                           start = StartDate,
                           end = EndDate)
  
}  # Creating the outcome for forecast horizon of 2

head(OUTCOME.2)
# Conclusion: I checked the results with the data and it worked correctly


OUTCOME.4 <- ts(matrix(NA, # Creating an empty shell to fill the outcome with for horizon 4
                       nrow = 28,
                       ncol = 4),
                start = c(2012,2),
                frequency = 4)

colnames(OUTCOME.4) <- paste("Horizon", 1:4, sep = " ")

for(i in 1:28){                                               
  StartDate <- 2012 + 1/4 + (i - 1) / 4           
  EndDate  <- StartDate + 3/4
  
  
  OUTCOME.4[i, ] <- window(var.data[,'Pi'],
                           start = StartDate,
                           end = EndDate)
} # Creating the outcome for forecast horizon of 4

head(OUTCOME.4)
# Conclusion: Worked correctly

OUTCOME.8 <- ts(matrix(NA, # Creating an empty shell to fill the outcome with for horizon 8
                       nrow = 24,
                       ncol = 8),
                start = c(2012,2),
                frequency = 4)

colnames(OUTCOME.8) <- paste("Horizon", 1:8, sep = " ")

for(i in 1:24){                                               
  StartDate <- 2012 + 1/4 + (i - 1) / 4           
  EndDate  <- StartDate + 7/4
  
  
  OUTCOME.8[i, ] <- window(var.data[,'Pi'],
                           start = StartDate,
                           end = EndDate)
}  # Creating the outcome for forecast horizon of 8

head(OUTCOME.8)
# Conclusion: Worked correctly

# Creating the Forecast errors

############################################################################################
########### 3.7.3.2 Forecast Errors Recursive ##############################################
############################################################################################
## Creating the forecast errors as forecast - actual

# Forecast errors of the ARIMA (1,1,0)

AR.RE.2.FE <- AR.RE[,1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
AR.RE.4.FE <- AR.RE[,1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the first four columns
AR.RE.8.FE <- AR.RE[,1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

# VAR model 1

VAR1.RE.2.FE <-  VAR1.RE[,1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
VAR1.RE.4.FE <-  VAR1.RE[,1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the first four columns
VAR1.RE.8.FE <- VAR1.RE[, 1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

# VAR model 2

VAR2.RE.2.FE <- VAR2.RE[, 1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
VAR2.RE.4.FE <- VAR2.RE[, 1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the four two columns
VAR2.RE.8.FE <- VAR2.RE[, 1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

# VAR model 3

VAR3.RE.2.FE <- VAR3.RE[, 1:2] - OUTCOME.2 # 1:2 because h=2 and we need only the first two columns
VAR3.RE.4.FE <- VAR3.RE[, 1:4] - OUTCOME.4 # 1:4 because h=4 and we need only the four two columns
VAR3.RE.8.FE <- VAR3.RE[, 1:8] - OUTCOME.8 # 1:8 because h=8 and we need the eight two columns

############################################################################################
########### 3.7.3.3 Forecast Errors Rolling ################################################
############################################################################################
# The forecasts errors are calculated as forecast - actual
# Forecast errors of the ARIMA (1,1,0)

AR.RO.2.FE <- AR.RO[,1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
AR.RO.4.FE <- AR.RO[,1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the first four columns
AR.RO.8.FE <- AR.RO[,1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

# Forecast errors of the VAR1

VAR1.RO.2.FE <-  VAR1.RO[,1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
VAR1.RO.4.FE <-  VAR1.RO[,1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the first four columns
VAR1.RO.8.FE <- VAR1.RO[, 1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

# Forecast errors of the VAR2

VAR2.RO.2.FE <- VAR2.RO[, 1:2] - OUTCOME.2  # 1:2 because h=2 and we need only the first two columns
VAR2.RO.4.FE <- VAR2.RO[, 1:4] - OUTCOME.4  # 1:4 because h=4 and we need only the first four columns
VAR2.RO.8.FE <- VAR2.RO[, 1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns 

# Forecast errors of the VAR3

VAR3.RO.2.FE <- VAR3.RO[, 1:2] - OUTCOME.2
VAR3.RO.4.FE <- VAR3.RO[, 1:4] - OUTCOME.4
VAR3.RO.8.FE <- VAR3.RO[, 1:8] - OUTCOME.8  # 1:8 because h=8 and we need the eight two columns

############################################################################################
########### 4. Forecast evaluation #########################################################
############################################################################################
########### 4.1 Plotting the forecasts #####################################################
############################################################################################

par(mfrow = c(2,2),
    cex = 0.8,
    mex = 0.8,
    cex.lab = 0.8,
    cex.main = 0.8,
    cex.axis = 0.8)

### Horizon 2

# Evaluation period
var.data.eval <- window(var.data,
                        start = c(2012,1),
                        end = c(2019,4),
                        frequency = 4)

colnames(var.data.eval) <- c('Pi','Pi_exp','RGDP', 'U') 

plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Recursive Forecast with Horizon 2',
     col = 'grey')

lines(AR.RE[, 'Horizon 2'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RE[, 'Horizon 2'],
      col = 'blue',
      lty = 2)

lines(VAR2.RE[, 'Horizon 2'],
      col = 'red',
      lty = 3)

lines(VAR3.RE[, 'Horizon 2'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex= 0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE) 

plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Rolling Forecast with Horizon 2',
     col = 'grey')

lines(AR.RO[, 'Horizon 2'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RO[, 'Horizon 2'],
      col = 'blue',
      lty = 2)

lines(VAR2.RO[, 'Horizon 2'],
      col = 'red',
      lty = 3)

lines(VAR3.RO[, 'Horizon 2'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex=0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE)

### Horizon 4

plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Recursive Forecast with Horizon 4',
     col = 'grey')

lines(AR.RE[, 'Horizon 4'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RE[, 'Horizon 4'],
      col = 'blue',
      lty = 2)

lines(VAR2.RE[, 'Horizon 4'],
      col = 'red',
      lty = 3)

lines(VAR3.RE[, 'Horizon 4'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex=0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE)


plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Rolling Forecast with Horizon 4',
     col = 'grey')

lines(AR.RO[, 'Horizon 4'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RO[, 'Horizon 4'],
      col = 'blue',
      lty = 2)

lines(VAR2.RO[, 'Horizon 4'],
      col = 'red',
      lty = 3)

lines(VAR3.RO[, 'Horizon 4'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex=0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE)

### Horizon 8

plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Recursive Forecast with Horizon 8',
     col = 'grey')

lines(AR.RE[, 'Horizon 8'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RE[, 'Horizon 8'],
      col = 'blue',
      lty = 2)

lines(VAR2.RE[, 'Horizon 8'],
      col = 'red',
      lty = 3)

lines(VAR3.RE[, 'Horizon 8'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex=0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE)

plot(var.data.eval[, 'Pi'],
     xlab = ' Time',
     ylab = 'Pi',
     main = 'Rolling Forecast with Horizon 8',
     col = 'grey')

lines(AR.RO[, 'Horizon 8'],
      col = 'yellow',
      lty = 1) 

lines(VAR1.RO[, 'Horizon 8'],
      col = 'blue',
      lty = 2)

lines(VAR2.RO[, 'Horizon 8'],
      col = 'red',
      lty = 3)

lines(VAR3.RO[, 'Horizon 8'],
      col = 'green',
      lty= 4)

legend('topright', 
       legend=c('Inflation',"AR 1", "VAR 1", 'VAR 2', 'VAR 3'),
       col = c('grey',"yellow", "blue", 'red', 'green'),
       lty = c(1,1,2,3,4),
       cex=0.5,
       box.lty=0, 
       box.lwd=1, 
       bg="transparent",
       horiz=FALSE)

# Conclusion:
# The plot seem that each model is performing their job to forecast the inflation rate
# The VAR2, VAR3 and ARIMA seem to perform equally good. The VAR1, however, does
# not as goood as the other models. It seems to be biased because it tends to underestimate
# the actual inflation rate. It has also less movement in comparison to the other models.
# This could hint that the model might be biased.

############################################################################################
########### 4.2 Bias Test  #################################################################
############################################################################################
########### 4.2.1 Test statistics  #########################################################
############################################################################################

# Testing the forecast for biasedness
# The H0 is that the forecasts are unbiased. Hence, not rejecting the null is
# favourable for the forecast because then the forecast error is the expected 
# forecast error (e=0) and it has no intercept a=0, which would lead to systematically
# over- or underestimation of the inflation rate.
# H0:  E[FE] = a + e , a=0 and hence the expected error e = 0

#### Horizon 2  

AR.RE.2p <- apply(t(1:2), # Perform the bias test
                  2, 
                  bias.test, 
                  AR.RE.2.FE)
names(AR.RE.2p) <- 1:2

VAR1.RE.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR1.RE.2.FE)
names(VAR1.RE.2p) <- 1:2

VAR2.RE.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR2.RE.2.FE)
names(VAR2.RE.2p) <- 1:2

VAR3.RE.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR3.RE.2.FE)
names(VAR3.RE.2p) <- 1:2

AR.RO.2p <- apply(t(1:2),
                  2, 
                  bias.test, 
                  AR.RO.2.FE)
names(AR.RO.2p) <- 1:2

VAR1.RO.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR1.RO.2.FE)
names(VAR1.RO.2p) <- 1:2

VAR2.RO.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR2.RO.2.FE)
names(VAR2.RO.2p) <- 1:2

VAR3.RO.2p <- apply(t(1:2),
                    2, 
                    bias.test, 
                    VAR3.RO.2.FE)
names(VAR3.RO.2p) <- 1:2

bias2 <- rbind(AR.RE.2p, # binding all test results in one matrix
              AR.RO.2p,
              VAR1.RE.2p,
              VAR1.RO.2p,
              VAR2.RE.2p,
              VAR2.RO.2p,
              VAR3.RE.2p,
              VAR3.RO.2p)

colnames(bias2) <- paste('Horizon', 1:2, sep = ' ')

### Horizon 4.

AR.RE.4p <- apply(t(1:4), #Perform the bias test
                  2, 
                  bias.test, 
                  AR.RE.4.FE)
names(AR.RE.4p) <- 1:4

VAR1.RE.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR1.RE.4.FE)
names(VAR1.RE.4p) <- 1:4

VAR2.RE.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR2.RE.4.FE)
names(VAR2.RE.4p) <- 1:4

VAR3.RE.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR3.RE.4.FE)
names(VAR3.RE.4p) <- 1:4

AR.RO.4p <- apply(t(1:4),
                  2, 
                  bias.test, 
                  AR.RO.4.FE)
names(AR.RO.4p) <- 1:4

VAR1.RO.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR1.RO.4.FE)
names(VAR1.RO.4p) <- 1:4

VAR2.RO.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR2.RO.4.FE)
names(VAR2.RO.4p) <- 1:4

VAR3.RO.4p <- apply(t(1:4),
                    2, 
                    bias.test, 
                    VAR3.RO.4.FE)
names(VAR3.RO.4p) <- 1:4

bias4<- rbind(AR.RE.4p,  # binding all test results in one matrix
              AR.RO.4p,
              VAR1.RE.4p,
              VAR1.RO.4p,
              VAR2.RE.4p,
              VAR2.RO.4p,
              VAR3.RE.4p,
              VAR3.RO.4p)

colnames(bias4) <- paste('Horizon', 
                         1:4, 
                         sep = ' ')


### Horizon 8

AR.RE.8p <- apply(t(1:8), #Perform the bias test
                  2, 
                  bias.test, 
                  AR.RE.8.FE)
names(AR.RE.8p) <- 1:8

VAR1.RE.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR1.RE.8.FE)
names(VAR1.RE.8p) <- 1:8

VAR2.RE.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR2.RE.8.FE)
names(VAR2.RE.8p) <- 1:8

VAR3.RE.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR3.RE.8.FE)
names(VAR3.RE.8p) <- 1:8

AR.RO.8p <- apply(t(1:8),
                  2, 
                  bias.test, 
                  AR.RO.8.FE)
names(AR.RO.8p) <- 1:8

VAR1.RO.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR1.RO.8.FE)
names(VAR1.RO.8p) <- 1:8

VAR2.RO.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR2.RO.8.FE)
names(VAR2.RO.8p) <- 1:8

VAR3.RO.8p <- apply(t(1:8),
                    2, 
                    bias.test, 
                    VAR3.RO.8.FE)
names(VAR3.RO.8p) <- 1:8

bias8 <- rbind(AR.RE.8p, # binding all test results in one matrix
              AR.RO.8p,
              VAR1.RE.8p,
              VAR1.RO.8p,
              VAR2.RE.8p,
              VAR2.RO.8p,
              VAR3.RE.8p,
              VAR3.RO.8p)


colnames(bias8) <- paste('Horizon', 
                         1:8, 
                         sep = ' ')


bias <- rbind(bias2[,'Horizon 2'], # This paper just compares the relevant horizons of the bias test.
              bias4[, 'Horizon 4'], 
              bias8[,'Horizon 8'])

rownames(bias) <- paste('Horizon', 
                        c(2, 4, 8), 
                        sep = ' ')

bias

#               AR.RE.2p AR.RO.2p VAR1.RE.2p VAR1.RO.2p VAR2.RE.2p VAR2.RO.2p VAR3.RE.2p VAR3.RO.2p
# Horizon 2     0.50     0.41       0.37       0.14       0.43       0.56       0.32       0.58
# Horizon 4     0.32     0.24       0.17       0.03       0.24       0.36       0.16       0.40
# Horizon 8     0.35     0.27       0.19       0.03       0.27       0.39       0.17       0.42             AR.RE.2p AR.RO.2p VAR1.RE.2p VAR1.RO.2p VAR2.RE.2p VAR2.RO.2p VAR3.RE.2p VAR3.RO.2p



# Conclusion: 
# Only the VAR1 model used rolling forecast technique to estimate the forecast is 
# biased in the horizon 4 and horizon 8. It makes sense from the plot, VAR 1 
# seems to be biased because it systematically underestimates the future
# inflation rate. The other models are unbiased. Thus, these models fullfil the 
# property of an optimal forecast that the expected forecast error is zero
# given that inflation rate probably has a symmetric cost function. Overshooting 
# and undershooting are equally painful. Considering a central bank that expect the
# inflation rate to decrease by a lot due to biased forecast and wrongly promotes a 
# dovish policy and thereby increasing unwillingly the inflation rate which might 
# already be sufficiently high.


############################################################################################
########### 4.3 Diebold-Mariano-Test #######################################################
############################################################################################
########### 4.3.2 Test statistic ###########################################################
############################################################################################
# In this section, we test for forecasting accuracy and thus by pairwaise comparing the 
# model and testing th DM test, which has a H0 of equal precision.
# H0: C[FE^a]) = C[FE^b], the loss difference of two forecasts are the same.

prognose2 <- c('AR.RE.2.FE', 
               'AR.RO.2.FE',
               'VAR1.RE.2.FE',
               'VAR1.RO.2.FE',
               'VAR2.RE.2.FE',
               'VAR2.RO.2.FE',
               'VAR3.RE.2.FE',
               'VAR3.RO.2.FE') 
testlist2 <- combn(prognose2, 2) # Creating the recombination of all models for pairwise comparison

prognose4 <- c('AR.RE.4.FE',
               'AR.RO.4.FE',
               'VAR1.RE.4.FE',
               'VAR1.RO.4.FE',
               'VAR2.RE.4.FE',
               'VAR2.RO.4.FE',
               'VAR3.RE.4.FE',
               'VAR3.RO.4.FE') 
testlist4 <- combn(prognose4, 2)

prognose8 <- c('AR.RE.8.FE',
               'AR.RO.8.FE',
               'VAR1.RE.8.FE',
               'VAR1.RO.8.FE',
               'VAR2.RE.8.FE',
               'VAR2.RO.8.FE',
               'VAR3.RE.8.FE',
               'VAR3.RO.8.FE') 
(testlist8 <- combn(prognose8, 2))


DM<- cbind(apply(testlist2, 2, DM.TEST, 2),
           apply(testlist4, 2, DM.TEST, 4),
           apply(testlist8, 2, DM.TEST, 8))

colnames(DM) <- paste('Horizon',
                      c(2,4,8),
                      sep =' ')

# Combination according to testlist to name the rows appropriately
rownames(DM) <- c('AR.RE vs. AR.RO', 'AR.RE vs. VAR1.RE', 'AR.RE vs. VAR1.RO',
                  'AR.RE vs. VAR2.RE', 'AR.RE vs. VAR2.RO', 'AR.RE vs. VAR3.RE',
                  'AR.RE vs. VAR3.RO', 'AR.RO vs. VAR1.RE', 'AR.RO vs. VAR1.RO',
                  'AR.RO vs. VAR2.RE', 'AR.RO vs. VAR2.RO', 'AR.RO vs. VAR3.RE',
                  'AR.RO vs. VAR3.RO', 'VAR1.RE vs. VAR1.RO', 'VAR1.RE vs. VAR2.RE',
                  'VAR1.RE vs. VAR2.RO', 'VAR1.RE vs. VAR3.RE','VAR1.RE vs. VAR3.RO',
                  'VAR1.RO vs. VAR2.RE', 'VAR1.RO vs. VAR2.RO','VAR1.RO vs. VAR3.RE',
                  'VAR1.RO vs. VAR3.RO','VAR2.RE vs. VAR2.RO', 'VAR2.RE vs. VAR3.RE',
                  'VAR2.RE vs. VAR3.RO','VAR2.RO vs. VAR3.RE', 'VAR2.RO vs. VAR3.RO',
                  'VAR3.RE vs. VAR3.RO')

DM # Visualising the DM test
# Note that VAR1.RO is biased.

#                     Horizon 2 Horizon 4 Horizon 8
# AR.RE vs. AR.RO          0.43      0.27      0.79  # No rejection of the null of equal loss difference
# AR.RE vs. VAR1.RE        0.06      0.00      0.00  # Rejection of the null of equal loss difference for h=4, 8 
# AR.RE vs. VAR1.RO        0.04      0.00      0.00  # Rejection of the null of equal loss difference
# AR.RE vs. VAR2.RE        0.56      0.51      0.57  # No rejection of the null of equal loss difference
# AR.RE vs. VAR2.RO        0.71      0.89      0.42  # No rejection of the null of equal loss difference
# AR.RE vs. VAR3.RE        0.80      0.80      0.60  # No rejection of the null of equal loss difference
# AR.RE vs. VAR3.RO        0.17      0.16      0.59  # No rejection of the null of equal loss difference
# AR.RO vs. VAR1.RE        0.03      0.00      0.01  # Rejection of the null of equal loss difference
# AR.RO vs. VAR1.RO        0.01      0.01      0.00  # Rejection of the null of equal loss difference for h=2,4,8
# AR.RO vs. VAR2.RE        0.34      0.82      0.77  # No rejection of the null of equal loss difference
# AR.RO vs. VAR2.RO        0.48      0.43      0.74  # No rejection of the null of equal loss difference
# AR.RO vs. VAR3.RE        0.60      0.61      0.22  # No rejection of the null of equal loss difference
# AR.RO vs. VAR3.RO        0.72      0.04      0.16  # Rejection of the null of equal loss difference for h=4
  
# VAR1.RE vs. VAR1.RO      0.49      0.83      0.31  # No rejection of the null of equal loss difference
# VAR1.RE vs. VAR2.RE      0.13      0.01      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RE vs. VAR2.RO      0.07      0.00      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RE vs. VAR3.RE      0.16      0.02      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RE vs. VAR3.RO      0.02      0.00      0.00  # Rejection of the null of equal loss difference for h=2,4,8
# VAR1.RO vs. VAR2.RE      0.09      0.00      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RO vs. VAR2.RO      0.06      0.00      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RO vs. VAR3.RE      0.12      0.00      0.00  # Rejection of the null of equal loss difference for h=4,8
# VAR1.RO vs. VAR3.RO      0.02      0.00      0.00  # Rejection of the null of equal loss difference for h=2,4,8
# VAR2.RE vs. VAR2.RO      0.57      0.55      0.92  # No rejection of the null of equal loss difference
# VAR2.RE vs. VAR3.RE      0.57      0.21      0.20  # No rejection of the null of equal loss difference
# VAR2.RE vs. VAR3.RO      0.18      0.29      0.43  # No rejection of the null of equal loss difference
# VAR2.RO vs. VAR3.RE      0.90      0.87      0.41  # No rejection of the null of equal loss difference
# VAR2.RO vs. VAR3.RO      0.01      0.09      0.07  # Rejection of the null of equal loss difference for h=2
# VAR3.RE vs. VAR3.RO      0.47      0.81      0.92  # No rejection of the null of equal loss difference


### Conclusion:
# ARIMA vs other models is the most interesting for this paper because not rejecting the
# H0 for those would infer that including more predictors did not yield in more accuracy

## ARIMA vs. VAR1 is rejected at nearly any horizon. Thus, the two models do not have the 
# same loss difference. Hence, the forecasts are not equally accurate. However, note 
# that VAR1.RO for horizon 4 and 8 is biased and might infer with the results. The forecasts plots suggests
# that ARIMA seems to be closer to the actual inflation rate and has probably the lower
# loss difference.


## ARIMA vs. VAR 2: in no instance is the H0 rejected independent of horizon and independent
# of forecasting technique. This highlights that ARIMA and the VAR2 model have the same accuracy.
# Thus, including two variables do not lead to a rejecting of the H0 and the ARIMA with one 
# predictors did an equally good job to forecast the inflation rate

## ARIMA vs. VAR 3: Only in one instance is the H0 rejected. AR.RO vs. VAR3.RO for horizon 4.
# Except that one instance, the H0 of equal forecast precision is not rejected meaning that
# performing a forecast with one predictor yielded in the same accuracy than the forecast
# with three additional predictors. No graphical inference is possible from the plots because 
# for the given horizon the two do not look different.


## VAR1 vs. VAR2
# The H0 of equal accuracy is independent of the forecasting technique not rejected for 
# the horizon 2. Thus, for two horizons forecast, the additional variable did not yield
# in additional precision. However, for horizons 4 and 8 is the H0 rejected suggesting 
# that the models did not equally good forecasts in longer horizons. The forecasts 
# plots suggest that the VAR2 model has a lesser distance to the actual infaltion rate 
# indicating that the VAR2 model might have the lower loss difference.
# Please not that VAR1.RO for horizon 4 and 8 is biased.


## VAR1 vs. VAR3
# The H0 is not rejected for horizon 2 comparing the recursive technique with the recursive
# technique and the rolling with the rolling technique. In any other instance is the H0 
# rejected. VAR1.RO is biase for horizon 4 and 8.
# The plots show that for horizon 2 

## VAR2 vs. VAR3
# The H0 is only rejected for VAR2.RO vs. VAR3.RO for horizon 2. In the other instances
# the H0 is not rejected. Hence, indicating that the increase in the additional variable 
# did predominantly yield in no higher accuracy.




############################################################################################
############################################################################################
############################################################################################
