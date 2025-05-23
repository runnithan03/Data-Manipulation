# Example: US Polio Data

library( "gamlss.data" )
data( "polio" )
uspolio <- as.data.frame( matrix( c( 1:168, t( polio ) ), ncol = 2 ) )
colnames( uspolio ) <- c("time", "cases")

plot( uspolio, type = "h" )


# Fit Poisson model with linear time trend
polio.glm <- glm( cases ~ time, family = poisson( link = log ), data = uspolio )

# Look at the summary
summary( polio.glm )

# Plot the model
plot(1970 + ((uspolio$time - 1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time - 1)/12), polio.glm$fitted)


# Fit Poisson model with linear trend and seasonal (annual) component
polio1.glm <- glm(cases~time + I(cos(2*pi*time/12)) + I(sin(2*pi*time/12)), 
                  family=poisson(link=log), data=uspolio)

summary(polio1.glm)

plot(1970 + ((uspolio$time - 1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time - 1)/12), polio1.glm$fitted, col=2)



# Fit Poisson model with linear trend and seasonal (annual + sixmonthly) component
polio2.glm <- glm(cases~time + I(cos(2*pi*time/12)) + I(sin(2*pi*time/12)) 
                  + I(cos(2*pi*time/6)) + I(sin(2*pi*time/6)),
                  family=poisson(link=log), data=uspolio)

summary(polio2.glm)

plot(1970 + ((uspolio$time - 1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time - 1)/12), polio2.glm$fitted, col=3)


# Add in temperature data.
# average annual temperature data over the 14 years.
temp_data <- rep(c(5.195, 5.138, 5.316, 5.242, 5.094, 5.108, 5.260, 5.153, 
                   5.155, 5.231, 5.234, 5.142, 5.173, 5.167), each = 12 )

# Scale the data so that it plots nicely.
scaled_temp = 10 * (temp_data - min(temp_data))/(max(temp_data) - min(temp_data))
uspolio$temp = scaled_temp

# Plot temperature data against cases data to see interest.
plot(1970 + ((uspolio$time - 1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time - 1)/12), uspolio$temp, col="red")

# Construct GLM.
polio3.glm <- glm(cases~time + temp + I(cos(2*pi*time/12)) + I(sin(2*pi*time/12)) 
                  + I(cos(2*pi*time/6)) + I(sin(2*pi*time/6)),
                  family=poisson(link=log), data=uspolio)

summary(polio3.glm)

plot(1970 + ((uspolio$time - 1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time - 1)/12), polio3.glm$fitted, col="red")


###########################################

# Example: Estimation of $\phi$ and prediction on Hospital Stay Data

data(hosp, package="npmlreg")

hosp.glm <- glm(duration~age+temp1, data=hosp, family=Gamma(link=log))

summary(hosp.glm)

# Compute dispersion by hand
1/(hosp.glm$df.res)*sum( (hosp$duration-hosp.glm$fitted)^2/(hosp.glm$fitted^2) )

# Predict duration of stay for a new individual with age 60 and temperature 99
predict(hosp.glm, newdata=data.frame(age=60, temp1=99), type="response")

# Example: Hospital Stay Data

# Load data
data(hosp, package="npmlreg")

# Fit a GLM with Gamma family to the data
hosp.glm <- glm(duration~age+temp1, data=hosp, family=Gamma(link=log))
summary(hosp.glm)


# Predict duration of stay for a new individual with age 60 and temperature 99

# First attempt. This would not work.
predict(hosp.glm, newdata=data.frame(age=60, temp1=99))

# This works.
exp(predict(hosp.glm, newdata=data.frame(age=60, temp1=99)))

# Another way to make prediction
predict(hosp.glm, newdata=data.frame(age=60, temp1=99), type="response")


# Compute 95% confidence interval for the expected mean of the new individual

# First attempt. This would not work.
predict(hosp.glm, newdata=data.frame(age=60, temp1=99), type="response", interval="confidence")

# Do it manually
# Compute the predicted linear predictor as above
lphat  <- predict(hosp.glm, newdata=data.frame(age=60, temp1=99))

# Extract the covariance
varhat <- summary(hosp.glm)$cov.scaled # = F^(-1)(betahat)

# Define new data point
x0 = c(1, 60, 99)

# Compute the width of the interval for the linear predictor
span <- qnorm(0.975) * sqrt( x0 %*% varhat %*% x0 )

# Compute the interval for the mean
c(exp(lphat-span), exp(lphat+span))

# Note that this is quite large, as the dataset is small!


########################

# Covariance of the estimated parameters
( varhat <- summary(hosp.glm)$cov.scaled )

# Original model
summary(hosp.glm)

# Assume the dispersion is known
summary(hosp.glm, dispersion = 0.2690233)

##### Separation issue #####

# Separation on 1-d data

y  <- c(0, 0, 0, 0, 1, 1,  1,  1)
x1 <- c(1, 2, 3, 4, 5, 6, 10, 11)
plot(y~x1, col=c("red","blue")[factor(y)])

m1 <- glm(y~x1, family=binomial)
summary(m1)


# Separation on 2-d data

y  <- c(0, 0,  0,  0, 1,  1,  1,  1)
x1 <- c(1, 2,  3,  4, 5,  6, 10, 11)
x2 <- c(5, 4, -1, -1, 1, -1,  2,  3)
plot(x2~x1, col=c("red","blue")[factor(y)])

m2 <- glm(y~x1+x2, family=binomial)
summary(m2)


# Check infinite estimates using the 'detectseparation' package on real data
install.packages("detectseparation")
# Get the endometrial cancer data and fit a logistic regression model
library("detectseparation")
data("endometrial", package = "detectseparation")
endo_glm <- glm(HG ~ NV + PI + EH, family = binomial(), data = endometrial)
summary(endo_glm)

# Check and plot infinite estimates
(inf_check <- check_infinite_estimates(endo_glm))
plot(inf_check)

# Using detect_separation method
endo_sep <- glm(HG ~ NV + PI + EH, data = endometrial,
                family = binomial("logit"),
                method = "detect_separation")
endo_sep


##############################
##### Hauck–Donner Effect ####

n <- 50  # Half the number of data points
x <- seq(from = 0, to = 1, length.out = 2*n)
y <- rep(0, length.out = 2*n)

t <- 30  # Number of data points to flip response from 0 to 1
y1 <- y
y1[(n+1):(n+t)] <- 1
plot(y1~x, col=c("red","blue")[factor(y1)])

model <- glm(y1~x, family = binomial)
summary(model)


# Compute and plot z values (Wald statistics) for all t from 1 to 49

w = c()
t_range <- 1:49

for (t in t_range) {
  y1 <- y
  y1[(n+1):(n+t)] <- 1
  
  model <- glm(y1~x, family = binomial)
  summary(model)
  
  wt = coef(summary(model))[2, "z value"] 
  w <- append(w, wt)
}

w[30] # Just to double check

plot(w~t_range, type="l")

# Test if a model is a good fit

# Example: US polio data

# Load and plot the data
library("gamlss.data")
data("polio")
uspolio <- as.data.frame(matrix(c(1:168, t(polio)), ncol=2))
colnames(uspolio) <- c("time", "cases")

plot(uspolio, type="h")

# Create Poisson GLM which includes time, six-month cycles 
# and twelve-month cycles (see Estimation chapter)
polio2.glm<- glm(cases~time + I(cos(2*pi*time/12))+I(sin(2*pi*time/12))
                 + I(cos(2*pi*time/6)) + I(sin(2*pi*time/6)), 
                 family=poisson(link=log), data=uspolio)
summary(polio2.glm)

plot(1970 + ((uspolio$time-1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time-1)/12), polio2.glm$fitted,col=4)

# Deviance
polio2.glm$dev

# Pearson statistic
sum((uspolio$cases-polio2.glm$fitted)^2 / polio2.glm$fitted)

# Either way, critical value at 5% level
qchisq(0.95, 162)

# So whether using deviance or Pearson statistic, we reject the null hypothesis of a good fit


# What happens if we add in the temperature data?

# Read in temperature data and scale it
temp_data <- rep(c(5.195, 5.138, 5.316, 5.242, 5.094, 5.108, 5.260, 5.153, 
                   5.155, 5.231, 5.234, 5.142, 5.173, 5.167), each = 12 )
scaled_temp = 10 * (temp_data - min(temp_data))/(max(temp_data) - min(temp_data))
uspolio$temp = scaled_temp

# Construct GLM
polio3.glm <- glm(cases~time + temp 
                  + I(cos(2*pi*time/12)) + I(sin(2*pi*time/12))
                  + I(cos(2*pi*time/6)) + I(sin(2*pi*time/6)), 
                  family=poisson(link=log),data=uspolio)
summary(polio3.glm)

plot(1970 + ((uspolio$time-1)/12), uspolio$cases, type="h")
lines(1970 + ((uspolio$time-1)/12), polio3.glm$fitted, col="red")

# Deviance
polio3.glm$dev

# Pearson statistic
sum((uspolio$cases-polio3.glm$fitted)^2 / polio3.glm$fitted)

# Now, critical value at 5% level
qchisq(0.95, 161)

# Still reject, though note deviance and Pearson statistic both reduced



###################################

# Inspecting residual plots

# Example: US polio data

par(mfrow=c(2,1))
plot(uspolio$time, uspolio$cases, type="h")
lines(uspolio$time, polio2.glm$fitted, col="blue")
plot(uspolio$time, residuals(polio2.glm, type="deviance"), type="b")
abline(a=0,b=0)

cor(residuals(polio2.glm, type="deviance")[2:168], residuals(polio2.glm, type="deviance")[1:167])

# Small but noticeable autocorrelation

par(mfrow=c(2,1))
plot(uspolio$time, uspolio$cases, type="h")
lines(uspolio$time, polio3.glm$fitted, col="red")
plot(uspolio$time, residuals(polio3.glm, type="deviance"), type="b")
abline(a=0,b=0)

cor(residuals(polio3.glm, type="deviance")[2:168], residuals(polio3.glm, type="deviance")[1:167])

# Autocorrelation now smaller. Small enough?


# Example: Hospital stay data

data(hosp, package="npmlreg")
hosp.glm <- glm(duration~age+temp1, data=hosp, family=Gamma(link=log))

par(mfrow=c(2,2))
plot(residuals(hosp.glm, type="deviance"))
plot(hosp.glm$fitted, residuals(hosp.glm, type="pearson"))
plot(hosp$age, residuals(hosp.glm, type="deviance"))
plot(hosp$temp1, residuals(hosp.glm, type="deviance"))

# Look OK?

# Autocorrelations
cor(residuals(hosp.glm, type="deviance")[1:24], residuals(hosp.glm, type="deviance")[2:25])
cor(residuals(hosp.glm, type="pearson")[1:24], residuals(hosp.glm, type="pearson")[2:25])

# Little high


###################################

# Analysis of deviance

# Example: Hospital stay data

data(hosp, package="npmlreg")

# Full model
fit1 <- glm(duration~age+temp1+wbc1+antib+bact+serv, data=hosp,
            family=Gamma(link=log))
summary(fit1)

# Get the dispersion
summary(fit1)$dispersion

# Get the deviance table
anova(fit1)

# p value for test problem 1
1-pchisq((8.1722-5.1200)/0.2661922, 6)

# p value for test problem 2
1-pchisq(1.00299/0.2662922, 1)

# p value for test problem 3
1-pchisq((5.7849-5.12)/0.2662922, 4)

# Add chi square p-values
anova(fit1, test="Chisq")


# Fit the model with a different order of variables
fit2 <- glm(duration~wbc1+antib+bact+serv+temp1+age, data=hosp,
            family=Gamma(link=log))

summary(fit2)
summary(fit2)$dispersion
anova(fit2, test="Chisq")
