install.packages("npmlreg")
library(npmlreg)

data(irlsuicide)
irlsuicide[1,]
dim(irlsuicide)

#(a)
par(mfrow =c(1,3))
boxplot(irlsuicide$age, xlab = "age", ylab = "frequency")
boxplot(irlsuicide$sex, xlab = "sex", ylab = "frequency")
boxplot(irlsuicide$Region, xlab = "Region", ylab = "frequency")

#(b)
#(i)
?glm
irlsuicide1 <- irlsuicide
irlsuicide1$suicide_rate <- irlsuicide1$death/irlsuicide1$pop #define suicide rate as death/pop
irlsuicide1 <- irlsuicide1[,-c(4,7,8)] #removes smr, expected, death and pop
dimnames(irlsuicide1)

#scaling
#irlsuicide1$suicide_rate<- scale(irlsuicide1$suicide_rate)

model <- glm(suicide_rate~.-pop, family = binomial( link = log ), data = irlsuicide1, weights = pop) #for binomial case, weights are the number of trials
summ <- summary(model)
summ
#ID2,...,ID13 are NA because they are the same as the Region
#All Regions are significant except RegionLim. and RegionSEHB

#(ii)
irlsuicide1$ID <- factor(irlsuicide1$ID)
model_factor <- glm(suicide_rate~sex+age+ID, family = binomial( link = log ),data = irlsuicide1, weights = pop)
sumf <- summary(model_factor)
sumf
#ID2,...,ID13 are NA because they are the same as the Region
#All Regions are significant except RegionLim. and RegionSEHB

#(iii)
newdata = data.frame(ID = "3", sex = "1", age = "3")
predict(model_factor, newdata = newdata, type = "response")
#0.002239629

#(c)
#(i)
irlsuicide2 <- irlsuicide
dimnames(irlsuicide2)
model2 <- glm(death~ID+sex+age, family = poisson(link = log), data = irlsuicide2, offset = log(pop))


#(ii)
#test:
#irlsuicide2[,"3",,,"1","3",,], trying to slice this to get pop = 4989, why does it not work?
newdata = data.frame(ID = "3", sex = "1", age = "3", pop = 4989) #using the pop which corresponds to ID = "3", sex = "1", age = "3", pop = 4989
death <- predict(model2, newdata = newdata, type = "response")
death/4989
#0.00223995
#similar to b.iii

#(iii)
#how to estimate the dispersion present in the data?
#lecture notes: there is a formula for the dispersion

# Compute dispersion by hand
1/(model_factor$df.res)*sum( (irlsuicide1$suicide_rate-model_factor$fitted)^2/(model_factor$fitted^2) )
#[1] 0.1432886
1/(model2$df.res)*sum( (irlsuicide2$death-model2$fitted)^2/(model2$fitted^2) )
#[1] 0.1432757
#Very similar dispersions but binomial model has a slightly higher dispersion so slightly lower accuracy
#Makes sense because the data is continuous, so a binomial model would be harder to fit
#as it can only take discrete values

#(d)
b <- irlsuicide$death/irlsuicide$pop
c <- predict(model_factor, newdata = irlsuicide, type = "response")
d <- predict(model2, newdata = irlsuicide, type = "response")/irlsuicide$pop
matrix_obs <- cbind(b,c,d)
matrix_obs
plot(matrix_obs)

c-b #indicates the fitted values are very close to the raw observed rates

#the fitted suicide rates using models (b) and (c) are 
#very similar

#(e)
#don't think it is a generally sound approach to include a large 
#number of regional indicators into models of this type because 
#we don't have many observations












