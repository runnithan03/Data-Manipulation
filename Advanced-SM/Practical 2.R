# Load the data
Football <- read.csv("FootballData.csv")
View(Football)

##### (a) #####
# Helpfully, R splits factors into columns, so we can use it to create an
# indicator variable for each team both when home and away

# Creating model matrices
#the model.matrix() function creates binary (0/1) indicator variables (dummy variables) for each level of the factor, 
#except for the first level which is omitted to avoid multicollinearity (this is known as the "reference level").

X_home <- model.matrix( ~ 0+HomeTeam, data=Football)
X_away <- model.matrix( ~ 0+AwayTeam, data=Football)


#Each column in X_home corresponds to a different home team, and each row corresponds to a match. 
#A "1" in a column indicates that the team associated with that column was playing at home in that match, and "0" otherwise.
#The - 1 in the formula removes the intercept term, ensuring that a column is created for each team without omitting one for the intercept.

#A similar idea is used for X_away but we have a disadvantage here so by multiplying by -1, the model can interpret the away team's 
#participation as a negative influence on the probability of a home win, consistent with the model's assumption.

#The final step involves combining these two matrices to form a single matrix X that represents the home and away status for all teams 
#across all matches. This is done by adding the X_home and X_away matrices together. In the resulting matrix, for any given match (row), 
#there will be a "1" in the column corresponding to the home team, a "-1" in the column corresponding to the away team, and "0"s elsewhere.

# Combining matrices and converting to a data frame
X <- data.frame(X_home + X_away)

#To avoid a singular (non-invertible) matrix, which cannot be used in regression analysis due to perfect multicollinearity, 
#one team's column is typically dropped. This team serves as the reference category against which the strengths of all other teams are compared. 
#In the provided example, the first team alphabetically (e.g., Arsenal) is dropped.
# Dropping the first team (e.g., Arsenal)

View(X_home)
View(X_away)

# X_home has 1 when home, and X_away has 1 when away, so we want the difference
# Convert the difference into a data frame so variables (columns) are named
X <- as.data.frame(X_home - X_away)

nrow(X)
ncol(X)

# Optionally, we can pretty up the names by removing the "HomeTeam" part 
# at the start
names(X) <- substring(names(X), 9)
# Use the "stringr" package to replace the empty space in team names with a dot.
library(stringr)
names(X) <- str_replace_all(names(X), c(" " = ".", "," = ""))

# To avoid a singular solution due to collinearity of indicators, we simply
# drop the first team, Arsenal
X <- X[,-1]
dim(X)

##### (b) #####
# Finally, create vector y of home team wins and make data frame
y <- ifelse(Football$FTR == "H", 1, 0)

# Alternative ways to get the response vector
# y <- as.integer(Football$FTR == "H")
# y <- 1*(Football$FTR == "H")

matchdata <- cbind(y = y, X)

View(matchdata)


##### (c) #####
# Fit Bradley-Terry model via logistic regression
fit1 <- glm(y ~ 0+., matchdata, family = binomial(link=logit))
#the use of 0+ means you don't include the intercept in the model
summary(fit1)

# Look at the coefficients sorted (remember that Arsenal is the reference team
# so its estimate is 0, and any teams with a estimate greater than 0 are stronger
# than Arsenal and any teams with a estimate less than 0 are weaker than Arsenal
# in terms of a high home-win probability)
sort(coef(fit1), decreasing = TRUE)

# Tottenham seems to be stronger under the model than Man City.
# Notice that looking at the away games for Tottenham and Man City,
# they achieved several draws which count as a 'victory' when playing away
# with this model. Hence, probability of away victories likely inflated.
Football[Football$AwayTeam=="Tottenham", 2:7]
Football[Football$AwayTeam=="Man City", 2:7]


##### (d) #####
# We expect a logistic regression to have dispersion parameter 1
# Estimating it from the model fit, we get
fit1$deviance / fit1$df.residual
# which looks not too far from 1, close to 1 suggests no overdispersion


##### (e) #####
# A home team advantage could be modelled as a constant multiplicative effect
# on the odds of a home team win. This is exactly the effect an intercept term
# would have.
# So, if we drop the 0+ to allow R to put in an intercept:
fit2 <- glm(y ~ ., matchdata, family = binomial(link=logit))
summary(fit2)
# Intercept is actually negative! However, in fact there's quite a lot of
# uncertainty here, so it's just saying there isn't enough data for us to
# distinguish it
confint(fit2, "(Intercept)")



##### (f) #####
# First, get the full Fisher information matrix
Finv_betahat <- summary(fit2)$cov.scaled
# For this confidence region, we just want the submatrix involving 
# Man City and Liverpool.
# Also, we want the inverse for the Mahalanobis distance
F_lmc <- solve(Finv_betahat[c(11, 13), c(11, 13)])
F_lmc

# Get the MLEs for these two teams
betahat_lmc <- coef(fit2)[c(11, 13)]

# Let's setup a grid of strengths which we'll check the Mahalanobis distance
# against the chi-squared critical value
liverpool <- seq(-4, 4, length.out = 300)
man_city <- seq(-4, 4, length.out = 300)

# Use the outer function to evaluate over a grid
HessCR <- outer(liverpool, man_city, 
                Vectorize(function(beta_l, beta_mc) {
                  beta <- c(beta_l, beta_mc)
                  t(betahat_lmc - beta) %*% F_lmc %*% (betahat_lmc - beta)
                }
                ))

# The image function now lets us colour the part we're interested in
image(liverpool, man_city, HessCR > qchisq(0.95, 2))
# and mark the location of the MLE for reference
points(betahat_lmc[1], betahat_lmc[2], pch = 3, col = "black")

# You could also add the individual confidence intervals calculated by R as
# horizontal and vertical lines, but note it is doing something called profiling
# to get the confidence interval accounting for all other variables so you will
# notice a small discrepancy between the marginal and the extremes of the joint
abline(h = confint(fit2, "Man.City"))
abline(v = confint(fit2, "Liverpool"))



##### (g) #####
# We need to compute the model with the three teams fixed to have the same
# coefficient. To do this, simply tell R to remove the individual predictors
# and insert a new predictor formed from all three
fit3 <- glm(y ~ . - Everton - Burnley - Sheffield.United
            + I(Everton + Burnley + Sheffield.United),
            matchdata, family = binomial(link=logit))
summary(fit3)

# In logistic regression, dispersion is 1, so the likelihood ratio test
# statistic can be found by just taking the difference of the deviances
fit3$deviance - fit2$deviance
# Difference in degrees of freedom (should be 2 as we've replaced 3 coefficients
# with 1)
fit3$df.residual - fit2$df.residual
# Test is against chi-sq(v=2)
qchisq(0.95, 2)
# LR is not larger, therefore not enough evidence to reject ... plausibly all
# teams up for relegation are equally weak



##### (h) #####
# We need to predict a new match.  Easiest to pull a row from X, zero out and
# then set the home/away teams we need
new_match <- X[1,]
new_match[,1:19] <- 0
new_match$Man.City <- +1
new_match$Chelsea <- -1

View(new_match)

# Probability of Chelsea win is probability of away win, 
# so 1 minus predicted probability of home win
1-predict(fit3, new_match, type = "response")



##### (i) #####
# Split the data up
training <- matchdata[1:150,]
testing <- matchdata[-(1:150),]

# Fit the model on the training data only
fit4 <- glm(y ~ ., training, family = binomial(link=logit))
summary(fit4)

# Predict on the testing
pred <- predict(fit4, testing, type = "response")

# Now produce a table, where we want to compare the truth to our prediction
res <- table(truth = testing$y,
             prediction = ifelse(pred > 0.5, 1, 0))
res

# Hence overall accuracy (%) is
(res[1,1]+res[2,2])/sum(res)*100

# Note: not very much high prediction accuracy here, perhaps we should be a bit
# more careful in choosing the cut-off point (here 0.5) for prediction of win
# for the unseen data. Or we may use more data to train the model and test it.

##### (j) #####
# To use a rescaled Bionamial model, one can consider the number of goals
# scored in a match. For this, one can define the response variable to be
# the number of goals by home team (in column FTHG) while "n" would be the
# total goals scored in the match (i.e., FTHG+FTAG).


#---------------------------------------------------------------------------------
# One final thing: one might be interested to predict all the remaining matches
# considering the remaining home and away matches for all teams, and predict the
# winners to approximately calculate the total points for all teams at the end
# of this season. A more accurate way it to fit a "multinomial logistic" with
# three categories "win,lose,draw", which could be more interesting but this is
# beyond the level of our course.
#---------------------------------------------------------------------------------
