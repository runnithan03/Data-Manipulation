install.packages("ggplot2")
data("diamonds",package = "ggplot2")
library("ggplot2")

# 53940 observations
# 10 variables 
#within the data 

#the above are correct but use R...
nrow(diamonds) #aka observations
ncol(diamonds) #aka variables 

#factor type variables...
#number of columns?
str(diamonds)

#says there are 3 ordered factors! 

#to get the first few observations... 
head(diamonds)
#automatically does first 6 rows... 

summary(diamonds)

sum(diamonds[diamonds$cut == "Ideal" & diamonds$color == "D" & diamonds$depth <= 60,
             "price"])

#so yes, less than $200,000

diamonds$ppc <- diamonds$price/diamonds$carat 
median(diamonds$ppc)

#carat x axis, ppc y axis 
diamonds$ppc <- diamonds$price/diamonds$carat
plot(diamonds$carat,diamonds$ppc)

#carat x axis, ppc y axis 
diamonds$ppc <- diamonds$price/diamonds$carat
plot(diamonds$carat,diamonds$ppc)
abline (h=mean(diamonds$ppc)) #h gives horizontal line
abline(v=1:5) #v gives vertical line 

high.ppc <- diamonds[diamonds$ppc>10000 & diamonds$carat >= 1 & diamonds$carat < 2,]
normal.ppc <- diamonds[diamonds$ppc<=10000 & diamonds$carat >= 1 & diamonds$carat < 2,]

table(high.ppc$clarity)

pricy <- diamonds[order(diamonds$price,decreasing = TRUE)[1], c('price','carat','cut','clarity')]
pricy

#order(diamonds$price) provides the row numbers from lowest to highest price
#order(diamonds$price, decreasing = TRUE) switches this to highest to lowest price
#order(diamonds$price, decreasing = TRUE)[1] gets just the first element, in other words the row number of the highest priced diamond
#The above is passed as the row value in the [ , ] and the vector of column names is provided in the column values.

pricy <- diamonds[order(diamonds$price,decreasing = TRUE)[1], c('price','carat','cut','clarity')]
improve <- diamonds[diamonds$carat > pricy$carat & diamonds$cut > pricy$cut & diamonds$clarity > pricy$clarity,]
improve

#saving...
pricy$price - min(improve$price)