x <- data.frame(sugars = c(1.9,2.6,2.3,1.9), 
                pH = c(3.51,3.20,3.26,3.16), 
                alcohol = c(9.4,9.8,9.8,9.8))

y <- data.frame(sulphates = wq.red[,10],
                quality = wq.red[12])
plot(y)

z <- wq.red[,10]
k <- median(z)

#the below separates rows where sulphates are above and below the median
low.s <- wq.red[wq.red$sulphates <= median(wq.red$sulphates),]
high.s <- wq.red[wq.red$sulphates > median(wq.red$sulphates),]

#ecdf of sulphur quality 
low.ecdf <- ecdf(low.s$quality)
plot(low.ecdf)
high.ecdf <- ecdf(high.s$quality)
lines(high.ecdf,col = "red")

str(low.ecdf)
str(high.ecdf)

#to compute probabilities, remember ecdf automatically does less than or equal to the certain value
#so P(x is less than or equal to 5) is
low.ecdf(5)
#and then P(x is greater than or equal to 7)is 
1-high.ecdf(6)

?t.test()
t.test(low.s$quality)
t.test(high.s$quality)

t.test