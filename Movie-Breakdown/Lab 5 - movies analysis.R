install.packages("ggplot2movies")
data("movies", package = "ggplot2movies")

#exercise 5.44
#28819 movies in help files
#check the actual movies!!
nrow(movies)
#note the different values... 

spenny <- movies[order(movies$budget, decreasing = TRUE)[1],c("budget")]
spenny


# We can just run the max directly on the budget variable.
max(movies$budget)

#note: max(movies$budget) didn't work as the 
#budget doesn't appear for all the movie data 
#can use:
max(movies$budget,na.rm = TRUE)

movies.withbudget <- movies[!is.na(movies$budget),]
#whichmax only gives the first maximum budget (so only)
#spiderman appears as it is first
#movies.withbudget[which.max(movies.withbudget$budget),]

movies.withbudget[movies.withbudget$budget==200000000,]

rm(movies.withbudget)
?hist
#of length (in minutes) of all films
#using defaults for function plotting

hist(movies$length)

boxplot(movies$length,horizontal = TRUE)
movies[movies$length > 1000,]

rm(movies.withbudget)
?hist
#of length (in minutes) of all films
#using defaults for function plotting

hist(movies[movies$length <= 180, "length"], breaks = seq(0,180,1))

rm(movies.withbudget)
?hist
#of length (in minutes) of all films
#using defaults for function plotting

h <- hist(movies[movies$length <= 180, "length"], breaks = seq(0,180,1))
str(h)

# NOTE: check the help file ?order and you'll see order sorts ascending by
# default! So we need to change to descending.
desc.length <- order(h$counts, decreasing = TRUE)
h$breaks[desc.length[1:3]]

plot(movies$year, movies$length, ylim = c(0,500))
abline(lm(length ~ year, movies),
       col = "red")
lines(lowess(movies$year, movies$length),
      col = "green")


bootstrap <- function(x,ci = 0.99, B= 10000){
  S.star <- rep(0,B)
  for(b in 1:B) {
    x.star <- sample(x, replace = TRUE)
    S.star[b] <- median(x.star)
  }
  lower <- sort(S.star)[round((1-ci)/2*B)]
  upper <- sort(S.star)[round((ci+(1-ci)/2)*B)]
  return(c(lower,upper))
  
}
# part (i) above

install.packages("ggplot2movies")
data("movies", package = "ggplot2movies")


bootstrap <- function(x,ci = 0.99, B= 10000){
  S.star <- rep(0,B)
  for(b in 1:B) {
    x.star <- sample(x, replace = TRUE)
    S.star[b] <- median(x.star)
  }
  lower <- sort(S.star)[round((1-ci)/2*B)]
  upper <- sort(S.star)[round((ci+(1-ci)/2)*B)]
  return(c(lower,upper))
  
}
# part (i) above

range(movies$year)

med.CIs <- NULL
for(decade in seq(1890, 2000, 10)) {
  new.CI <- bootstrap(movies[movies$year >= decade & movies$year < decade+10,"length"])
  med.CIs <- rbind(med.CIs,
                   data.frame(decade = decade,
                              lower = new.CI[1],
                              upper = new.CI[2]))
}
med.CIs

install.packages("ggplot2movies")
data("movies", package = "ggplot2movies")


bootstrap <- function(x,ci = 0.99, B= 10000){
  S.star <- rep(0,B)
  for(b in 1:B) {
    x.star <- sample(x, replace = TRUE)
    S.star[b] <- median(x.star)
  }
  lower <- sort(S.star)[round((1-ci)/2*B)]
  upper <- sort(S.star)[round((ci+(1-ci)/2)*B)]
  return(c(lower,upper))
  
}
# part (i) above

range(movies$year)

med.CIs <- NULL
for(decade in seq(1890, 2000, 10)) {
  new.CI <- bootstrap(movies[movies$year >= decade & movies$year < decade+10,"length"])
  med.CIs <- rbind(med.CIs,
                   data.frame(decade = decade,
                              lower = new.CI[1],
                              upper = new.CI[2]))
}
med.CIs
# As explained in the lecture, we need to use plot first as this creates a new
# plot, so we will call with line type.
# Also we add 5 to centre the point in the middle of the decade
plot(med.CIs$decade+5, med.CIs$lower, type="l")
# Then we add the upper confidence interval by using lines to add to the plot
lines(med.CIs$decade+5, med.CIs$upper, type="l")