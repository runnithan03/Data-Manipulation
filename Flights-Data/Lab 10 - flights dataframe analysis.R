install.packages("tidyverse")
library("tidyverse")

install.packages("nycflights13")
#library("nycflights13")

#make sure to include what the package is called... 
data("flights", package = "nycflights13")

str(flights)

#my solution:
not_2013 <- filter(flights, year != 2013)
nrow(not_2013)

#their solution:
flights |> 
  filter(year != 2013) |> 
  nrow()

#my solution:
my_birthday <- filter(flights,flights$month == 3, flights$day == 4)
nrow(my_birthday)

#their solution:
flights |> 
  filter(month == 3 & day == 4) |> 
  nrow()

#my solution with their solution correction:
num.flights <- flights |>
  group_by(year,month,day) |> #need year month day here 
  summarise(count = n()) |> #n needs no argument inside! 
  arrange(desc(count))
#don't need to separate out year,month,day and count 
num.flights

#filter() focus on a subset of rows
#arrange() reorder the rows
#select() focus on a subset of variables (columns)
#mutate() create new derived variables
#summarise() create summary statistics (collapsing many rows) by groupings

flights_late <- flights |>
  filter(arr_delay >= 120)

#their solution: rtq because it says by carrier...
#so think use group and summarise! 

flights |>
  group_by(carrier) |>
  summarise(avg_dep_delay = mean(dep_delay,na.rm = TRUE), 
            sd_dep_delay = sd(dep_delay,na.rm = TRUE),
            count = n())
#so sd left, mean middle and count right... 
#Carrier code US seems to be the best, with the lowest average 
#departure delay and lowest standard deviation. We would clearly 
#prefer them, because they achieve this on a large volume of flights, 
#not just by having a small number of flights that were on time.

data("airlines",package = "nycflights13")

#was never gonna get this bit below... 
airlines |> filter(carrier == "US") |> pull(name)

#My preferred airline is US Airways Inc.

hist(flights$dep_delay)
boxplot(flights$dep_delay)

#was not getting this...
hist(flights |> filter(dep_delay <= 120) |> pull(dep_delay))
boxplot(flights |> filter(dep_delay <= 120) |> pull(dep_delay))

#preamble for 5.61
#ifelse() takes 3 arguments:
#first should be a vector of logical conditions
#second should be any vector of the same length as the first argument
#third should also be any vector of the same length as the first argument
#look at that lab for the explanation....

#example to check:
ifelse(c(TRUE, TRUE, FALSE, TRUE, FALSE),
       c(1, 2, 3, 4, 5),
       c(1000, 2000, 3000, 4000, 5000))

x <- rexp(10) #gives 10 random values from the exponential distribution
#10 corresponds to the lambda in the exponential distribution...
#i.e. f(x) = lambda*e^(-lambda*x)

x <- ifelse(x > 1,
            rep(1, 10),
            x)
x
#the above replaces any values over 1 generated in x with 1
#think about it like we generate the TRUE/FALSE values via x > 1
#then we sub in 1 or x accordingly... 

#use mutate as the variable is changing! 
flights |>
  mutate(dep_delay = ifelse(dep_delay < 0,0,dep_delay)) |>
  group_by(carrier) |>
  summarise(avg_dep_delay = mean(dep_delay,na.rm = TRUE), 
            sd_dep_delay = sd(dep_delay,na.rm = TRUE),
            count = n())

data("airlines",package = "nycflights13")

#was never gonna get this bit below... 
airlines |> filter(carrier == "US") |> pull(name)
#therefore no change in conclusion... 
# although the average delay is now doubled when we don’t allow 
#early departures to compensate for late ones.

#remember floor, ceiling and round only apply to decimals! 
clock_to_minutes <- function(time) {
  hour_minutes <- floor(time/100)*60 
  minutes <- time - 100*floor(time/100)
  return (minutes + hour_minutes)
}  
clock_to_minutes(1259)

clock_to_minutes <- function(time) {
  hour_minutes <- floor(time/100)*60 
  minutes <- time - 100*floor(time/100)
  return (minutes + hour_minutes)
}  

flights2 <- flights |>
  mutate(dep_time = clock_to_minutes(dep_time), 
         arr_time = clock_to_minutes(arr_time))

#just to see the first few rows...
head(flights2)

#remember about is.na ig... 
# Total flights cancelled
flights |>
  filter(is.na(dep_time)) |> 
  nrow()

#Create a table showing the total number of cancellations by the hour of 
#scheduled departure, sorted in descending order by the count.

#remember to mutate when changing a column.... 
#and add "|>" when doing this type of data manipulation 
#also use arrange when sorting from highest to lowest! 
flights <- flights |>
  filter(is.na(dep_time)) |>
  mutate(sched_dep_time = floor(sched_dep_time/100)) |>
  group_by(sched_dep_time) |>
  summarise(cancelled = n()) |>
  arrange(desc(cancelled))

#now do proportion of flights cancelled each hour of the day 
#remember we do not need to remove flights which have a dep_time 
flights <- flights |>
  mutate(sched_dep_time = floor(sched_dep_time/100)) |>
  group_by(sched_dep_time) |>
  #note use of is.na function...
  summarise(total = n(), cancelled = sum(is.na(dep_time))) |>
  #and mutate used to change a column by a specific factor 
  #or making a new column! 
  mutate(prop = cancelled/total) |>
  arrange(desc(prop))

library("tidyverse")
library("lubridate")

data("flights", package = "nycflights13")
?nycflights13::flights

#part (i)
flights <- flights |>
  mutate(sched_dep_hour = floor(sched_dep_time/100), 
         sched_dep_min = sched_dep_time - 100*sched_dep_hour)

#part (ii)
flights <- flights |>
  mutate(sched_dep = make_datetime(year,month,day,sched_dep_hour,sched_dep_min, tz = "America/New_York")) |>
  select(-sched_dep_hour,-sched_dep_min)
#note: the minus sign removes these columns! 

#part (i)
flights |> 
  mutate(day_of_week = wday(sched_dep, label = TRUE)) |> 
  #label = TRUE gives us the day of the week as an ordered factor of character strings
  #such as "Sunday"
  group_by(day_of_week) |> 
  #remember na.rm means we don't consider values = N/A
  summarise(avg_delay = mean(dep_delay, na.rm = TRUE))

#part (ii)
flights |> 
  #self-explanatory really! 
  filter(dep_delay >= 0) |>
  mutate(day_of_week = wday(sched_dep, label = TRUE)) |> 
  #label = TRUE gives us the day of the week a san ordered factor of character strings
  #such as "Sunday"
  group_by(day_of_week) |> 
  #remember na.rm means we don't consider values = N/A
  summarise(avg_delay = mean(dep_delay, na.rm = TRUE))

#prior exercise
flights |> 
  mutate(day_of_week = wday(sched_dep, label = TRUE)) |> 
  #label = TRUE gives us the day of the week a san ordered factor of character strings
  #such as "Sunday"
  group_by(day_of_week) |> 
  #remember na.rm means we don't consider values = N/A
  summarise(avg_delay = mean(dep_delay, na.rm = TRUE))

#I had no clue, here is the solution! 

# Function to perform bootstrap and return the standard error
bootstrap.se <- function(x) {
  B <- 200
  
  # Statistic
  S <- mean
  
  # Perform bootstrap
  S.star <- rep(0, B)
  for(b in 1:B) {
    x.star <- sample(x, replace = TRUE)
    S.star[b] <- S(x.star)
  }
  
  return(sd(S.star))
}

# Tidyverse pipeline adding calculation of the bootstrap standard error for
# each day of week grouping using above function, with a final mutate to 
# compute confidence interval
mean_ci <- flights |> 
  filter(dep_delay >= 0) |> 
  mutate(day_of_week = wday(sched_dep, label = TRUE)) |> 
  group_by(day_of_week) |> 
  summarise(mean_delay = mean(dep_delay, na.rm = TRUE),
            bootstrap_se = bootstrap.se(dep_delay)) |> 
  mutate(mean_ci_lower = mean_delay - 1.96*bootstrap_se,
         mean_ci_upper = mean_delay + 1.96*bootstrap_se)
mean_ci

library("tidyverse")
data("weather", package = "nycflights13")

#remember: sched_dep comes from exercise 5.94 
flights <- flights |>
  mutate(time_hour2 = floor_date(sched_dep, unit = "hour"))

#to check these are the same...
all(flights$time_hour == flights$time_hour2)