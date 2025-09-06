# Chady Gemayel
# 4/4/2019
# A (very) brief introduction to data.table



# Load packages
library(data.table)
library(zoo)
rm(list=ls())



# Load data
Chicken = data.table(ChickWeight)
Chicken[, Chick := as.numeric(levels(Chick))[Chick]] # I don't like working with factors. This is not necessary for the code below to work.
Chicken[, Diet := as.numeric(levels(Diet))[Diet]] # I don't like working with factors. This is not necessary for the code below to work.



# Average weight of chickens by age
Q1 = Chicken[, .(Average_weight = mean(weight, na.rm = T)), by = Time] # Define a new datatable with column Time and Average_weight, which calculates the average weight for each Time.



# Average change in weight in the first 10 days
## Assume we are only interested in 2 day changes
setorder(Chicken, Chick, Time) # Sort data, to make sure I use correct lag
Chicken[, delta_weight := weight - shift(weight), by = Chick] # Calculate the change in weight between this row and the last, by Chick
Chicken[, delta_Time := Time == shift(Time) + 2, by = Chick]# Calculate the change in Time between this row and the last, by Chick
Q2 = Chicken[delta_Time == T & Time <= 10, .(Average_weight_change = mean(delta_weight)), by = Chick] # Define a new datatable with column Chick and Average_weight_change, which calculates the average 2-day change in weight over the first ten days for each Chick.



# Standard deviation in weight for the first 20 days of each chicken's life
Chicken[, weight_sigma := shift(rollapplyr(weight, 11, sd, fill = NA))] # Not done by chick! So if we have unbalanced data, or look at the wrong value, this would not be doing what we want!
setorder(Chicken, Chick, Time)
Q3 = unique(Chicken, by = c('Chick'), fromLast = T)



# Fill in missing days in chicks
Chicken[, time_diff := as.integer(Time - shift(Time)), by = Chick] # Difference in time between current and previous observation by chick
incomplete_Chicken = unique(Chicken[time_diff != 1, .(Chick)]) # Mark all chicks that have time gaps
setkey(incomplete_Chicken, Chick)
Time = seq(from = min(Chicken$Time, na.rm = T), to = max(Chicken$Time, na.rm = T))
filled = incomplete_Chicken[CJ(unique(incomplete_Chicken$Chick), Time), allow.cartesian = T, roll = -Inf] # Build a datatable with all possible times
Chicken = merge(Chicken, filled, by.x = c('Chick', 'Time'), by.y = c('Chick', 'V2'), all = T) # Merge in the filled dataset into the original
setkey(Chicken)
setorder(Chicken, Chick, Time)




# How to connect to WRDS, download dataset, and disconnect
## Instructions: https://wrds-www.wharton.upenn.edu/pages/support/programming-wrds/programming-r/r-wrds-cloud/

## Load packages
library(RPostgres)

## Function that connects to WRDS
wrds <- dbConnect(Postgres(), host = 'wrds-pgdata.wharton.upenn.edu', port = 9737, user = 'Your username', password = 'Your password', sslmode = 'require', dbname = 'wrds')

## Import some Compustat data
res = dbSendQuery(wrds, "select gvkey, datadate, fyear, EXCHG from COMPM.FUNDA where consol = 'C' and indfmt = 'INDL' and datafmt = 'STD' and popsrc = 'D'")
Cstat = data.table(dbFetch(res))
dbClearResult(res)

## Disconnect from WRDS
dbDisconnect(wrds)
rm(res, wrds)