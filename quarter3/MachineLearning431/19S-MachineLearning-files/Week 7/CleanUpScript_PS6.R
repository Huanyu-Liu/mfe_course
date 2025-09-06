# clean-up R script for Problem Set 4

require(tm)
require(glmnet)
require(ggplot2)
require(pROC)
require(SnowballC)

# Read in the data, use the right file path for you in the below
data <- read.csv('D:/llochsto/Dropbox/Data Analytics/Data/StockNewsData/DJIA_Headline_News.csv', stringsAsFactors = FALSE)

# First, we will clean up the data and do some quick preprocessing. Let's also add a '<\\s>' token between
# the headlines. We don't want the last word of a headline and the first word of the next to be counted as a bigram.

# Make 'Date' column a Date object to make train/test splitting easier
data$Date <- as.Date(data$Date)

# Combine headlines into one text blob for each day and add sentence separation token
data$all <- paste(data$Top1, data$Top2, data$Top3, data$Top4, data$Top5, data$Top6,
                  data$Top7, data$Top8, data$Top9, data$Top10, data$Top11, data$Top12, 
                  data$Top13, data$Top14, data$Top15, data$Top16, data$Top17, data$Top18,
                  data$Top19, data$Top20, data$Top21, data$Top22, data$Top23, data$Top24,
                  data$Top25, sep=' <s> ')

# Get rid of those pesky b's and backslashes you see if you inspect the raw data
data$all <- gsub('b"|b\'|\\\\|\\"', "", data$all)

# Get rid of all punctuation except headline separators, alternative to cleaning done in tm-package
data$all <- gsub("([<>])|[[:punct:]]", "\\1", data$all)

# Reduce to only the three columns we need. 
data <- data[, c('Date', 'Label', 'all')]