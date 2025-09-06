#
# code snippets for Topic 7
#

rm(list = ls())

# set your working directory, where you have downloaded the Stata data file
# change the below to match your folder
setwd("D:/llochsto/Dropbox/Data Analytics/Data")

# need text mining package
require(NLP)
require(tm)

# Read in text files
text_dir <- file.path("D:/llochsto/Dropbox/Data Analytics/Data", "TextData")
main_corpus <- Corpus(DirSource(text_dir))
summary(main_corpus)

# let's process our text-file

# remove punctuation
main_corpus <- tm_map(main_corpus, removePunctuation)

# remove "\t" (example of specific use of gsub, and everything that is not a number or letter or space
# [^a-zA-Z0-9 ] is an example of a Regular Expression 
main_corpus <- tm_map(main_corpus, content_transformer(gsub), pattern = "\t", replacement = " ")
main_corpus <- tm_map(main_corpus, content_transformer(gsub), pattern = "[^a-zA-Z0-9 ]", replacement = " ")

# remove numbers (we have those better represented in CompuStat)
main_corpus <- tm_map(main_corpus, removeNumbers)

# convert all to lowercase, so word is recognized with arbitrary capitalization
main_corpus <- tm_map(main_corpus, tolower)

# remove particular words that you know are irrelevant noise
main_corpus <- tm_map(main_corpus, removeWords, c("table of contents","sec","securities exchange commission", "united states"))
main_corpus <- tm_map(main_corpus, removeWords, c("company","company's","financial", "september", "net", "securities", "including", "inc", "billion", "million", "assets", "operating", "statements", "tax"))
main_corpus <- tm_map(main_corpus, removeWords, c("may","notes","can","changes","cost","will","also","rate", "rates","equity","available","certain","results","relative"))

# remove "stopwords" (e.g., and, to, a, as, the, ...)
main_corpus <- tm_map(main_corpus, removeWords, stopwords("english"))

# stemming words, i.e., keep only the stem so as not to differentially count investing, invest, invests
# taking out common word endings such as 'ing', 'es', and 's'
require(SnowballC)
main_corpus <- tm_map(main_corpus, stemDocument)

# finally, let's get rid of all the extra white space in the document so all words are only 
# separated by one space
main_corpus <- tm_map(main_corpus, stripWhitespace)


# if you want to save as one big corpus, do the below
# main_corpus <- tm_map(main_corpus, PlainTextDocument)
# writeLines(as.character(main_corpus), con="TextData/mycorpus.txt")


# next, organize words into document term matrix, which then can be used for analysis
corpus_matrix <- DocumentTermMatrix(main_corpus)
inspect(corpus_matrix[1:10,1:10])

# organize words by frequency
freq <- colSums(as.matrix(corpus_matrix))
ord_corpus <- order(freq)

# see most common words
freq[tail(ord_corpus)]

#  identify words that appear frequently
# create a convenient data.frame
word_freq <- data.frame(word=names(freq), freq=freq)

# plot most frequent words along with frequency
require(ggplot2)
p <- ggplot(subset(word_freq, freq>2100), aes(word, freq))
p <- p + geom_bar(stat = "identity")
p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
p

# plot wordcloud, a net way to express the data, size of font relates to frequency
require(wordcloud)
# plot 100 most frequent words, and add some color
wordcloud(names(freq), freq, max.words=100, rot.per = 0.2, colors=brewer.pal(6, "Dark2"))

# so far, mainly a word salad... How can we extract information from this?


# let's search for frequency of a pre-set list of words we perceive as related to high growth expectations
freq[c("invest","growth","grow","high","strong","lead","good")]

freq[c("loss","weak","low","poor","uncertain","under","disappoint")]


# let's loop through the 
higrowth_words = NULL
logrowth_words = NULL

for (j in 1:23)
{
  corpus_matrix <- DocumentTermMatrix(main_corpus[j])
  
  # organize words by frequency
  freq <- colSums(as.matrix(corpus_matrix))  
  higrowth_words = rbind(higrowth_words,freq[c("invest","growth","grow","high","strong","lead")])
  
  logrowth_words = rbind(logrowth_words,freq[c("loss","weak","low","poor","uncertain","under")])
}

# add rows to get a score by year
hi_growth <- rowSums(higrowth_words, na.rm = TRUE)
lo_growth <- rowSums(logrowth_words, na.rm = TRUE)

# display scores
hi_growth
lo_growth

# create hi over lo metric
growth_score = hi_growth / lo_growth
growth_score

# plot growth score versus year
year <- c(1994:2016)
qplot(year,growth_score, col=I("blue"), geom = "line") + theme_bw()


# plot growth score in year t versus lnROE of Apple in year t + 1
aapl_lnROE <- c(0.0324513, 0.1050476, 0.1230894, -0.2899466, -.5851619,.2052859,.2108069,.1922505,-.0059545,.0151465,.016474,.0559076,.1884116,.2130643,.2623312,.2519465,.2787039,.3133555,.3649631,.3712058,.2820229)
reg_data <- data.frame(year = year[1:18], growth_score = growth_score[1:18] / 2 - 0.5, series  = "growth_score")
reg_data <- rbind(reg_data, data.frame(year = year[1:18], growth_score = aapl_lnROE[2:19], series  = "aapl_lnROE"))
qplot(year,growth_score,data=reg_data,facets=series~.,col=I("blue"), geom="line") + theme_bw()

# regress to see if a statistically significant forecasting relation
reg <- lm(aapl_lnROE[2:19]~growth_score[1:18])
summary(reg)





# LDA topic allocation modeling, Topic 7b
# Clean workspace
rm(list = ls())

# Load packages
library(data.table)
library(tm)
library(glmnet)
library(ggplot2)
library(pROC)
library(SnowballC)
library(dplyr)
library(wordcloud)
library(ROCR)
library(stargazer)
library(sandwich)
library(lmtest)
library(lda)


# Load data
data <- read.csv('DJIA_Headline_News.csv', stringsAsFactors = FALSE)

# just learn topics of top 2 headlines per day (to save time)
data <- rbind(data$Top1, data$Top2)

## Get rid of those pesky b's and backslashes you see if you inspect the raw data
data <- gsub('b"|b\'|\\\\|\\"', "", data)

## Get rid of all punctuation except headline separators, alternative to cleaning done in tm-package
data <- gsub("([<>])|[[:punct:]]", "\\1", data)

# lda routine requires documents to be in specific list format:
# A list whose length is equal to the number of documents, D. Each element of
# documents is an integer matrix with two rows.  Each column of
# documents[[i]] (i.e., document i) represents a word occurring in the document.
# documents[[i]][1, j] is a 0-indexed word identifier for the jth word in document
# i.  That is,  this should be an index - 1 into vocab. documents[[i]][2,  j] is an integer specifying the number of
# times that word appears in the document.

# create large list for lda routine, split (tokenize) on space
doc.list <- strsplit(data, "[[:space:]]+")
doc.list <- sapply(doc.list, tolower)

# create a table of terms
term.table <- table(unlist(doc.list))
term.table <- sort(term.table, decreasing = TRUE)
head(term.table)

# remove terms that are stop_words or occur less than 100 times
del <- names(term.table) %in% stopwords(kind = 'en') | term.table < 2
term.table <- term.table[!del]
head(term.table)

# list of words in corpus, for use by lda
vocab = names(term.table)

# define function that helps getting the text data in right format for lda
get.terms <- function(x) {
  index <- match(x, vocab)
  index <- index[!is.na(index)]
  rbind(as.integer(index - 1), as.integer(rep(1,length(index))))
}

# create list version of data
documents <- lapply(doc.list, get.terms)

# Compute some statistics related to the data set
D <- length(documents)  # number of documents (1989)
W <- length(vocab)  # number of terms in the vocab (8848)
doc.length <- sapply(documents, function(x) sum(x[2, ]))  # number of tokens per document [11,3, 4...]
head(doc.length,20)

N <- sum(doc.length)  # total number of tokens/terms in the data
term.frequency <- as.integer(term.table) # Frequency of each of the terms
term.frequency

# lda model tuning parameters
K <- 5
G <- 5000     # Number of iterations to arrive at convergence 
alpha <- 0.02
eta <- 0.02

# Fit the model
set.seed(357)
# Capture start time from your system
t1 <- Sys.time()

# Begin lda execution
lda_fit <- lda.collapsed.gibbs.sampler(documents = documents, K = K, vocab = vocab, 
                                       num.iterations = G, alpha = alpha, 
                                       eta = eta, initial = NULL, burnin = 0,
                                       compute.log.likelihood = TRUE)

# Capture end-time from your system
t2 <- Sys.time()
# Notice the time it took to complete lda on this data (Use this technique whenever you want to check the execution time for an algorithm/piece of code)
t2 - t1  # outputs time it took to fit model


# show top topic words
top.words <- top.topic.words(lda_fit$topics, num.words = 10, by.score = FALSE)
top.words

# show top documents per topic
top.documents <- top.topic.documents(lda_fit$document_sums, num.documents = 10, alpha = alpha)
top.documents

# look at top documents in each topic (most representative)
data[top.documents[1,1]]
data[top.documents[1,2]]
data[top.documents[1,3]]
data[top.documents[1,4]]
data[top.documents[1,5]]


# get the predictive model probability that a word will show up in document
pred.dist <- predictive.distribution(lda_fit$document_sums, lda_fit$topics, alpha = alpha, eta = eta)
# words in rows, document in columns
pred.dist[1:10,1:8]

# note that this can easily be used as a similarity measure to gauge similarity of two documents
# e.g.
sim_measure12 <- t(pred.dist[,1]) %*% pred.dist[,2]
sim_measure12
sim_measure13 <- t(pred.dist[,1]) %*% pred.dist[,3]
sim_measure13

# get probability that each document belongs to a certain topic
## Number of documents to display
N <- 5
topic.proportions <- t(lda_fit$document_sums) / colSums(lda_fit$document_sums)
topic.proportions <- topic.proportions[sample(1:dim(topic.proportions)[1], N),]
topic.proportions[is.na(topic.proportions)] <-  1 / K
colnames(topic.proportions) <- apply(top.words, 2, paste, collapse=" ")
topic.proportions.df <- melt(cbind(data.frame(topic.proportions),document=factor(1:N)),variable.name="topic",id.vars = "document")  
head(topic.proportions.df)




