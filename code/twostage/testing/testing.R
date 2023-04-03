# Poisson ----------------------------------------------------------------------
# Entire trial
to_hist_antici <- c()
to_hist_react <- c()
for (i in 1:1000) {
  to_hist_antici <- c(to_hist_antici, read.csv(paste0('~/netVax/scratch/', i, '_poisanticipatory.csv'))[1,2])
  to_hist_react <- c(to_hist_react, read.csv(paste0('~/netVax/scratch/', i, '_poisreactionary.csv'))[1,2])
}
hist(to_hist_antici, col='blue', breaks=20, xlim=c(0, 0.3))
hist(to_hist_react, add=TRUE, col='red', breaks=20)
mean(to_hist_antici)
sd(to_hist_antici)
mean(to_hist_react)
sd(to_hist_react)

# Control
to_hist_antici <- c()
to_hist_react <- c()
for (i in 1:2000) {
  to_hist_antici <- c(to_hist_antici, read.csv(paste0('~/netVax/scratch/', i, '_poisanticipatory.csv'))[2,2])
  to_hist_react <- c(to_hist_react, read.csv(paste0('~/netVax/scratch/', i, '_poisreactionary.csv'))[2,2])
}
hist(to_hist_antici, col='blue', breaks=20, xlim=c(0, 0.3))
hist(to_hist_react, add=TRUE, col='red', breaks=20)
mean(to_hist_antici)
sd(to_hist_antici)
mean(to_hist_react)
sd(to_hist_react)

# Overdispersion ---------------------------------------------------------------
# Entire trial
to_hist_antici <- c()
to_hist_react <- c()
for (i in 1:1000) {
  to_hist_antici <- c(to_hist_antici, read.csv(paste0('~/netVax/scratch/', i, '_overanticipatory.csv'))[1,2])
  to_hist_react <- c(to_hist_react, read.csv(paste0('~/netVax/scratch/', i, '_overreactionary.csv'))[1,2])
}
hist(to_hist_antici, col='blue', breaks=20, xlim=c(0, 0.3))
hist(to_hist_react, add=TRUE, col='red', breaks=20)
mean(to_hist_antici)
sd(to_hist_antici)
mean(to_hist_react)
sd(to_hist_react)

# Control
to_hist_antici <- c()
to_hist_react <- c()
for (i in 1:1000) {
  to_hist_antici <- c(to_hist_antici, read.csv(paste0('~/netVax/scratch/', i, '_overanticipatory.csv'))[2,2])
  to_hist_react <- c(to_hist_react, read.csv(paste0('~/netVax/scratch/', i, '_overreactionary.csv'))[2,2])
}
hist(to_hist_antici, col='blue', breaks=20, xlim=c(0, 0.3))
hist(to_hist_react, add=TRUE, col='red', breaks=20)
mean(to_hist_antici)
sd(to_hist_antici)
mean(to_hist_react)
sd(to_hist_react)

