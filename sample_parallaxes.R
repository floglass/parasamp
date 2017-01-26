#! usr/bin/env Rscript

#############
### TO DO ###
#############
### 3)Plot histogram of samples
#############

### Run with: Rscript sample_parallaxes.R from terminal

### Create synthetic data of a cluster at distance "true.distance", with a
### relative standard deviation of "true.relative.sigma".

### The code computes three different distances from sample as shown in the diagram below.
### parallax [arcsecond] = 1/distance [parsec]
### noise is gaussian, with a standard deviation of "true.sigma"
###
##                     +--------+
##                     |Distance|
##                     +---+----+
##                         |
##                         |
##                     +---v----+
##                     |Parallax|
##                     +---+----+
##                         |
##                         |
##                    +----v----+
##         +----------+Add Noise+--------------------------------------+
##         |          +----+----+                                      |
##         |               +--------+                                  |
##         |                        |                                  |
## +-------v------------------+   +-v--------------+   +---------------v---------+
## | Average positive Distance|   |Average Parallax|   |Average positive Parallax|
## +--------------------------+   +-+--------------+   +---------------+---------+
##                                  |                                  |
##                                  |                                  |
##                                +-v------+                      +----v---+
##                                |Distance|                      |Distance|
##                                +--------+                      +--------+
##
## code by Florian Glass, 2016

library(ggplot2)
library(reshape2)

### Model parameters
### General
iterations <- 10000  # number of iteration
sample.size <- 100  # size of sample
cat("Iterations:", iterations, "\n")
cat("Samples size:", sample.size, "\n")

### Specific
true.distance <- 134  # pc -- 134 pc is distance to Pleiades
true.relative.sigma <- 0.75  # sigma_parallax/parallax
true.parallax <- 1/true.distance  # arcsec
true.sigma <- true.relative.sigma * true.parallax
cat("Relative parallax error: <", true.relative.sigma, "\n")
###

#### Main loop
count <- 0
cat("Main loop..\n")

average <- vector("numeric", iterations)
median <- vector("numeric", iterations)
average.distance <- vector("numeric", iterations)
distance.from.average.parallax <- vector("numeric", iterations)
average.positive.parallax <- vector("numeric", iterations)
distance.from.positive.parallax <- vector("numeric", iterations)
negative.parallaxes.number <- vector("numeric", iterations)

while (count < iterations){
  count <- count + 1
  pi <- rnorm(sample.size, mean = true.parallax, sd = true.sigma)  # parallaxes + noise
  average[[count]] <- mean(pi)  # average parallax of the sample
  median[[count]] <- median(pi)  # median parallax of the sample

  ## distance from parallax sample
  distances <- vector("numeric", sample.size)
  for (i in seq_along(pi)){
    if (pi[[i]]>0.){
      distances[[i]] <- 1/pi[[i]]
    }
  }
  distances <- distances[which(distances!=0)]
  average.distance[[count]] <- mean(distances)  # average distance of the draws

  ## distance from average parallax
  distance.from.average.parallax[[count]] <- 1/average[[count]]

  ## get positive parallaxes
  positive.parallaxes <- vector("numeric", sample.size)
  neg.par <- 0
  for (i in seq_along(pi)){
    if (pi[[i]] > 0.){
      positive.parallaxes[[i]] <- pi[[i]]
    } else {
      neg.par <- neg.par + 1
    }
  }
  negative.parallaxes.number[[count]] <- neg.par
  positive.parallaxes <- positive.parallaxes[which(!is.na(positive.parallaxes))]
  average.positive.parallax[[count]] <- mean(positive.parallaxes)

  ## distance from average positive parallax
  distance.from.positive.parallax[[count]] <- 1/average.positive.parallax[[count]]
}
cat("Done..\n")
cat("----------\n")

### Compute the average of all three distances
final.average.distance <- mean(average.distance)
final.sd.average.distance <- sd(average.distance)/sqrt(iterations)

final.dist.parallaxes <- mean(distance.from.average.parallax)
final.sd.dist.parallaxes <- sd(distance.from.average.parallax)/sqrt(iterations)

final.dist.pos.parallaxes <- mean(distance.from.positive.parallax)
final.sd.dist.pos.parallaxes <- sd(distance.from.positive.parallax)/sqrt(iterations)
average.negative.parallaxes <- mean(negative.parallaxes.number)  # per sample

### Print results (True.distance, average.distance, distance.from.parallaxes)
cat("True Distance:", true.distance, "\n")
cat("Average positive distance:", final.average.distance, "+/-",
    final.sd.average.distance, "\n")
cat("Distance from average parallaxes:", final.dist.parallaxes, "+/-",
    final.sd.dist.parallaxes, "\n")
cat("Distance from average positive parallaxes:", final.dist.pos.parallaxes,
    "+/-", final.sd.dist.pos.parallaxes, "\n")
cat("Average negative parallaxes per sample:",
            average.negative.parallaxes, "\n")
cat("----------\n")

avg.paral <- melt(distance.from.average.parallax)  # create data.frame from vector
pos.paral <- melt(distance.from.positive.parallax)
avg.dist <- melt(average.distance)
hist.data <- cbind(avg.paral, pos.paral, avg.dist)  # combine two data.frames column-wise
names(hist.data) <- c("dist from average parallaxes", "dist from positive parallaxes", "average distances")
write.csv(hist.data, file="distances.csv", row.names=FALSE)

### Plot histograms
### ggplot2 only works with data.frame --> has to create some with library reshape2 melt function
### R is really not practical for plotting. --> better use SM or python.
# mean.average.parallaxes <- mean(hist.data$`average parallaxes`)
# mean.positive.parallaxes <- mean(hist.data$`positive parallaxes`)
# p <- ggplot(data=hist.data)+
#   geom_histogram(aes(x=hist.data$`average parallaxes`),color="black", fill="blue", alpha=0.2)+
#   geom_histogram(aes(x=hist.data$`positive parallaxes`),color="black", fill="red", alpha=0.2)+
#   geom_vline(aes(xintercept=mean.average.parallaxes),color="blue", linetype="dashed", size=1)+
#   geom_vline(aes(xintercept=mean.positive.parallaxes),color="red", linetype="dashed", size=1)+
#   xlab("distance")+
#   geom_text(aes(mean.average.parallaxes,0),label = round(mean.average.parallaxes,2), vjust = 1.5, hjust=-0.1, colour="blue")+
#   geom_text(aes(mean.positive.parallaxes,0), label=round(mean.positive.parallaxes,2), vjust=1.5, hjust=1.1, colour="red")+
#   #scale_colour_manual(values =c('blue','red'), labels = c('average parallaxes','positive parallaxes'))+
#   #guides(color=guide_legend("Computed distances from:"))
# print(p)
