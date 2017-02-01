#! usr/bin/env Rscript

### Florian Glass 2017
### Monte Carlo simulation on parallaxes from De Ridder et al (2016)
### Creation of synthetic asteroseismic parallaxes matching TGAS', then adding gaussian noise to both.
### The goal is to look at the influence of negative parallaxes using an error-in-variable method.

library(leiv)  # computes the slope and intercept of a dataset using error-in-variables
library(ggplot2)

theme_set(theme_bw(14))  # define the "black and white" ggplot2 theme with font.size of 14 for the whole script

###########################
### Function definition ###

#' Remove negative parallaxes from data.frame
#' 
#' @param to.subset The data.frame to subset
#' @return subsetted data.frame
remove.negatives <- function(to.subset){
  if (length(to.subset) > 1){
    to.subset <- to.subset[[1]]
  }
  neg.para <- which(to.subset < 0.)
  to.subset.pos <- to.subset[-c(neg.para),1]
  return(to.subset.pos)
}

### parameters and data handling
iterations <- 3
DeRidder.data <- read.csv(file="KeplerGiants_inputForGeert2.txt", header=TRUE, sep=" ")
TGAS.parallaxes <- DeRidder.data[2]
TGAS.errors <- DeRidder.data[3]
seismic.parallaxes <- DeRidder.data[4]
seismic.errors <- DeRidder.data[5]

#####################
### Control plots ###
#####################
### plot the uncertainties vs parallaxes. Fig.3 in J.De Ridder et al (2016) "Asteroseismic versus Gaia distances"
matplot(cbind(TGAS.parallaxes[[1]], seismic.parallaxes[[1]]), cbind(TGAS.errors[[1]], seismic.errors[[1]]), pch=1,
        xlab='parallaxes [mas]', ylab='parallaxes errors [mas]')
### plot the original data
Odata <- data.frame(TGAS.parallaxes[[1]], seismic.parallaxes[[1]])
p <-ggplot(Odata, aes(x=Odata[[1]], y=Odata[[2]])) +
  geom_errorbar(aes(ymin = Odata[[2]]-seismic.errors, ymax = Odata[[2]]+seismic.errors), col='gray') + 
  geom_errorbarh(aes(xmin = Odata[[1]]-TGAS.errors, xmax = Odata[[1]]+TGAS.errors), col='gray') + geom_point() + 
  xlab("TGAS parallaxes [mas]") + ylab("seismic parallaxes [mas]") + ggtitle("Original data") + theme(plot.title=element_text(hjust=0.5))
print(p)
############

### main loop
synthetic.seismic.parallaxes <- TGAS.parallaxes
names(synthetic.seismic.parallaxes) <- c("synthetic seismic parallaxes")

#noised <- vector('numeric', nrow(synthetic.seismic.parallaxes))
synthetic.noised <- list()
TGAS.noised <- list()

### pick "iterations" element on each sample point for computation efficiency.
row.number <- as.numeric(rownames(synthetic.seismic.parallaxes))
for (i in row.number){
  synthetic.noised[[i]] <- rnorm(iterations, mean=synthetic.seismic.parallaxes[i,1], sd=seismic.errors[i,1])
  TGAS.noised[[i]] <- rnorm(iterations, mean=TGAS.parallaxes[i,1], sd=TGAS.errors[i,1])
}

### Need to get the 100 iterations of 939 samples, NOT 939 samples of 100 iterations as the output of the for loop
data.frame.synthetic.noised <- data.frame(synthetic.noised)
names(data.frame.synthetic.noised) <- row.number
data.frame.TGAS.noised <- data.frame(TGAS.noised)
names(data.frame.TGAS.noised) <- row.number

### Create data.frame containing 100 iterations of the sample of 939 objects, also select a positive subset
astero.synth <- data.frame(t(data.frame.synthetic.noised))
TGAS.synth <- data.frame(t(data.frame.TGAS.noised))

#################################################
### Error-in-variables for all the iterations ###
#################################################
### Create fit.results data.frame 
fit.slope <- vector('numeric', iterations)
fit.intercept <- vector('numeric', iterations)
cat("Processing column:\n")
for (i in seq(iterations)){
  col.name <- paste('X', i, sep='')
  cat(col.name, '\n')
  fit <- leiv(astero.synth[[col.name]]~TGAS.synth[[col.name]], abs.tol=1e-5)  # must specify abs.tol, since data is nearly linear (see leiv documentation)
  fit.slope[[i]] <- fit@slope
  fit.intercept[[i]] <- fit@intercept
}
#################################################

fit.results <- data.frame(cbind(fit.slope, fit.intercept))
names(fit.results) <- c('slope', 'intercept')

average.slope <- mean(fit.results$slope)
sd.slope <- sd(fit.results$slope)
error.slope <- sd.slope/sqrt(iterations)
average.intercept <- mean(fit.results$intercept)
sd.intercept <- sd(fit.results$intercept)
error.intercept <- sd.intercept/sqrt(iterations)

####################
### Control plot ###
####################
### Plot the data. Similar to Fig. 4 in J. De Ridder et al (2016) "Asteroseismic versus Gaia distances"
X1 <- data.frame(TGAS.synth[['X1']], astero.synth[['X1']])
p2 <-ggplot(X1, aes(x=X1[[1]], y=X1[[2]])) +
  #geom_errorbar(aes(ymin = X1[[2]]-seismic.errors, ymax = X1[[2]]+seismic.errors), col='gray') + 
  #geom_errorbarh(aes(xmin = X1[[1]]-TGAS.errors, xmax = X1[[1]]+TGAS.errors), col='gray') + 
  geom_point() + 
  geom_abline(intercept=fit.results[[2]][[1]], slope=fit.results[[1]][[1]], colour='red', linetype='dashed', size=1) +
  geom_abline(intercept=average.intercept, slope=average.slope, colour='blue', size=1) +  # averaged regression of the whole iteration dataset
  geom_vline(xintercept=0, colour='black', linetype='dashed') + 
  #geom_hline(yintercept=fit.results[[1]][[2]], colour='black') +
  xlab("synthetic TGAS parallaxes [mas]") + ylab("synthetic seismic parallaxes [mas]") + ggtitle("Synthetic data") + 
  theme(plot.title=element_text(hjust=0.5))
print(p2)
cat('Error-in-variable for sample X1:\nslope:', fit.results[[1]][[1]], 'intercept:', fit.results[[2]][[1]],'\n')
############

cat('Error-in-variable for whole sample:\nslope:', average.slope,'+/-',error.slope,'intercept:',average.intercept,'+/-',error.intercept,'\n')