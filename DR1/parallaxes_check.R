#! usr/bin/env Rscript

### Florian Glass 2017
### Monte Carlo simulation on parallaxes from De Ridder et al (2016)
### Creation of synthetic asteroseismic parallaxes matching TGAS', then adding gaussian noise to both.
### The goal is to look at the influence of negative parallaxes using an error-in-variable method.

library(leiv)  # computes the slope and intercept of a dataset using error-in-variables
library(ggplot2)
library(grid)

theme_set(theme_bw(14))  # define the "black and white" ggplot2 theme with font.size of 14 for the whole script

### parameters and data loading
iterations <- 2
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
#####################

### Creation of synthetic data:
synthetic.seismic.parallaxes <- TGAS.parallaxes
names(synthetic.seismic.parallaxes) <- c("synthetic seismic parallaxes")

synthetic.noised <- list()
TGAS.noised <- list()

### pick "iterations" element on each sample point for computation efficiency. Reshape data later.
row.number <- as.numeric(rownames(synthetic.seismic.parallaxes))
for (i in row.number){
  synthetic.noised[[i]] <- rnorm(iterations, mean=synthetic.seismic.parallaxes[i,1], sd=seismic.errors[i,1])
  TGAS.noised[[i]] <- rnorm(iterations, mean=TGAS.parallaxes[i,1], sd=TGAS.errors[i,1])
}

data.frame.synthetic.noised <- data.frame(synthetic.noised)
names(data.frame.synthetic.noised) <- row.number
data.frame.TGAS.noised <- data.frame(TGAS.noised)
names(data.frame.TGAS.noised) <- row.number

astero.synth <- data.frame(t(data.frame.synthetic.noised))
TGAS.synth <- data.frame(t(data.frame.TGAS.noised))

#################################################
### Error-in-variables for all the iterations ###
#################################################
### leiv.function buffers the data handling. We have a data.frame (synthetic data), but need to get
### subsets of different sizes from that (positive data) which data.frame doesn't accept.
leiv.function <- function(astero, TGAS){
  fit <- leiv(astero~TGAS, abs.tol=5e-6)
  return(c(fit@slope, fit@intercept))
}

fit.results <- data.frame()
fit.pos.results <- data.frame()
cat("Processing sample set:\n")
for (i in seq(iterations)){
  col.name <- paste('X', i, sep='')
  cat(col.name, '\n')
  fit.results <- rbind(fit.results,leiv.function(astero.synth[[col.name]], TGAS.synth[[col.name]]))
  ### also compute slope and intercept on positive TGAS parallax subset
  TGAS.positive <- TGAS.synth[which(TGAS.synth[[col.name]] > 0.), i]
  astero.positive <- astero.synth[which(TGAS.synth[[col.name]] > 0.), i]
  fit.pos.results <- rbind(fit.pos.results, leiv.function(astero.positive, TGAS.positive))
}
names(fit.results) <- c('slope', 'intercept')
names(fit.pos.results) <- c('slope', 'intercept')

### Whole sample
average.slope <- mean(fit.results$slope)
sd.slope <- sd(fit.results$slope)
error.slope <- sd.slope/sqrt(iterations)
average.intercept <- mean(fit.results$intercept)
sd.intercept <- sd(fit.results$intercept)
error.intercept <- sd.intercept/sqrt(iterations)

### Positive sample
average.pos.slope <- mean(fit.pos.results$slope)
sd.pos.slope <- sd(fit.pos.results$slope)
error.pos.slope <- sd.pos.slope/sqrt(iterations)
average.pos.intercept <- mean(fit.pos.results$intercept)
sd.pos.intercept <- sd(fit.pos.results$intercept)
error.pos.intercept <- sd.pos.intercept/sqrt(iterations)

####################
### Control plot ###
####################
### Plot the data. Similar to Fig. 4 in J. De Ridder et al (2016) "Asteroseismic versus Gaia distances"
X1 <- data.frame(TGAS.synth[['X1']], astero.synth[['X1']])
grob <- grobTree(textGrob(sprintf('Current sample, slope: %s, intercept: %s',round(fit.results[[1]][[1]],3), round(fit.results[[2]][[1]],3)), x=0.01,  y=0.98, hjust=0,
                          gp=gpar(col="black", fontsize=13)))
grob2 <- grobTree(textGrob(sprintf("Average all samples, slope: %s±%s, intercept: %s±%s", round(average.slope,3), round(error.slope,3), round(average.intercept,3),round(error.intercept,3)), x=0.01, y=0.95, hjust=0, 
                           gp=gpar(col="blue", fontsize=13)))
grob3 <- grobTree(textGrob(sprintf("Average all positive values, slope: %s±%s, intercept: %s±%s", round(average.pos.slope,3), round(error.pos.slope,3), round(average.pos.intercept,3), round(error.pos.intercept,3)), x=0.01, y=0.92, hjust=0,
                           gp=gpar(col="red", fontsize=13)))
p2 <-ggplot(X1, aes(x=X1[[1]], y=X1[[2]])) +
  #geom_errorbar(aes(ymin = X1[[2]]-seismic.errors, ymax = X1[[2]]+seismic.errors), col='gray') +  # should errorbars be shown ?
  #geom_errorbarh(aes(xmin = X1[[1]]-TGAS.errors, xmax = X1[[1]]+TGAS.errors), col='gray') + 
  geom_point() + 
  geom_abline(intercept=fit.results[[2]][[1]], slope=fit.results[[1]][[1]], colour='black', linetype='dashed', size=1) +
  geom_abline(intercept=average.intercept, slope=average.slope, colour='blue', size=1) +  # averaged regression of the whole iteration dataset
  geom_abline(intercept=average.pos.intercept, slope=average.pos.slope, colour='red', size=1) +
  geom_vline(xintercept=0, colour='black', linetype='dashed') + 
  xlab("synthetic TGAS parallaxes [mas]") + ylab("synthetic seismic parallaxes [mas]") + ggtitle("Synthetic data (showing one sample)") + 
  theme(plot.title=element_text(hjust=0.5)) +
  scale_colour_manual(name = '', values=c("Current sample", "All samples", "Positive samples")) +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
print(p2)

### Terminal output
cat('Error-in-variable on sample 1:\nslope:', fit.results[[1]][[1]], 'intercept:', fit.results[[2]][[1]],'\n')
cat('Error-in-variable averaged on whole sample:\nslope:', average.slope,'±',error.slope,'intercept:',average.intercept,'±',error.intercept,'\n')
cat('Error-in-variable averaged on positive samples:\nslope:', average.pos.slope,'±',error.pos.slope,'intercept:',average.pos.intercept,'±',error.pos.intercept,'\n')