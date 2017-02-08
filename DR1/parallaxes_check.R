#! usr/bin/env Rscript

### Florian Glass 2017
### Monte Carlo simulation on parallaxes from De Ridder et al (2016)
### Creation of synthetic asteroseismic parallaxes matching TGAS', then adding gaussian noise to both.
### The goal is to look at the influence of negative parallaxes using an error-in-variable method.

### TO DO:
### DONE Change synthetic data from asteroseismic to TGAS (astero have smaller errors, so it should be the basis)
### DONE Add noise to the *distance*, instead of the parallax for asteroseismic parallaxes
### (Error-in-variables on Distances vs TGAS parallaxes)
### DONE Histogram of slopes and intercepts of regression model

library(leiv)  # computes the slope and intercept of a dataset using error-in-variables
library(ggplot2)  # plotting package
library(grid)  # annotation on plot (grobTree function)

theme_set(theme_bw(14))  # define the "black and white" ggplot2 theme with font.size of 14 for the whole script

### Function definition:
leiv.function <- function(astero, TGAS){
  fit <- leiv(astero~TGAS, abs.tol=1e-6)  # linear fit: estimate the slope and intercept of a bivariate linear relationship when both variables are observed with error.
  return(c(fit@slope, fit@intercept))
}

get.csv.data <- function(){
  new.data <- read.csv(file="synthetic_data.csv", sep=',', header=TRUE)
  fit <- leiv.function(new.data$astero.synth, new.data$TGAS.synth)
  
  grob <- grobTree(textGrob(sprintf('slope: %s, intercept: %s',round(fit[1],3), round(fit[2],3)), x=0.01,  y=0.98, hjust=0,
                            gp=gpar(col="red", fontsize=13)))
  
  plot1 <- ggplot(new.data, aes(x=new.data$TGAS.synth, y=new.data$astero.synth)) +
    geom_point() + 
    geom_abline(slope=1, intercept=0, colour='gray', size=1) +
    geom_abline(slope=fit[1], intercept=fit[2], colour='red', size=1) +
    xlab("synthetic TGAS parallaxes [mas]") + ylab("synthetic seismic parallaxes [mas]") +
    annotation_custom(grob)
  print(plot1)
  return(c(new.data, fit))
}

### parameters and data loading
iterations <- 1
TGAS.offset <- 0.

DeRidder.data <- read.csv(file="KeplerGiants_inputForGeert2.txt", header=TRUE, sep=" ")
TGAS.parallaxes <- DeRidder.data[2]
TGAS.errors <- DeRidder.data[3]
seismic.parallaxes <- DeRidder.data[4]
seismic.errors <- DeRidder.data[5]

####################
### Control plot ###
####################
### plot the uncertainties vs parallaxes. Fig.3 in J.De Ridder et al (2016) "Asteroseismic versus Gaia distances"
matplot(cbind(TGAS.parallaxes[[1]], seismic.parallaxes[[1]]), cbind(TGAS.errors[[1]], seismic.errors[[1]]), pch=1,
        xlab='parallaxes [mas]', ylab='parallaxes errors [mas]')
### plot the original data
Odata <- data.frame(TGAS.parallaxes[[1]], seismic.parallaxes[[1]])
Original.fit <- leiv.function(Odata[[1]], Odata[[2]])

# Odata.pos <- data.frame(TGAS.parallaxes[which(TGAS.parallaxes>0.)], seismic.parallaxes[which(seismic.parallaxes>0.)])
# Original.pos.fit <- leiv.function(Odata.pos[[1]], Odata[[2]])

p <-ggplot(Odata, aes(x=Odata[[1]], y=Odata[[2]])) +
  geom_errorbar(aes(ymin = Odata[[2]]-seismic.errors, ymax = Odata[[2]]+seismic.errors), col='gray') + 
  geom_errorbarh(aes(xmin = Odata[[1]]-TGAS.errors, xmax = Odata[[1]]+TGAS.errors), col='gray') + geom_point() + 
  xlab("TGAS parallaxes [mas]") + ylab("seismic parallaxes [mas]") + ggtitle("Original data") + theme(plot.title=element_text(hjust=0.5)) +
  geom_abline(slope=Original.fit[[1]], intercept=Original.fit[[2]], colour='blue') +
  geom_abline(slope = 1, intercept = 0, colour = 'gray', linetype='dashed', size = 1)
print(p)
cat('Original data:\nslope:',Original.fit[1],'intercept:', Original.fit[2],'\n')
#####################

### Creation of synthetic data:
synthetic.TGAS.parallaxes <- seismic.parallaxes
names(synthetic.TGAS.parallaxes) <- c("synthetic TGAS parallaxes")

synthetic.noised <- list()
TGAS.noised <- list()

### pick "iterations" element on each sample point for computation efficiency. Reshape data later.
row.number <- as.numeric(rownames(synthetic.TGAS.parallaxes))
for (i in row.number){
  synthetic.noised[[i]] <- 1./rnorm(iterations, mean=1./seismic.parallaxes[i,1], sd=seismic.errors[i,1])  # add error on the distance. BUT error on distance won't be gaussian !
  TGAS.noised[[i]] <- rnorm(iterations, mean=synthetic.TGAS.parallaxes[i,1], sd=TGAS.errors[i,1])
}

data.frame.synthetic.noised <- data.frame(synthetic.noised)
names(data.frame.synthetic.noised) <- row.number
data.frame.TGAS.noised <- data.frame(TGAS.noised)
names(data.frame.TGAS.noised) <- row.number

astero.synth <- data.frame(t(data.frame.synthetic.noised))
TGAS.synth <- data.frame(t(data.frame.TGAS.noised))

### Add offset to TGAS.synth. Default is 0, ±0.3 is suggested too.
TGAS.synth <- TGAS.synth + TGAS.offset

### Fix column names if iterations==1
if (iterations==1){
  names(astero.synth) <- "X1"
  names(TGAS.synth) <- "X1"
}

#################################################
### Error-in-variables for all the iterations ###
#################################################
### leiv.function buffers the data handling. We have a data.frame (synthetic data), but need to get
### subsets of different sizes from that (positive data) which data.frame doesn't accept.

fit.results <- data.frame()
fit.pos.results <- data.frame()
cat("Processing sample set:\n")
for (i in seq(iterations)){
  col.name <- paste('X', i, sep='')
  cat(col.name, '\n')
  fit.results <- rbind(fit.results, leiv.function(astero.synth[[col.name]], TGAS.synth[[col.name]]))
  ### also compute slope and intercept on positive TGAS parallax subset
  TGAS.positive <- TGAS.synth[which(TGAS.synth[[col.name]] > 0.), i]
  astero.positive <- astero.synth[which(TGAS.synth[[col.name]] > 0.), i]
  fit.pos.results <- rbind(fit.pos.results, leiv.function(astero.positive, TGAS.positive))
}
names(fit.results) <- c('slope', 'intercept')
fit.results$partial <- 'all'
names(fit.pos.results) <- c('slope', 'intercept')
fit.pos.results$partial <- 'positive'

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
### Histograms of slopes and intercepts
tot.results <- data.frame(rbind(fit.results, fit.pos.results))
hist.slopes <- ggplot(tot.results, aes(x=slope, fill=partial)) +
  geom_histogram(data= subset(tot.results,partial=='all')[1], fill= 'blue', alpha=0.2) +
  geom_histogram(data= subset(tot.results,partial=='positive')[1], fill= 'red', alpha=0.2) +
  geom_vline(xintercept=average.slope, colour='blue') + geom_vline(xintercept=average.pos.slope, colour='red')
print(hist.slopes)

hist.intercepts <- ggplot(tot.results, aes(x=intercept, fill=partial)) +
  geom_histogram(data= subset(tot.results, partial=='all')[2], fill = 'blue', alpha=0.2) +
  geom_histogram(data= subset(tot.results, partial=='positive')[2], fill = 'red', alpha=0.2) +
  geom_vline(xintercept=average.intercept, colour='blue') + geom_vline(xintercept=average.pos.intercept, colour='red')
print(hist.intercepts)

### Plot the data similar to Fig. 4 in J. De Ridder et al (2016) "Asteroseismic versus Gaia distances"
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
  geom_abline(slope = 1, intercept = 0, colour = 'gray', linetype='dashed', size = 1) +
  geom_vline(xintercept=0, colour='black', linetype='dashed') + 
  xlab("synthetic TGAS parallaxes [mas]") + ylab("synthetic seismic parallaxes [mas]") + ggtitle(sprintf("Synthetic data (showing 1 sample of %s). TGAS offset: %s mas",iterations,TGAS.offset)) + 
  theme(plot.title=element_text(hjust=0.5)) +
  scale_colour_manual(name = '', values=c("Current sample", "All samples", "Positive samples")) +
  annotation_custom(grob) + annotation_custom(grob2) + annotation_custom(grob3)
print(p2)

cat('TGAS offset:',TGAS.offset,'mas\n')
cat('Error-in-variable on sample 1:\nslope:', fit.results[[1]][[1]], 'intercept:', fit.results[[2]][[1]],'\n')
cat('Error-in-variable averaged on whole sample:\nslope:', average.slope,'±',error.slope,'intercept:',average.intercept,'±',error.intercept,'\n')
cat('Error-in-variable averaged on positive samples:\nslope:', average.pos.slope,'±',error.pos.slope,'intercept:',average.pos.intercept,'±',error.pos.intercept,'\n')