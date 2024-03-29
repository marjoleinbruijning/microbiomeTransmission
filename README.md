# Simulate selection on vertical microbiome transmission

R-code to perform analysis as presented in the manuscript: **When the microbiome defines the host phenotype: selection on transmission in varying environments.** (Bruijning *et al.*, in prep).


See [here](https://marjoleinbruijning.shinyapps.io/simulhostmicrobiome/) for an interactive tool.

## Example code

The example code to run simulations can be found in R/scenarios.R. Start by loading required packages and functions.

```r
# Load packages
require(parallel)

# Load functions
source('R/functions.R')

```

Then create a data frame that includes the settings of all simulations you want to run. An example is given below:

```r
# Settings
# See runSimul.R for default settings of other variables
scenarios <- expand.grid(varGm=c(0.01),
                         minTrans=seq(0,1,length.out=9),
                         varEnv=c(0,2),
                         timeWithin=1,
                         selcoef=1,
                         nhost=100)
scenarios$maxTrans <- scenarios$minTrans # transmission same for all hosts and microbes
reps <- 10 # number of replicates
scenarios <- scenarios[rep(seq_len(nrow(scenarios)),reps),]
scenarios$rep <- rep(1:reps,each=nrow(scenarios)/reps)
scenarios$setseedEnv <- scenarios$rep  # each replicate different seed
```

Note that in file runSimul.R the default settings of some other variables are given. These can also be added to the data frame if they need to be altered.

Finally, run the simulation. Code below illustrates how this could be done using a for loop, or in parallel.

```r
## Run simulations (can take a while)

## Parallel
out <- unlist(mclapply(1:nrow(scenarios), function(r) {
  source('R/runSimul.R',local=TRUE)
  meanfit <- mean(log(fitness[(timeBetween*0.5):(timeBetween-1)]))
  varpheno <- mean(varpheno[(timeBetween*0.5):(timeBetween-1)])
  phenodev <- mean(colMeans(pheno[,(timeBetween*0.5):(timeBetween-1)])^2)
  return(cbind(meanfit=meanfit,varpheno=varpheno,phenodev=phenodev))
}, mc.cores=3)) # number of cores
out <- matrix(out,nrow=nrow(scenarios),byrow=TRUE)
colnames(out) <- c('meanfit','varpheno','phenodev')

## For loop
out <- matrix(NA,nrow=nrow(scenarios),ncol=3)
colnames(out) <- c('meanfit','varpheno','phenodev')

for(r in 1:nrow(scenarios)) {
  source('R/runSimul.R')
  out[r,1] <- mean(log(fitness[(timeBetween*0.5):(timeBetween-1)]))
  out[r,2] <- mean(varpheno[(timeBetween*0.5):(timeBetween-1)])
  out[r,3] <- mean(colMeans(pheno[,(timeBetween*0.5):(timeBetween-1)])^2)  
}

```

Results can be plotted using:

```r
## Plot results
par(mfrow=c(1,3))
cols <- rainbow(2,alpha=.6)

mu <- tapply(out[,'varpheno'],list(scenarios$minTrans,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1],bg=cols[1],
     xlab='Transmission fidelity',ylab='Phenotypic variance',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2],bg=cols[2],
       pch=21,cex=2)

mu <- tapply(out[,'phenodev'],list(scenarios$minTrans,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1],bg=cols[1],ylim=c(0,0.3),
     xlab='Transmission fidelity',ylab='Average deviation from P=0',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2],bg=cols[2],
       pch=21,cex=2)

mu <- tapply(out[,'meanfit'],list(scenarios$minTrans,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1]-mu[1,1],bg=cols[1],ylim=c(-.3,0.5),
     xlab='Transmission fidelity',ylab='Relative long-term fitness',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2]-mu[1,2],bg=cols[2],
       pch=21,cex=2)
```       
