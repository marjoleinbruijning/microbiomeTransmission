# Simulate selection on vertical microbiome transmission

R-code to perform analysis as presented in the manuscript: **Natural selection for imprecise vertical transmission in host-microbiota systems** (Bruijning *et al.*, in prep).


See [here](https://marjoleinbruijning.shinyapps.io/simulhostmicrobiome/) for an interactive tool.

## Example code

The example code to run simulations can be found in R/scenarios.R. Start by loading required packages and functions.

```r

# Load functions
source('R/functions.R')

```

Then create a data frame that includes the settings of all simulations you want to run. An example is given below:

```r
# Settings
# See runSimul.R for default settings of other variables
scenarios <- expand.grid(varGm=c(0.01),
                         tau=seq(0,1,length.out=10),
                         varEnv=c(0,2),
                         timeWithin=1,
                         selcoef=1,
                         nhost=500)
reps <- 10 # number of replicates
scenarios <- scenarios[rep(seq_len(nrow(scenarios)),reps),]
scenarios$rep <- rep(1:reps,each=nrow(scenarios)/reps)
scenarios$setseedEnv <- scenarios$rep  # each replicate different seed
```

Note that in file runSimul.R the default settings of some other variables are given. These can also be added to the data frame (scenarios) if they need to be altered.

Finally, run the simulation. Code below illustrates how this could be done using a for loop.

```r
## Run simulations (can take a while)
out <- matrix(NA,nrow=nrow(scenarios),ncol=3)
colnames(out) <- c('meanfit','varpheno','phenodev')

exl <- 1:500 # remove first 100 timesteps

for(r in 1:nrow(scenarios)) {
  source('R/runSimul.R')
  out[r,'meanfit'] <- mean(log(fitness[-exl]))
  out[r,'varpheno'] <- mean(varpheno[-exl])
  out[r,'phenodev'] <- mean(colMeans(pheno[,-exl])^2)  
}

```

Results can be plotted using:

```r
## Plot results
par(mfrow=c(1,3))
cols <- rainbow(2,alpha=.6)

mu <- tapply(out[,'varpheno'],list(scenarios$tau,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1],bg=cols[1],
     xlab='Transmission fidelity',ylab='Phenotypic variance',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2],bg=cols[2],
       pch=21,cex=2)

mu <- tapply(out[,'phenodev'],list(scenarios$tau,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1],bg=cols[1],ylim=c(0,0.3),
     xlab='Transmission fidelity',ylab='Average deviation from P=0',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2],bg=cols[2],
       pch=21,cex=2)

mu <- tapply(out[,'meanfit'],list(scenarios$tau,scenarios$varEnv),mean)
plot(as.numeric(rownames(mu)),mu[,1]-mu[1,1],bg=cols[1],ylim=c(-.3,0.5),
     xlab='Transmission fidelity',ylab='Relative long-term fitness',pch=21,cex=2)
points(as.numeric(rownames(mu)),mu[,2]-mu[1,2],bg=cols[2],
       pch=21,cex=2)
```       
