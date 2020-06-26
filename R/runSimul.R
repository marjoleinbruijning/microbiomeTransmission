
####################################################################################
########## Code for manuscript: When the microbiome defines the host ###############
######### phenotype: selection on transmission in varying environments #############
####################################################################################
############################## Marjolein Bruijning #################################
####################################################################################
################################### 2020-05-27 #####################################
####################################################################################

###################### Code to run one specific simulation #########################

# host
timeBetween <- 1000 # number of host generations

if('nhost' %in% colnames(scenarios)) {
  nhost <- scenarios$nhost[r]
} else {nhost <- 500} # Number of different hosts

## microbes
timeWithin <- scenarios$timeWithin[r] # timesteps within host generation
nmicrobe <- 100 # Number of microbes within host
deathrate <- 1 # rate at which microbes die (when 1, corresponds to 1 generation)
if('imm' %in% colnames(scenarios)) { # new immigrant probabilitiy ('horizontal transmission')
  migr <- scenarios$imm[r]
} else {migr <- 0}
smicrobe <- 1000
selcoef <- scenarios$selcoef[r]

# genetic/environmental variances
if('varGh' %in% colnames(scenarios)) {
  varGh <- scenarios$varGh[r]
} else {varGh <- 0}
varGm <- scenarios$varGm[r]
if('varRes' %in% colnames(scenarios)) {
  varE <- scenarios$varRes[r]
} else {varE <- 0}

# draw host and microbe genetic values
Gmicro <- rnorm(smicrobe,0,sqrt(varGm))
Gmicro <- Gmicro - mean(Gmicro) # mean at exactly 0
Ghost <- rnorm(nhost,0,sqrt(varGh))

# transmission fidelity
# same within each host
trans <- matrix(runif(nhost,scenarios$minTrans[r],scenarios$maxTrans[r]),
                nrow=nhost,ncol=smicrobe)

# colonization prob
if('colonize' %in% colnames(scenarios)) {
  tmp <- scenarios$colonize[r]
} else {tmp <- 0}
colonize <- rbeta(smicrobe,1,tmp) ; colonize <- colonize/sum(colonize)

# proliferation prob
if('vgrowth' %in% colnames(scenarios)) {
  tmp <- scenarios$vgrowth[r]
} else {tmp <- 0}
growth <- rbeta(smicrobe,1,tmp) ; growth <- growth/sum(growth)

# sexual reproduction
if('sex' %in% colnames(scenarios)) {
  sex <- scenarios$sex[r]
} else {sex <- FALSE}

# set mean phenotype to 0 every generation
if('correctP' %in% colnames(scenarios)) {
  correctP <- scenarios$correctP[r]
} else {correctP <- FALSE}

# set phenotypic variance to 0 every generation
if('correctVar' %in% colnames(scenarios)) {
  correctVar <- scenarios$correctVar[r]
} else {correctVar <- FALSE}

## Simulate environments
# temporal autocorrelation
if('corrEnv' %in% colnames(scenarios)) {
  corrEnv <- scenarios$corrEnv[r]
} else {corrEnv <- 1E-8}

if('setseedEnv' %in% colnames(scenarios)) {
  set.seed(scenarios$setseedEnv[r])
} else {set.seed(as.numeric(Sys.time()))}
allEnv <- arima.sim(list(order=c(1,0,0), ar=corrEnv), n=timeBetween)
allEnv <- sqrt(scenarios$varEnv[r]) * ((allEnv-mean(allEnv)) / (sd(allEnv)))


## Create variables to store output
microb <- array(NA,dim=c(nmicrobe,timeWithin,timeBetween,nhost))
# random microbial compositions
microb[,1,1,] <- sample(1:smicrobe,size=nmicrobe*nhost,
                         prob=rep(1,smicrobe),
                         replace=TRUE)

# start with every host
host <- array(NA,dim=c(nhost,timeBetween))
host[,1] <- 1:nhost
pheno <- array(NA,dim=c(nhost,timeBetween))
fitness <- varpheno <- rep(NA,timeBetween)

for (j in 1:(timeBetween-1)) {
  for (h in 1:nhost) { # for all hosts
    # dynamics within host
    if (timeWithin > 1) {
    for (i in 1:(timeWithin-1)) {
        microb[,i+1,j,h] <- runWithin(microbStart=microb[,i,j,h],
                                      colonize=colonize,
                                      growth=growth,
                                      nmicrobe=nmicrobe,
                                      deathrate=deathrate,
                                      smicrobe=smicrobe,
                                      migr=migr)
      }
    }

  }

  # dynamics between host generations
  # (selection and vertical transmission)
  simm <- runBetween(microbStart=microb[,timeWithin,j,],
                    hostStart=host[,j],
                    trans=trans,
                    colonize=colonize,
                    Gmicro=Gmicro,
                    Ghost=Ghost,
                    env=allEnv[j],
                    nhost=nhost,
                    nmicrobe=nmicrobe,
                    smicrobe=smicrobe,
                    varE=varE,
                    selcoef=selcoef,
                    sex=sex,
                    correctP=correctP,
                    correctVar=correctVar)

  microb[,1,j+1,] <- simm$microbEnd
  host[,j+1] <- simm$hostEnd

  fitness[j] <- simm$fitness # fitness for this timestep
  pheno[,j] <- simm$phenotypes # phenotypes
  varpheno[j] <- var(simm$phenotypes) # variance in phenotypes for this timestep

  cat('\r \t', round((j/(timeBetween-1))*100),'%','\t of round',r,'/',nrow(scenarios))
}
