
####################################################################################
########## Code for manuscript: When the microbiome defines the host ###############
######### phenotype: selection on transmission in varying environments #############
####################################################################################
############################## Marjolein Bruijning #################################
####################################################################################
################################### 2020-05-27 #####################################
####################################################################################

############################# Functions to be loaded ###############################


## Function runWithin
# Returns vector with present microbe species 1 time step later
runWithin <- function (microbStart,colonize,growth,
                       nmicrobe,smicrobe,deathrate,migr) {

  microbEnd <- rep(NA,length(microbStart))

  # mortality (when deathrate=1, every timestep is 1 microbe generation)
  disappear <- rbinom(nmicrobe,1,deathrate) == 1

  # individuals replaced by immigrants entering?
  imm <- rbinom(sum(disappear),1,migr)

  # replace disappearing individuals with individual from source
  microbEnd[disappear][imm == 1] <- sample(1:smicrobe,sum(imm == 1),
                                                prob=colonize,replace=TRUE)

  # replace disappearing individuals with individual from community
  # (reproduction)
  probs <- growth[microbStart]
  microbEnd[disappear][imm == 0] <- sample(microbStart,
                                           sum(imm == 0),
                                           prob=probs,replace=TRUE)

  microbEnd[!disappear] <- microbStart[!disappear]

  return(microbEnd)

}


## Function runBetween
# Returns list with microbes, hosts, fitness and phenotypes
runBetween <- function (hostStart,
                        microbStart,
                        trans,
                        colonize,
                        Gmicro,Ghost,
                        nhost,
                        nmicrobe,smicrobe,
                        varE,
                        env,
                        selcoef,
                        sex,
                        correctP,
                        correctVar) {

  # store output
  hostEnd <- rep(NA,nhost)
  microbEnd <- matrix(NA,nrow=nmicrobe,ncol=nhost)
  P <- rep(NA,nhost)

  ## Obtain host phenotypes
  if (!sex) {
    for (h in 1:nhost) {
      P[h] <- Ghost[hostStart[h]] +
              sum(Gmicro[microbStart[,h]]) +
              rnorm(1,0,sqrt(varE))
    }
  } else if (sex) {
    for (h in 1:nhost) {
      P[h] <- 0.5*Ghost[hostStart[h]] + 0.5*rnorm(1,0,sqrt(1)) +
                    sum(Gmicro[microbStart[,h]]) +
                    rnorm(1,0,sqrt(varE))
    }
  }

   ## set mean phenotype to 0
   if (correctP) {P <- P - mean(P)}

   ## set phenotypic variance to 1
   if (correctVar) {P <- rnorm(nhost,mean(P),1)}

  # Phenotypic selection based on environment
  psurv <- calcfitness(x=P,opt=env,w=selcoef)
  fitness <- mean(psurv)

  # which hosts to keep, keep total host population size equal
  inc <- sample(1:nhost,nhost,replace=TRUE,prob=psurv)

  Ne <- length(unique(inc))

  for (h in 1:length(inc)) {

    hostEnd[h] <- hostStart[inc[h]]

    ## Transmission of microbiome
    # Faithfully transmitted?
    faithTrans <- rbinom(nmicrobe,1,
                         trans[ hostEnd[h], microbStart[,inc[h]] ]
                         ) == 1

    microbEnd[faithTrans,h] <- microbStart[faithTrans,inc[h]]
    microbEnd[!faithTrans,h] <- sample(1:smicrobe,sum(!faithTrans),
                                       prob=colonize,replace=TRUE)
  }

  return(list(microbEnd=microbEnd,hostEnd=hostEnd,fitness=fitness,
              phenotypes=P,Ne=Ne,
              repr=inc))
}

## Function calcfitness
# Calculate fitness given trait, environment and selection
calcfitness <- function(x,opt=1,w=1) {
  exp(-(x-opt)^2/(2*w))
}

