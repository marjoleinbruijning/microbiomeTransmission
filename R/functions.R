
####################################################################################
########## Code for manuscript: Natural selection for imprecise vertical 
###################### transmission in host-microbiota systems  ####################
####################################################################################
############################## Marjolein Bruijning #################################
####################################################################################
################################### 2021-09-15 #####################################
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
                        correctP) {

  # store output
  hostEnd <- rep(NA,nhost)
  microbEnd <- matrix(NA,nrow=nmicrobe,ncol=nhost)
  P <- rep(NA,nhost)

  ## Obtain host phenotypes
  for (h in 1:nhost) {
    P[h] <- Ghost[hostStart[h]] +
            sum(Gmicro[microbStart[,h]]) +
            rnorm(1,0,sqrt(varE))
  }

  ## set mean phenotype to 0?
  if (correctP) {P <- P - mean(P)}

  # Phenotypic selection based on environment
  psurv <- calcfitness(x=P,opt=env,w=selcoef)
  fitness <- mean(psurv) # average host fitness for this timestep

  # which hosts to keep, keep total host population size equal
  inc <- sample(1:nhost,nhost,replace=TRUE,prob=psurv)

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
              phenotypes=P,
              repr=inc))
}

## Function calcfitness
# Calculate fitness given trait, environment and selection
calcfitness <- function(x,opt=1,w=1) {
  exp(-(x-opt)^2/(2*w))
}
