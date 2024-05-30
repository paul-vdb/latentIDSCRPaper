#################################
## Latent ID SCR
## Fisher Example
#################################
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(ggplot2)

source("NimbleFunctions.R")

## Load the Fisher data
###--- set directory and load files
load("../data/FisherData.Rda")

buffer <- 6
traps <- fisher.data[["traps"]]

xlim <- range(traps[,1])+c(-buffer,buffer)
ylim <- range(traps[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100

obs <- fisher.data[["observations"]]
collars <- fisher.data$collarInfo

## Start: January 20, 2016
start.date <- as.Date("01-01-2016", format = "%d-%m-%Y")
end.date <- as.Date("05-03-2016", format = "%d-%m-%Y")
StudyPeriod <- as.numeric(end.date - start.date)

## 207 detections.
obs <- fisher.data[["observations"]]
obs <- obs[as.Date(obs$DateTime) >= start.date & end.date > as.Date(obs$DateTime),]
nrow(obs)

omega <- obs$TrapNumber ## Trap Mark
M <- 400		# Data augmentation super pop
J <- nrow(traps)

## Chandler and Royle Spatial Count model, Algorithm 2
## Same as MPP marginalized over all ID. 
counts <- numeric(J)
for(j in 1:J)
{
	counts[j] <- sum(obs$TrapNumber == j)
}

Model1 <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 50)
	lambda ~ dunif(0, 20)
	psi ~ dbeta(1,1)
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau2 )*z[k]
	}

	### Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time
		counts[j] ~ dpois(Hj[j])
	}
	
	N <- sum(z[1:M])
	D <- N/area
} )

# Initialize model:
#-------------------------
init1 <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	psi <- rbeta(1,1,1)
	z <- rbinom(M, prob=psi, size = 1)
	z[1:nrow(traps)] <- 1
	X[1:nrow(traps),] <- as.matrix(traps) + cbind(rnorm(nrow(traps), 0, 0.1), rnorm(nrow(traps), 0, 0.1))
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z
    )
}

constants1 <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
	area = diff(xlim)*diff(ylim)/100
	)

data1 <- list(
	z =  rep(NA, M),
	counts = counts
	)

Rmodel <- nimbleModel(Model1, constants1, data1, inits = init1())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M)	{
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE)	# Adaptive is okay for this model.
}
conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out1 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init1(), init1(), init1()), samplesAsCodaMCMC = TRUE)

# plot(out1[,c("sigma", "lambda", "psi")])
# plot(out1[,c("sigma", "N", "D")])
# summary(out1)
# effectiveSize(out1)

## Marked Poisson process model incorporating the 
## latent ID variable as unknown. 
##------------------------------------------------
ModelMPP <- nimbleCode({
    lambda ~ dunif(0, 20)   # Detection rate at distance 0
    psi ~ dbeta(1, 1)       # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...

    for(k in 1:M) {
        ## Standard SCR data augmentation, animal in pop, loc and dist to traps.
		z[k] ~ dbern(psi)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2

		## Now the rate for the Poisson process:
		## Don't put lambda here for efficiency of updates.
        pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)
        ## Hazard rate for animal across all traps.
		Hk[k] <- sum(pkj[k,1:J])*Time*lambda*z[k]
		zeros[k] ~ dpois(Hk[k])
   }

	pID[1:M] <- z[1:M]*lambda	# Puts a prior on ID of z==1, adds lambda to the detections likelihood for efficiency.

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Trap probability given ID: Returns pkj[ID_i, omega_i]
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])
		# ID probability, returns pID[ID_i], does not check for proper probabilities.
		ID[i] ~ dID(pID = pID[1:M])	
    }
	
    # Estimated population size
    N <- sum(z[1:M])
	D <- N/area
})

# Run the same model from van Dam-Bates et al. Marked Poisson Process on fisher.
#----------------------------------------------------------------------
constants.mpp <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100
	)

data.mpp <- list(
    zeros = rep(0, M),
	z =  rep(NA, M),
	ID = rep(NA, length(omega)),
	omega = omega	
)

# Need to initialize this model as the stochastic node for ID is kind of wrong...
inits.mpp <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 2)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	psi <- rbeta(1,1,1)
	z <- rbinom(M, size = 1, prob = psi)
	ID <- do.call('c', lapply(omega, FUN = function(x) {sample(1:M, 1, prob = z*hkj[,x])}))
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID
    )
}

## Latent ID SCR as MPP
######################################
Rmodel <- nimbleModel(ModelMPP, constants.mpp, data.mpp, inits = inits.mpp())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))

# Use a block update on locations. Saves time.
# Turn off adaptive sampling and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, 
		control = list(scale = 0.5, xlim = xlim, ylim = ylim, temp = 0.2))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))	

## Optimized z sampler, it is a small change from standard as it does a check if ID is matched.
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', targetByNode = TRUE)
## van Dam-Bates categorical sampler, although Nimble has updated their categorical function so they are very similar.
conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', targetByNode = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out2 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(inits.mpp(), inits.mpp(), inits.mpp()), samplesAsCodaMCMC = TRUE)

# plot(out2[, c('sigma', 'D')])
# plot(out2[, c('lambda', 'N')])
# summary(out2[, c('lambda, 'sigma', 'N', 'D')])

Cmcmc$run(niter = 10000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])


####################################
# Model 3: Sex as a covariate
####################################
Model3 <- nimbleCode({
    lambda ~ dunif(0, 20)
    sigma ~ dunif(0, 50)
    tau2 <- 1/(2*sigma^2)

    psi ~ dbeta(1, 1)      	# prior on data augmentation Bernoulli vec.

	psex ~ dbeta(1, 1)		# Sex Ratio 

	for(k in 1:M) {
		z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)

		# Hazard rate for animal across all traps.
		Hk[k] <- sum(pkj[k,1:J])*Time*z[k]*lambda
		zeros[k] ~ dpois(Hk[k])
	}

	pID[1:M] <- z[1:M]*lambda

    ## Trap history model and unobserved animal ID.
    for(i in 1:n_obs) {
		# Hard match/no match info.
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		keep[i] ~ dbern(pSex[i])
		
		# Trap Prob.
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

		# ID prior based on z, with lambda mixed in.
		ID[i] ~ dID(pID = pID[1:M])
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

obsSex <- obs$sex
constants3 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	obsSex = obsSex
)

data3 <- list(
    zeros = rep(0, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  rep(NA, M),
	ID = rep(NA, length(omega)),
	sex = rep(NA, M)	# Note 0 is male and 1 is female.
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init3 <- function(){
	N <- floor(runif(1, 1, M/2))
	psi <- rbeta(1, N, M-N)	# NA inits...	
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, 1, 1)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sex <- rbinom(M, size = 1, prob = psex)
	ID <- do.call('c',lapply(1:length(omega), FUN = function(x){sample(1:M, 1, prob = hkj[1:M,omega[x]]*(sex+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	Hk <- rowSums(hkj)*StudyPeriod*lambda
	p <- exp(-Hk)*psi/(exp(-Hk)*psi + (1-psi))
	z <- rbinom(M, size = 1, prob=p)
	z[ID] <- 1
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID,
		sex = sex,
		psex = psex
    )
}

###################################
## Run the Sex fisher model MCMC
###################################
Rmodel <- nimbleModel(Model3, constants3, data3, inits = init3())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'z', 'sex', 'psex')) # 'psex', 'sex',

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, 
		control = list(scale = 0.5, xlim = xlim, ylim = ylim, temp = 0.2))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))	

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', targetByNode = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', targetByNode = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# The full run...
out3 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init3(), init3(), init3()), samplesAsCodaMCMC = TRUE)

# plot(out3[, c('sigma', 'lambda', 'psex')])
# plot(out3[, c('D', 'N')])
# summary(out3[, c('lambda, 'sigma', 'psex', 'N', 'D')])

####################################
## Model 4: Sex and collar as covariates
####################################
Model4 <- nimbleCode({
	
    lambda ~ dunif(0, 20)

	# Home range parameter m=0/f=1:
    sigma ~ dunif(0, 50)	

	# convert now to make it simpler in half-normal computation
    tau2 <- 1/(2*sigma^2)	

    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.

	psex ~ dbeta(1, 1)

	# For the collared individuals, we observe the sex too.
	# So we have observed sex 1:14 and z = 1 for those animals.
	for(k in 1:M) {
		z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
		
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2

		pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*lambda	## In this case I've put lambda here... this was for sex specific rates if necessary.
		# Hazard rate for animal across all traps.
		Hk[k] <- sum(pkj[k,1:J])*Time*z[k]

		# Zero count for all.
		zeros[k] ~ dpois(Hk[k])
	}

	# Prior on animal ID is just 1|z = 1 or 0|z = 0. 
	# The reason is in homogeneous time things cancel nicely and this is easiest.
	pID[1:M] <- z[1:M]

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
		# Hard match/no match info.
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		# We add matchCollar here to make sure collar matches.
		pMatch[i] <- matchCollar[ID[i],i]*pSex[i]
		# Hard 1/0 for animal match to sex and collar.
		keep[i] ~ dbern(pMatch[i])
		
		# Trap Prob:
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

		# ID prior based on z:
		ID[i] ~ dID(pID = pID[1:M])
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

# Process collar information per animal.
########################################
collars <- fisher.data$collarInfo
matchCollar <- matrix(1, nrow = M, ncol = length(omega))
collars$Sex <- 2 - (collars$Sex == "M")
dates <- obs$DateTime
dates <- as.Date(dates)
obsCollar <- obs$collar

# Process collar matrix:
for(i in 1:length(omega))
{
	indx.match <- collars$Start < dates[i] & collars$Finish >= dates[i]
	if(obsCollar[i] == 1)
	{
		matchCollar[(nrow(collars)+1):M,i] <- 0
		matchCollar[1:nrow(collars), i] <- indx.match*1
	}else{
		if(obsCollar[i] == 2) matchCollar[1:nrow(collars),i] <- (1-indx.match)		
	}
	if(dates[i] > as.Date("2016-01-25", format = "%Y-%m-%d")){
		matchCollar[which(collars$FisherID == "M02"),i] <- 0
	}
}

constants4 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	obsSex = obs$sex,
	obsMark = obs$collar
	)

data4 <- list(
    zeros = rep(0, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  c(rep(1, nrow(collars)), rep(NA, M-nrow(collars))),
	ID = rep(NA, length(omega)),
	sex = c(collars$Sex - 1, rep(NA, M - nrow(collars))),	# Note 0 is male and 1 is female.
	matchCollar = matchCollar	# needs to be data as we can't legally dynamically index constants...
)

init4 <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, 9, 14)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sexCollar <- collars$Sex-1
	sex <- c(sexCollar, rbinom(M- nrow(collars), size = 1, prob = psex))
	ID <- do.call('c',lapply(1:length(omega), 
		FUN = function(x){sample(1:M, 1, prob = hkj[,omega[x]]*matchCollar[,x]*(obs$sex[x] == 0 | (obs$sex[x] - 1) == sex))}))
	sex[1:nrow(collars)] <- NA
	z <- rep(0, M)
	z[ID] <- 1
	z[1:nrow(collars)] <- NA
	psi <- rbeta(1, sum(z, na.rm = TRUE) + nrow(collars), M - sum(1-z, na.rm = TRUE))	# NA inits...
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID,
		sex = sex,
		psex = psex
    )
}

###################################
# Chandler and Royle Sampler:
###################################
Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex'))

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, control = list(scale = 1, xlim = xlim, ylim = ylim, temp = 0.2))
}

conf$removeSamplers('sigma')
# sigma slice sampling.
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, targetByNode = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', targetByNode = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', targetByNode = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# The full run...
out4 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init4(), init4(), init4()), samplesAsCodaMCMC = TRUE)

plot(out4[, c("N", "D")])
plot(out4[, c("sigma", "lambda", "psex")])


