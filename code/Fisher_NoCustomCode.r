#################################
## Latent ID SCR
## Fisher Example
## Marked Poisson Process SCR
## but no custom distributions or samplers.
#################################
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(ggplot2)

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
        Hkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*Time*lambda
        ## Hazard rate for animal across all traps.
		Hk[k] <- sum(Hkj[k,1:J])*z[k]
		zeros[k] ~ dpois(Hk[k])
   }

	H. <- sum(Hk[1:M])
	pID[1:M] <- Hk[1:M]/H.	# Puts a prior on ID of z==1.

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Trap probability given ID:
		pobs[i, 1:J] <- Hkj[ID[i], 1:J]/Hk[ID[i]]
		## Minimum speed up to make this a ones trick instead as this is just indexing the probability
        omega[i] ~ dcat(pobs[i,1:J])	
		# ID probability:
		ID[i] ~ dcat(pID[1:M])
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
		type = 'RW_block', silent = TRUE, 
		control = list(adaptive = FALSE, scale = 0.75, xlim = xlim, ylim = ylim))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))	

## Optimized z sampler, it is a small change from standard as it does a check if ID is matched.
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## Painfully Slow. That' why things were customized.
Cmcmc$run(niter = 10000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out.mpp <- mcmc(samples[-(1:5000),])

# outmpp <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	# thin = 1, inits = list(inits.mpp(), inits.mpp(), inits.mpp()), samplesAsCodaMCMC = TRUE)

# plot(outmpp[, c('sigma', 'D')])
# plot(outmpp[, c('lambda', 'N')])
# summary(outmpp[, c('lambda, 'sigma', 'N', 'D')])
