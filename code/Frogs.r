#################################
# Latent ID ASCR
# Frog Example with Time of Arrival.
#################################
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)
library(ggplot2)
library(secr)
library(ascr)

source("NimbleFunctions.R")

load("../data/stacked-lightfooti.Rdata")

# Change the mask to finer resolution for the CD MLE.
traps2 <- convert.traps(traps)
mask <- make.mask(traps = traps2, buffer = 15, spacing = 0.2, type = "trapbuffer")
A <- attr(mask, "area")
mask <- as.matrix(mask)
attr(mask, "area") <- A
attr(mask, "buffer") <- 15

###########
# Latent ID ASCR from Stevenson 2021
###########
inits <- function(){
	ID <- numeric(length(occ))
	ID[occ == 1] <- 1:sum(occ == 1)
	ID[occ == 2] <- 1:sum(occ == 2)
	z <- matrix(0, ncol = 2, nrow = M)
	z[ID[occ==1],1] <- 1
	z[ID[occ==2],2] <- 1
    p <- runif(1, 0.1, 0.3)
	lambda = runif(1, 0.1, 2)
	sigma = runif(1, 3, 5)
	g0 = runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
	X <- array(c(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]), 
		runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])), c(M,2,2))		
	
	for(k in 1:2){
		for(i in 1:M){
			if(sum(ID[occ == k] == i) == 0) next;
			if(sum(ID[occ == k] == i) == 1) pID <- panimal[which(ID[occ == k] == i), ]
			if(sum(ID[occ == k] == i) > 1) pID <- colSums(panimal[which(ID[occ == k] == i), ])
			mpt <- sample(ncol(panimal), 1, prob = exp(pID))
			X[i,,k] <- mask[mpt,]
		}
	}

	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = runif(1, 0.01, 1),
		g0 = g0,
		X = X,
        ID = ID,
		z = z
    )
}

code <- nimbleCode({
	lambda ~ dunif(0, 10) # Detection rate at distance 0
	psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
	sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
	tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0, 1)
	g0 ~ dunif(0, 20)
	for(v in 1:n_occ) {
		for(k in 1:M) {
			z[k,v] ~ dbern(psi)
			X[k, 1, v] ~ dunif(xlim[1], xlim[2])
			X[k, 2, v] ~ dunif(ylim[1], ylim[2])
			d2[k,1:J, v] <- (X[k,1,v]-traps[1:J,1])^2 + (X[k,2,v]-traps[1:J,2])^2
			expTime[k, 1:J, v] <- sqrt(d2[k,1:J,v])/nu
			pkj[k,1:J,v] <- (1-exp(-g0*exp(-d2[k,1:J,v]*tau2)))
			# Hazard rate for animal across all traps.
			Hk[k,v] <-(1-prod(1-pkj[k,1:J,v]))*lambda*Time*z[k,v]
			zeros[k,v] ~ dpois(Hk[k,v])
		}
		# Predicted population size for each occ
		N[v] <- sum(z[1:M,v])
	}
	
	pID[1:M, 1:n_occ] <- z[1:M,1:n_occ]*lambda
	
	# Trap history model.
	# and unobserved animal ID.
	for(i in 1:n_obs) {
		# Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J, occ[i]])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i], 1:J, occ[i]], sd = sigmatoa, y = y[i,1:J])
		# The likelihood needs to be multiplied by lambda for each detection and
		# I need ID to be a stochastic node. 2 birds...
		ID[i] ~ dID(pID = pID[1:M, occ[i]])
	}

	# Derived Variables.
	EN <- psi*M
	D <- EN/area*10000
})

xlim <- range(mask[,1])
ylim <- range(mask[,2])
area <- diff(xlim)*diff(ylim)
toa <- capt.all$toa
capt <- capt.all$bincapt[,1:6]
tmin <- apply(toa, 1, max)
occ <- 1+(tmin > 1200)
ID <- capt.all$bincapt[, 7]
ID <- as.integer(as.factor(ID))
IDBen <- ID

# Constants:
M <- 200
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30

constants <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
    n_obs = nrow(capt),
	trials = rep(1, J),
	nu = nu,
	area = area,
	n_occ = 2,
	occ = occ,
	mi = rowSums(capt))

data <- list(
	zeros = cbind(rep(0,M), rep(0,M)),
    y = capt,
	toa = toa,
	z = cbind(rep(NA, M), rep(NA, M)),
	ID = rep(NA, nrow(capt))
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D', 'ID'))

conf$removeSamplers('X')
for(v in 1:2){
	for(i in 1:M) {
		conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), 
			type = 'myJAM', silent = TRUE, 
			control = list(scale = 1, xlim = xlim, ylim = ylim, temp = 0.2, occ = 2))
	}
}

conf$removeSamplers('g0')
conf$addSampler(target = 'g0', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', targetByNode = TRUE, 
	control = list('Noccasion' = 2, 'IDoccasion' = occ))

conf$removeSamplers('ID')
# Sampler from Chandler and Royle 2013
conf$addSampler('ID', type = 'myCategorical', targetByNode = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out <- runMCMC(Cmcmc, nburnin = 20000, niter = 50000, nchains = 3, 
	inits = list(inits(), inits(), inits()), samplesAsCodaMCMC = TRUE)	

summary(out)

###########
# Known ID ASCR from Stevenson 2021
###########
initsID <- function(){
    p <- runif(1, 0.1, 0.3)
	z <- matrix(rbinom(M*2, 1, p), ncol = 2, nrow = M)
	z[IDBen[occ==1],1] <- NA
	z[IDBen[occ==2],2] <- NA
	lambda = runif(1, 0.1, 2)
	sigma = runif(1, 3, 5)
	g0 = runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
	X <- array(c(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]), 
		runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])), c(M,2,2))		
	
	for(k in 1:2){
		for(i in 1:M){
			if(sum(IDBen[occ == k] == i) == 0) next;
			if(sum(IDBen[occ == k] == i) == 1) pID <- panimal[which(IDBen[occ == k] == i), ]
			if(sum(IDBen[occ == k] == i) > 1) pID <- colSums(panimal[which(IDBen[occ == k] == i), ])
			mpt <- sample(ncol(panimal), 1, prob = exp(pID))
			X[i,,k] <- mask[mpt,]
		}
	}

	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = runif(1, 0.01, 1),
		g0 = g0,
		X = X,
		z = z
    )
}

ID <- capt.all$bincapt[, 7]
IDBen[occ == 1] <- as.integer(as.factor(ID[occ == 1]))
IDBen[occ == 2] <- as.integer(as.factor(ID[occ == 2]))
K <- c(max(IDBen[occ == 1]), max(IDBen[occ == 2]))

constants.id <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = Time,
    M = M,
    n_obs = nrow(capt),
	trials = rep(1, J),
	nu = nu,
	area = area,
	n_occ = 2,
	occ = occ
)

data.id <- list(
	zeros = matrix(0, nrow = M, ncol = 2),
    y = capt,
	toa = toa,
	z = cbind(c(rep(1, K[1]), rep(NA, M-K[1])),
		c(rep(1, K[2]), rep(NA, M-K[2]))),
	ID = IDBen
)

Rmodel <- nimbleModel(code, constants.id, data.id, inits = initsID())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D'))

conf$removeSamplers('X')
for(v in 1:2){
	for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), type = 'RW_block', silent = TRUE)
}

conf$removeSamplers('g0')
conf$addSampler(target = 'g0', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

out.id <- runMCMC(Cmcmc, nburnin = 20000, niter = 50000, nchains = 3, 
	inits = list(initsID(), initsID(), initsID()), samplesAsCodaMCMC = TRUE)	

summary(out.id)

############################
# Now run it with ASCR, 
# but we need to increase the mask:
############################
multi.capt <- list(list("bincapt" = capt[occ==1,], "toa" = toa[occ==1,]), list("bincapt" = capt[occ ==2,], "toa" = toa[occ==2,]))
ascr.res <- fit.ascr(capt = multi.capt, list(traps, traps), mask = list(mask, mask), 
			detfn = "hhn", survey.length = c(30, 30))
# citation("ascr")
summary(ascr.res)
confint(ascr.res)
# Detection function: Hazard halfnormal 
# Information types: Times of arrival
 
 # Parameters: 
            # Estimate Std. Error
# D         1.2489e+02    10.1940
# lambda0   6.5932e+00     1.0356
# sigma     2.2916e+00     0.0943
# sigma.toa 9.6946e-04     0.0001
# ---                            
                               
# esa.1     2.4556e-02     0.0009
# esa.2     2.4556e-02     0.0009
# > confint(ascr.res)
                 # 2.5 %       97.5 %
# D         1.049065e+02 1.448663e+02
# lambda0   4.563476e+00 8.622954e+00
# sigma     2.106834e+00 2.476330e+00
# sigma.toa 7.833584e-04 1.155556e-03
