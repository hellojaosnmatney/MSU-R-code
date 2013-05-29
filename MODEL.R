rm(list=ls())
options(width = 200)
setwd('/home/jason/Documents/PEF_Paper')
library(spBayes)
library(geoR)
library(gstat)
library(MBA)

## Holdout Data #########################
y.hold <- read.table("Model_Data/y.hold")
x.hold <- read.table("Model_Data/x.hold")
coords.hold <- as.matrix(read.table("Model_Data/coords.hold"))

## Model Data ###########################
y.mod <- read.table("Model_Data/y.model")
x.mod <- read.table("Model_Data/x.model")
coords.mod <- as.matrix(read.table("Model_Data/coords.model" ))

#########################################
## TrueHt
#########################################

X <-  cbind(x.mod[,2],x.mod[,3]) ## 50 & 95
coords <- coords.mod
mod <- lm(y.mod[,1]~X)
summary(mod)
vario <- variog(coords = coords,
                data = resid(mod),
                uvec=(seq(0, 500, length=10)))

variomod <- variofit(vario, c(0.04, 150))   ## partial sill and range
plot(vario); lines(variomod, col = "red")

tau <- variomod$nugget
sigma <- variomod$cov.pars[1]
phi <- 3/variomod$cov.pars[2]

#########################################
## CLAVG
#########################################

X <- cbind(x.mod[,3])  ##  95
mod <- lm(y.mod[,2]~X)
summary(mod)
vario <- variog(coords = coords, data = resid(mod),  uvec=(seq(0, 500, length=10)))
variomod <- variofit(vario, c(0.05, 100))   ## partial sill and range
plot(vario);lines(variomod, col = "red")

#########################################
## CRAVG
#########################################
X <-  cbind(x.mod[,3]) ## 95
mod <- lm(y.mod[,3]~X)
summary(mod)
vario <- variog(coords = coords, data = resid(mod),  uvec=(seq(0, 500, length=10)))
variomod <- variofit(vario, c(0.1, 200))   ## partial sill and range
plot(vario);lines(variomod, col = "red")

#########################################
## dbh
#########################################
X <-  cbind(x.mod[,2],x.mod[,3])  ## 50 & 95
mod <- lm(y.mod[,4]~X)
summary(mod)
vario <- variog(coords = coords, data = resid(mod),  uvec=(seq(0, 500, length=10)))
variomod <- variofit(vario, c(0.1, 300))   ##partial sill and range
plot(vario)
lines(variomod, col = "red")

########################################
## FIT THE SPATIAL MODEL
########################################

set.seed(1)
spmod <- spLM(y.mod[,1] ~X , coords=coords,
              starting=list("phi"=phi,"sigma.sq"=sigma, "tau.sq"=tau),
              sp.tuning=list("phi"=0.5, "sigma.sq"=0.03, "tau.sq"=0.03),
              priors=list("phi.Unif"=c(3/500, 3/10), "sigma.sq.IG"=c(2, sigma),
                "tau.sq.IG"=c(2, tau)),cov.model="exponential",
              n.samples=1000, verbose=TRUE, n.report=10)

plot(spmod$p.samples,density=F)

##Saved stuff that Chad wrote######################

samps  <- read.table('/home/jason/Documents/PEF-LiDAR/Samples') ## Samples = samps + colnames
w.samps <- read.table('/home/jason/Documents/PEF-LiDAR/w.samps')

## investigate spmod #######################

plot(spmod$p.samples) 

## what are the beta values? ###############
## ie intercept and slope parameters. ######

B <- apply(spmod$p.samples,2,mean)

## the burn-in period of 2000 samples ######
samps <- spmod$p.samples[2001:10000,]
B <- apply(samps,2,mean)

## this fits 1000 different values for each beta
## (we take the means of those 1000 values as our estimate of beta)
## we also get estimates for our nugget, psill, and range here.

## what does w look like?
## For each location 'spLM' fit 1000 values.
## We take the mean of the 1000 values at each
## location to be the estimate of w(s).

w <- rowMeans(spmod$sp.effects)
## burn in w #################################
w.samps <- spmod$sp.effects[,2001:10000]
w <- rowMeans(w.samps)


## write out samps and w.samps so we dont have to run spmod again.
write.table(samps,"samps",row.names = FALSE, col.names = FALSE)
write.table(w.samps, "w.samps", row.names = FALSE, col.names = FALSE)

## what does our new 'e' look like? ##########
## spmod doesnt return this directly #########
## so we have to calculate it. ###############

## y = xb+w+e so, y-xb-w = e. ################
b <- matrix(B[1:2], ncol = 1)
resid.spmod <- Y - X%*%b - w

## better way to get e
y.hat.samps <- matrix(nrow = 494, ncol = 8000)
for(j in 1:494){
  for(i in 1:8000){
    y.hat.samps[j,i] <- rnorm(1,samps[i,1] + samps[i,2]*X[j,2] + w.samps[j,i], samps[j,4])
  }
}

y.hat <- rowMeans(y.hat.samps)


SSy <- sum((Y-mean(Y))^2)
SSe <- sum((Y - (X%*%b + w))^2)
SSr <- SSy - SSe

SSr/SSy ## R^2 is one.

Y.hat = X%*%b + w

plot(Y,Y.hat)

## we can be fancier about this ##############
## (by factoring in all the samples) #########
## but this will get the point across. #######

newvario <- variog(coords = coords, data = resid.spmod, max.dist=250)
plot(newvario)

## notice that the variogram psill
## looks like it should be zero.
## (this means there is no spatial dependence in our residuals now).
## this is what we wanted.
##

vario  <- variog(coords = coords, data = resmod, max.dist = 250)
j2 <- fit.variogram(j, )
plot(j)
lines(j2)

v_model <- vgm(mean(samples$sigma.sq), "Exp", mean(samples$phi), mean(samples$tau.sq))

v_model
______________________________
  model        psill     range
1   Nug 9.029489e-05 0.0000000
2   Exp 2.806277e+01 0.1333455



                    
## > mean(samples$sigma.sq)
## [1] 28.06277
## > mean(samples$phi)
## [1] (3/0.1333455) = 22.49963
## > mean(samples$tau.sq)
## [1] 9.029489e-05

vario  <- variog(coords = coords, data = resmod, max.dist = 250)
fit.vario <- vgm(28.06277, "Exp", 3/0.1333455, nugget = 9.029489e-05)

plot(vario, ylim = c(0,50))

lines.variomodel(cov.model = "exp", cov.pars = c(28.06277, 3/0.1333455), nug = 9.029489e-05)

abline(h = mean(samps[,4]), col = "blue")
abline(h = mean(samps[,3]) + mean(samps[,4]), col = "green")
abline(v = 3/mean(samps[,5]), col = "red3")

vario.DBH.resid <- variog
                   (coords=coords,
                    data=DBH.resid,
                    uvec=(seq(0, max.dist, length=bins)))

fit.DBH.resid <-variofit(vario.DBH.resid, ini.cov.pars=c(300, 200/-log(0.05)),
                         cov.model="exponential", minimisation.function="nls", weights="equal")
