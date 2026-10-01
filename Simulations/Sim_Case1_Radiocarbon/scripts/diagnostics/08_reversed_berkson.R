# Reversing Case 1: the calendar date as the response, not the predictor
#
# Case 1 puts the date on the x axis, so it is a predictor measured with error
# (classical errors in variables) and the slope attenuates. Here the date is the
# response of a covariate d known exactly, following measurement_error_example.R
# in ercrema/statistical_modelling_review. Error in the response should leave
# the slope alone, in theory, but only if it does not depend on the predictor.
# Radiocarbon error does, and the slope halves: see
# Notes/reversed_berkson_2026-09-22.md.
#
# Follows Crema's example closely, same variable names and order. Additions:
# lab error on top of the curve error, a slopes/sigmas table, and the agreement
# index between each calibrated date and its posterior. No Stan version: it did
# not fit correctly and Crema's own code is Nimble.
#
#   Rscript Simulations/Sim_Case1_Radiocarbon/scripts/diagnostics/08_reversed_berkson.R

# Load Relevant R Packages ----
library(rcarbon)
library(nimbleCarbon)
library(coda)
library(here)

out_dir <- here('Simulations', 'Sim_Case1_Radiocarbon', 'output', 'diagnostics', 'reversed')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Simulate Regression Data ----
set.seed(12332) #Random seed
true.alpha  <- 2700 #True Intercept
true.beta  <- 1/3.7 #True Slope
true.sigma <- 70 #True Error
n <- 100 #Sample Size
d <- runif(n=n,min=0,max=1000) |> round()
calendar.dates <- rnorm(n=n,mean=true.alpha - true.beta*d,sd=true.sigma) |> round()
c14ages.error <- rep(30,n) #assign errors
# Crema back-calibrates with uncalibrate()[,4], which carries the curve error
# only. Case 1 adds the lab error on top, because calibrate() assumes both.
curve <- uncalibrate(calendar.dates,verbose=FALSE)
c14ages <- round(rnorm(n,curve$ccCRA,sqrt(c14ages.error^2+curve$ccError^2)))
median.calibrated <- calibrate(c14ages,c14ages.error,verbose=FALSE) |> medCal() #calibrate and compute median calibrated dates

# Run regression analysis on calendar dates ----
fitcal  <- lm(calendar.dates~d,data=data.frame(d=d,calendar.dates=calendar.dates))

# Run regression analyses on median calibrated dates ----
fitmed <- lm(median.calibrated~d,data=data.frame(d=d,median.calibrated=median.calibrated))
confint(fitmed)

# Run Bayesian EIV via Nimble ----
# Prepare data list:
dat <- list()
dat$cra <- c14ages
dat$cra.error <- c14ages.error
# Prepare constant list:
constants <- list()
constants$d <- d
data(intcal20)
constants$calBP <- intcal20$CalBP
constants$C14BP  <- intcal20$C14Age
constants$C14err <- intcal20$C14Age.sigma
constants$n <- n
# Prepare initialisation values
inits <- list()
inits$theta <- median.calibrated
inits$alpha <- 3000
inits$beta <- 1/2
inits$sigma <- 50

# Regression model in nimble
model <- nimbleCode({
	for (i in 1:n)
	{
		mu[i]  <- alpha +  beta * d[i]
		theta[i] ~ dnorm(mean=mu[i],sd=sigma)
		c14age[i] <- interpLin(z=theta[i],x=calBP[],y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i],x=calBP[],y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i])
	}
	alpha ~ dunif(1000,5000)
	beta ~ dnorm(0,1)
	sigma ~ dexp(0.1)
})

# Run MCMC
out  <- nimbleMCMC(code=model,constants=constants,inits=inits,data=dat,nchains=4,niter=100000,nburnin=50000,samplesAsCodaMCMC=TRUE,monitors=c('theta','alpha','beta','sigma'))

# Check Convergence in MCMC
rhats <- coda::gelman.diag(out,multivariate=FALSE)
which(rhats[[1]][,1]>1.01) #Only theta, can be ignored as calibrated 14C dates are not normally distributed

# Extract Posterior
posterior <- do.call(rbind.data.frame,out)

# Compute HPD interval
HPDinterval(as.mcmc(posterior[,'alpha']))
HPDinterval(as.mcmc(posterior[,'beta']))

# Extract Mean Posterior for Theta
theta.median <- apply(posterior[,grep('theta',colnames(posterior))],2,median)

# Collect the three slopes ----
# The model is written as alpha + beta * d with dates in cal BP, which run
# backwards, so every fitted slope is negative where true.beta is positive. The
# signs are flipped here so the rows can be read against true.beta.
beta.hpd <- HPDinterval(as.mcmc(posterior[,'beta']))
slopes <- data.frame(
	source = c('true calendar dates (lm)','median calibrated dates (lm)','EIV, Nimble'),
	beta = -c(coef(fitcal)[2],coef(fitmed)[2],mean(posterior$beta)),
	lower = -c(confint(fitcal)[2,2],confint(fitmed)[2,2],beta.hpd[2]),
	upper = -c(confint(fitcal)[2,1],confint(fitmed)[2,1],beta.hpd[1]))
slopes$ratio <- slopes$beta/true.beta

sigmas <- data.frame(
	source = slopes$source,
	sigma = c(summary(fitcal)$sigma,summary(fitmed)$sigma,mean(posterior$sigma)))

write.csv(slopes,file.path(out_dir,'slopes.csv'),row.names=FALSE)
write.csv(sigmas,file.path(out_dir,'sigmas.csv'),row.names=FALSE)
print(slopes)
print(sigmas)

# Agreement between the calibrated dates and the posterior dates ----
# agreementIndex is the index OxCal reports, where 60 is the usual threshold.
theta.posterior <- as.matrix(posterior[,grep('theta',colnames(posterior))])
agr <- agreementIndex(c14ages,c14ages.error,theta=theta.posterior,verbose=FALSE)

agreement <- data.frame(date=1:n,agreement=agr$agreement,
                        cal.median=median.calibrated,post.median=theta.median,
                        true.date=calendar.dates)
agreement$shift <- agreement$post.median - agreement$cal.median
write.csv(agreement,file.path(out_dir,'agreement.csv'),row.names=FALSE)

cat('Overall agreement:',round(agr$overall.agreement,1),'\n')
cat('Dates below 60:',sum(agr$agreement<60),'\n')

# Make a data.frame for plotting the results ----
pred <- data.frame(d=-100:1200)

# Regression on True Calendar Dates:
pred$true  <- predict(fitcal,newdata=pred)
pred$true.lo  <- predict(fitcal,newdata=pred,interval='confidence')[,2]
pred$true.hi  <- predict(fitcal,newdata=pred,interval='confidence')[,3]

# Regression on Median Calibrated Dates:
pred$m  <- predict(fitmed,newdata=pred)
pred$m.lo95  <- predict(fitmed,newdata=pred,interval='confidence')[,2]
pred$m.hi95  <- predict(fitmed,newdata=pred,interval='confidence')[,3]

# Bayesian Regression on 14C Dates:
predmatrix <- matrix(NA,nrow=nrow(posterior),ncol=length(-100:1200))
for (i in 1:nrow(posterior))
{
	predmatrix[i,]  <- posterior$alpha[i] + posterior$beta[i] * c(-100:1200)
}
pred$b <- apply(predmatrix,2,mean)
pred$b.lo95 <- apply(predmatrix,2,function(x){HPDinterval(mcmc(x))[1]})
pred$b.hi95 <- apply(predmatrix,2,function(x){HPDinterval(mcmc(x))[2]})

# Extract values for error plot ----
cl <- calibrate(c14ages,c14ages.error,calMatrix=TRUE,timeRange=c(3000,2000),verbose=FALSE)
caldd <- vector('list',length=100)
sc  <- 1000
for (i in 1:100)
{
	xx  <- c(d[i] + cl$calmatrix[,i]*sc,d[i] - rev(cl$calmatrix[,i])*sc)
	yy  <- c(3000:2000,2000:3000)
	caldd[[i]] <- data.frame(xx=xx,yy=yy)
}

# Plot Results ----
png(file.path(out_dir,'reversed_regression.png'),height=3.5,width=9,units='in',res=300)
par(mfrow=c(1,3),mar=c(4,4,2,1))
par(lend=2)
plot(d,calendar.dates,xlim=c(0,1000),ylim=c(2850,2100),pch=19,ylab='BP',xlab='x',col=adjustcolor('black',0.6))
abline(a=true.alpha,b=-true.beta,lty=2,lwd=2)
polygon(x=c(pred$d,rev(pred$d)),y=c(pred$true.hi,rev(pred$true.lo)),border=NA,col=adjustcolor('firebrick',0.3))
lines(pred$d,pred$true,lwd=2,col='firebrick')
legend(x=20,y=2100,legend=c('Calendar Date','True Relationship','Regression on Calendar Date'),pch=c(19,NA,NA),lwd=c(NA,1,8),col=c('black','black','firebrick'),bty='n',cex=0.9,lty=c(NA,2,1))

plot(NA,xlim=c(0,1000),ylim=c(2850,2100),ylab='BP',xlab='x')
for(i in 1:100)
{
	polygon(caldd[[i]]$xx,caldd[[i]]$yy,border=NA,col=adjustcolor('lightblue',1))
}

for (i in 1:length(calendar.dates))
{
	lines(x=c(d[i],d[i]),y=c(calendar.dates[i],median.calibrated[i]),lty=3,col='black')
}

points(d,calendar.dates,pch=20)
points(d,median.calibrated,pch=4)
legend(x=20,y=2100,legend=c('Calendar Date','Median Calibrated Date','Calibrated Distribution'),pch=c(19,4,NA),lwd=c(NA,NA,8),col=c('black','black','lightblue'),bty='n',cex=0.9)
plot(NA,xlim=c(0,1000),ylim=c(2850,2100),pch=19,ylab='BP',xlab='x')
abline(a=true.alpha,b=-true.beta,lty=2,lwd=2)
lines(pred$d,pred$m,lwd=2,col='darkorange')
polygon(x=c(pred$d,rev(pred$d)),y=c(pred$m.lo95,rev(pred$m.hi95)),border=NA,col=adjustcolor('darkorange',0.3))
lines(pred$d,pred$b,lwd=2,col='darkgreen')
polygon(x=c(pred$d,rev(pred$d)),y=c(pred$b.lo95,rev(pred$b.hi95)),border=NA,col=adjustcolor('darkgreen',0.3))
legend(x=20,y=2100,legend=c('True Relationship','Regression on Median Date','Bayesian EIV Model'),lwd=c(1,8,8),col=c('black','darkorange','darkgreen'),bty='n',cex=0.9,lty=c(2,1,1))

dev.off()

# The six dates whose shape changed most ----
png(file.path(out_dir,'date_shapes.png'),height=6,width=9,units='in',res=300)
par(mfrow=c(2,3),mar=c(4,4,3,1))
worst <- order(agr$agreement)[1:6]
for (i in worst)
{
	g <- cl$calmatrix[,i]
	yy <- 3000:2000
	dens <- density(theta.posterior[,i],from=2000,to=3000,n=1001)
	plot(yy,g/max(g),type='l',xlim=c(2900,2100),ylim=c(0,1.05),xlab='BP',ylab='scaled density',col='lightblue4',main=sprintf('Date %d, agreement %.0f',i,agr$agreement[i]))
	polygon(c(yy,rev(yy)),c(g/max(g),rep(0,length(yy))),border=NA,col=adjustcolor('lightblue',1))
	lines(dens$x,dens$y/max(dens$y),lwd=2,col='darkgreen')
	abline(v=calendar.dates[i],lty=2)
}
dev.off()

# Agreement against how far the date moved ----
png(file.path(out_dir,'agreement.png'),height=5,width=6,units='in',res=300)
par(mar=c(4,4,2,1))
plot(agreement$shift,agreement$agreement,pch=19,col=adjustcolor('darkgreen',0.6),xlab='posterior median minus calibrated median (yr)',ylab='agreement index',main='Do the modelled dates still agree')
abline(h=60,lty=2)
dev.off()

cat('Done. Output in',out_dir,'\n')
