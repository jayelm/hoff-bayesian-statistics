library(parallel)

# Initial stuff borrowed from peter hoff because this data munging is unclear
dtmp = as.matrix(read.table(url("http://www.stat.washington.edu/people/pdhoff/Book/Data/data/vostok.1999.temp.dat"), header = TRUE))
dco2 = as.matrix(read.table(url('http://www.stat.washington.edu/people/pdhoff/Book/Data/data/vostok.icecore.co2.dat'), header = TRUE))

dtmp[,2]<- -dtmp[,2]
dco2[,2]<- -dco2[,2]
library(nlme)

#### get evenly spaced temperature points
ymin<-max( c(min(dtmp[,2]),min(dco2[,2])))
ymax<-min( c(max(dtmp[,2]),max(dco2[,2])))
n<-200
syear<-seq(ymin,ymax,length=n)
dat<-NULL
for(i in 1:n) {
 tmp<-dtmp[ dtmp[,2]>=syear[i] ,]
 dat<-rbind(dat,  tmp[dim(tmp)[1],c(2,4)] )
               }
dat<-as.matrix(dat)
####

####
dct<-NULL
for(i in 1:n) {
  xc<-dco2[ dco2[,2] < dat[i,1] ,,drop=FALSE]
  xc<-xc[ 1, ]
  dct<-rbind(dct, c( xc[c(2,4)], dat[i,] ) )
               }



dct<-dct[,c(3,2,4)]
colnames(dct)<-c("year","co2","tmp")
rownames(dct)<-NULL
dct<-as.data.frame(dct)
# Finally, dct is the final data frame

# Functions
tr = function(m) sum(diag(m))
inv = solve

rmvnorm = function(n, mu, Sigma) {
  p = length(mu)
  res = matrix(0, nrow = n, ncol = p)
  if (n > 0 & p > 0) {
    E = matrix(rnorm(n * p), n, p)
    res = t(t(E %*% chol(Sigma)) + c(mu))
  }
  res
}


GLOBAL.S = 1000

# Data setup
n = dim(dct)[1]
y = dct[, 3] # i.e. predicted temperatures
# Feature vector: intercept (so x_1 = 1, always) and then co2 itself
X = cbind(rep(1, n), dct[, 2])
DY<-abs(outer( (1:n),(1:n) ,"-"))

# Initial estimates from fit
lmfit = lm(y~-1+X)
fit.gls = gls(y~X[,2], correlation=corARMA(p=1), method="ML")

beta.init = lmfit$coef
s2.init = summary(lmfit)$sigma^2
phi.init = acf(lmfit$res,plot=FALSE)$acf[2]

# Prior (diffuse)
nu0 = 1
s20 = 1
T0 = diag(1/1000,nrow=2)

# TODO: Implement thinning
do.mh.sample = function(S, burnin) {
  # Set initial guesses
  beta = beta.init
  s2 = s2.init
  phi = phi.init

  OUT = matrix(nrow = S - burnin, ncol = length(c(beta, s2, phi)))

  # Note: no acceptance ratio statistics
  for(s in 1:S) {

    # 1: Sample new beta
    Cor = phi^DY
    iCor = inv(Cor)
    V.beta = inv( t(X)%*%iCor%*%X/s2 + T0)
    E.beta = V.beta%*%( t(X)%*%iCor%*%y/s2  )
    beta = t(rmvnorm(1,E.beta,V.beta)  )

    # 2: sample of inverse gamma
    s2 = 1 / rgamma(1,(nu0 + n)/2,(nu0*s20+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )

    # Update p
    # How to use proposals with bounds??
    # 3a. proposal
    phi.p = abs(runif(1,phi-.1,phi+.1))
    phi.p = min( phi.p, 2-phi.p)

    # 3b. Compute acceptance ratio
    lr = -.5*( determinant(phi.p^DY,log=TRUE)$mod -
               determinant(phi^DY,log=TRUE)$mod  +
     tr( (y-X%*%beta)%*%t(y-X%*%beta)%*%(inv(phi.p^DY) -inv(phi^DY)) )/s2 )

    if (log(runif(1)) < lr) {
      phi = phi.p
    }

    if (s > burnin) {
      OUT[s - burnin, ] = c(beta, s2, phi)
    }
  }
  OUT
}

S.global = 25000

outs = mclapply(X = rep(S.global, 4), FUN = do.mh.sample, mc.cores = 4, burnin = 5000)

outs = do.call(rbind, outs)

save(outs, file = "./icecore_mcmc")
