library(DiceKriging)
library(viridis)

# Chargement des codes
source("C:/Users/barrya/Desktop/DONEcal/Codes/KL.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/SSW.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/Entropy.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/Utiles.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/SUR.R")

# Exemple en dimension 2
fcode <- function(x, theta) {
  o = ((6*x-2)**2)*sin(theta*x-4)
  return(o)
}
# Domaine
Dx = matrix(c(0,1), ncol=2)
Dtheta = matrix(c(5,15), ncol=2)

# Observations Physiques
n = 10
sigeps=0.01
Xobs = matrix(seq(0,1, len=n), ncol=1)
Yobs = apply(array(Xobs), 1, function(u){fcode(u,theta=12)}) + rnorm(n, mean=0, sd=sigeps)


# Plot observations physiques 
xx = seq(0,1, len=100)
yreal = apply(array(xx), 1, function(x){fcode(x,theta=12)})

setwd("C:/Users/barrya/Desktop/DONEcal/Graphiques")

png("Yobs2D.png",width  = 2400,height = 1800, res=300,type="cairo")
plot(xx, yreal,col  = "blue", lwd  = 2,type = "l", xlab = expression(x),
     ylab = expression(y[real](x)),  main = "")
points(Xobs, Yobs,col   = "red", pch   = 17,type  = "p")
legend("topleft",legend = c("Yreal", "Observations"),col    = c("blue", "red"),
       lty    = c(1, NA), pch    = c(NA, 17), lwd= c(2, NA), pt.cex = c(NA, 1), bty    = "n")    
dev.off()

# Plot code de calculs
tt = seq(5,15, len=100)
zz = t(apply(array(xx), 1, function(x){apply(array(tt), 1, function(theta){fcode(x,theta)})}))

png("Code2D.png",width  = 1800,height = 1800, res=300, type="cairo")
par(mfrow = c(1,1), mar = c(4,4,2,1))
image(xx, tt, zz,col = viridis(100), xlab = expression(x), 
      ylab = expression(theta),
      main = "Code de calculs")
contour(xx, tt, zz,levels = pretty(range(zz), 10), add = TRUE, lwd = 0.5)
dev.off()

# Observations numeriques initiales
Dinit = LHSD(20, rbind(Dx,Dtheta))
fDinit = apply(Dinit, 1, function(u){fcode(u[1],u[2])})

# Modele GP initial
model <- km(formula=~1+.+.^2, design=Dinit, response=fDinit, covtype="matern5_2")

# Distribution a priori
dprior <- function(theta, Dtheta){
  o <- prod(1 /  (Dtheta[,2] - Dtheta[,1]))
  return(o)
}

rprior <- function(N, Dtheta){
  o <- matrix(runif(N*nrow(Dtheta), rep(Dtheta[,1], each = N), 
                    rep(Dtheta[,2], each = N)), nrow = N, byrow = FALSE)
  return(o)
}

# Vraie densité a posteriori
logp <- function(theta,Yobs, Xobs, code=fcode, sigeps, Dtheta){
  if(all(theta >= Dtheta[, 1] & theta <= Dtheta[, 2])){
    n = nrow(Xobs)
    Dxtheta =  Dbind(Xobs,theta)
    DD = Yobs - apply(Dxtheta, 1, function(u){code(u[1], u[2]) })
    V = diag(sigeps^2, n)
    Cte = - 0.5*n*log(2*pi) - 0.5*log(det(V)+1e-16) + log(dprior(theta, Dtheta)) 
    R = Cte - 0.5*t(DD)%*%solve(V)%*%DD
  }else{
    R = -Inf
  }
  return(R)
}

SamplePost <- MetropInGibbs(logfd=function(theta){logp(theta,Yobs, Xobs, code=fcode, sigeps, Dtheta)}, 
                        maxiter=1e4, xinit=rowMeans(Dtheta), sig=0.25, Domaine=Dtheta)

TruePDF = density(SamplePost, from=5, to=15, bw=0.5)

# Planification d'experiences numeriques
# KL 
R = KLdone(M=20, code=fcode, model, Xobs, Yobs, Dx, Dtheta, sigeps=sigeps, 
           ctype="VAR", rprior, dprior, Lmc=1000, Lopt=1000, Lmcmc=1e4)

ApproxPDF = density(R$chain, from=5, to=15, bw=0.5)

png("KLPDF.png",width  = 1800,height = 1800, res=300, type="cairo")
plot(TruePDF$x,TruePDF$y, col='blue', lwd=2, type='l', main='KL',
     xlab=expression(theta), ylab="PDF")
lines(ApproxPDF$x, ApproxPDF$y, col='red', lwd=2, lty=2)
abline(v=12, col='black', lwd=2, lty=3)
legend("topleft",
       legend = c("Vraie densité", "Approximation", "Vraie valeur"),
       col    = c("blue", "red","black"),
       lty    = c(1,2, 3),
       lwd    = c(2,2,2),
       bty    = "n")
dev.off()

# SSW
R =  SSWdone(M=20, code=fcode, model, Xobs, Yobs, Dx, Dtheta, sigeps, ctype="VAR",
             rprior, dprior, Lopt=100, Lmcmc = 1e4)

ApproxPDF = density(R$chain, from=5, to=15, bw=0.5)

png("SSWPDF.png",width  = 1800,height = 1800, res=300, type="cairo")
plot(TruePDF$x,TruePDF$y, col='blue', lwd=2, type='l', main='SSW',
     xlab=expression(theta), ylab="PDF")
lines(ApproxPDF$x, ApproxPDF$y, col='red', lwd=2, lty=2)
abline(v=12, col='black', lwd=2, lty=3)
legend("topleft",
       legend = c("Vraie densité", "Approximation", "Vraie valeur"),
       col    = c("blue", "red","black"),
       lty    = c(1,2, 3),
       lwd    = c(2,2,2),
       bty    = "n")
dev.off()

# ENTROPY
R = ENTdone(M=20, code=fcode, model, Xobs, Yobs, Dx, Dtheta, sigeps, 
            ctype="VAR", rprior, dprior, Lopt=1000, Lmc=1000, Lmcmc=1e4)

ApproxPDF = density(R$chain, from=5, to=15, bw=0.5)

png("ENTPDF.png",width  = 1800,height = 1800, res=300, type="cairo")
plot(TruePDF$x,TruePDF$y, col='blue', lwd=2, type='l', main='ENT',
     xlab=expression(theta), ylab="PDF")
lines(ApproxPDF$x, ApproxPDF$y, col='red', lwd=2, lty=2)
abline(v=12, col='black', lwd=2, lty=3)
legend("topleft",
       legend = c("Vraie densité", "Approximation", "Vraie valeur"),
       col    = c("blue", "red","black"),
       lty    = c(1,2, 3),
       lwd    = c(2,2,2),
       bty    = "n")
dev.off()

# SUR 1
R = SUR(M=20, code, model, Xobs, Yobs, Dx, Dtheta, sigeps, SURtype="SUR1",
        rprior, dprior, Lmc=35, Lopt=1e3, Lmcmc=1e4)
ApproxPDF = density(R$chain, from=5, to=15, bw=0.5)

png("SURPDF.png",width  = 1800,height = 1800, res=300, type="cairo")
plot(TruePDF$x,TruePDF$y, col='blue', lwd=2, type='l', main='SUR',
     xlab=expression(theta), ylab="PDF")
lines(ApproxPDF$x, ApproxPDF$y, col='red', lwd=2, lty=2)
abline(v=12, col='black', lwd=2, lty=3)
legend("topleft",
       legend = c("Vraie densité", "Approximation", "Vraie valeur"),
       col    = c("blue", "red","black"),
       lty    = c(1,2, 3),
       lwd    = c(2,2,2),
       bty    = "n")
dev.off()

