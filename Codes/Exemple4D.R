library(DiceKriging)
library(viridis)
library(MASS)

# Chargement des codes
source("C:/Users/barrya/Desktop/DONEcal/Codes/KL.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/SSW.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/Entropy.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/Utiles.R")
source("C:/Users/barrya/Desktop/DONEcal/Codes/SUR.R")

# Exemple en dimension 4
# Code de calculs
fcode <- function(x, theta) {
  o = theta[1]*x[1]^2 + theta[2]*x[2]^2 + theta[1]*theta[2]*x[1]*x[2] + 5
  return(o)
}

# Domaine
Dx = matrix(c(-2,-2,2,2), ncol=2)
Dtheta = matrix(c(0,0,5,5), ncol=2)

# Observations physiques 
n=20
sigeps=0.5
Xobs = LHSD(n,Dx)
Yobs = apply(Xobs, 1, function(u){fcode(u,theta=c(2,3))}) + rnorm(n, mean=0, sd=sigeps)

# Observations numeriques initiales
Dinit = LHSD(25, rbind(Dx,Dtheta))
fDinit = apply(Dinit, 1, function(u){fcode(u[1:2],u[3:4])})

# Modele GP initial
model <- km(formula=~., design=Dinit, response=fDinit, covtype="matern5_2")

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
    DD = Yobs - apply(Dxtheta, 1, function(u){code(u[1:2], u[3:4]) })
    V = diag(sigeps^2, n)
    Cte = - 0.5*n*log(2*pi) - 0.5*log(det(V)+1e-16) + log(dprior(theta, Dtheta)) 
    R = Cte - 0.5*t(DD)%*%solve(V)%*%DD
  }else{
    R = -Inf
  }
  return(R)
}


Chain <- MetropInGibbs(logfd=function(theta){logp(theta,Yobs, Xobs, code=fcode, sigeps, Dtheta)}, 
                        maxiter=1e4, xinit=rowMeans(Dtheta), sig=0.25, Domaine=Dtheta)
TruePDF <- kde2d(Chain[,1], Chain[,2], n = 1000,
             lims = c(range(c(1.8,2.4)), range(c(2.4,3.4))))

png("4DTruePDF.png",width  = 2400,height = 1800, res=300,type="cairo")
image(TruePDF,xlim=c(1.8,2.4), ylim=c(2.4,3.4),
      col  = viridis(1000),
      xlab = expression(theta[1]),
      ylab = expression(theta[2]),
      main = "Vraie PDF")
contour(TruePDF,
        add    = TRUE,
        drawlabels = FALSE,
        nlevels   = 15,
        lwd       = 1)
dev.off()


# Planification d'expériences numeriques : approche KL 
R = KLdone(M=25, code=fcode, model, Xobs, Yobs, Dx, Dtheta, sigeps=sigeps, 
           ctype="VAR", rprior, dprior, Lmc=200, Lopt=200, Lmcmc=1e4)



ApproxPDF <- kde2d(R$chain[,1], R$chain[,2], n = 1000,
                 lims = c(range(c(1.8,2.4)), range(c(2.4,3.4))))

png("4DApproxPDF.png", width  = 2400, height = 1800, res=300,type="cairo")
image(ApproxPDF,xlim=c(1.8,2.4), ylim=c(2.4,3.4),
      col  = viridis(1000),
      xlab = expression(theta[1]),
      ylab = expression(theta[2]),
      main = "PDF approchée")
contour(ApproxPDF,
        add    = TRUE,
        drawlabels = FALSE,
        nlevels   = 15,
        lwd       = 1)
dev.off()


png("4D_2PDF.png", width  = 2400, height = 1800, res=300,type="cairo")
par(mfrow=c(1,2))
image(TruePDF,xlim=c(1.8,2.4), ylim=c(2.4,3.4),
      col  = viridis(1000),
      xlab = expression(theta[1]),
      ylab = expression(theta[2]),
      main = "Vraie PDF")
contour(TruePDF,
        add    = TRUE,
        drawlabels = FALSE,
        nlevels   = 15,
        lwd       = 1)
image(ApproxPDF,xlim=c(1.8,2.4), ylim=c(2.4,3.4),
      col  = viridis(1000),
      xlab = expression(theta[1]),
      ylab = expression(theta[2]),
      main = "PDF approchée")
contour(ApproxPDF,
        add    = TRUE,
        drawlabels = FALSE,
        nlevels   = 15,
        lwd       = 1)
dev.off()

# R$thetamap
# Chain = R$chain
# 
# colMeans(Chain)
# par(mfrow=c(1,2))
# plot(1:nrow(Chain), Chain[,1], main="theta1", col='blue', type='l', lwd=2)
# plot(1:nrow(Chain), Chain[,2], main="theta2", col='red', type='l', lwd=2)
