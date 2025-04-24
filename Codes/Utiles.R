# Algorithme de Metropolis dans Gibbs
MetropInGibbs <- function(logfd, maxiter, xinit, sig=0.5, Domaine=NULL) {
  d <- length(xinit)
  sig = rep(sig, d)
  if(!is.null(Domaine)){sig <- apply(Domaine, 1, function(u){return(sig*(u[2]-u[1]))})}
  samples <- matrix(NA, nrow = maxiter, ncol = d)
  samples[1, ] <- xinit
  for (i in 2:maxiter) {
    xi <- samples[i - 1, ]
    for (j in 1:d) {
      xprop <- xi
      xprop[j] <- rnorm(1, mean = xi[j], sd = sig[j])
      rr <- as.numeric(logfd(xprop) - logfd(xi))
      if (log(runif(1)) < rr) {
        xi[j] <- xprop[j]
      }
      
    }
    samples[i, ] <- xi
    if(i%%1000==0){cat("ITERATION MCMC :", i,"\n")}
  }
  
  samples <- samples[(round(0.2*maxiter) + 1):maxiter, ]
  samples <- matrix(samples, ncol=d)
  samples <- samples[seq(1, nrow(samples), by=3), ]
  samples <- matrix(samples, ncol=d)
  return(samples)
}
# Grille reguliere sur D
GridD <- function(L, D){
  R = matrix(0, ncol=nrow(D), nrow=L)
  for(i in 1:nrow(D)){
    R[,i] <- base::seq(D[i,1],D[i,2],len=L)
  }
  return(R)
}

# Plan LHSmaximin sur D
LHSD<- function(N, D){
  G <- lhs::maximinLHS(N,nrow(D))
  for(i in 1:nrow(D)){
    G[,i] <- (D[i,2]-D[i,1])*G[,i] + D[i,1]
  }
  return(as.matrix(G))
}

# Tirage uniforme dans D
Unif <- function(N, D){
  R <- matrix(0, nrow = N, ncol=nrow(D))
  for(i in 1:nrow(D)){
    R[,i] <- stats::runif(N, D[i,1], D[i,2])
  }
  return(R)
}

# Combinaison de D et x
Dbind <- function(D,x){
  mat <- cbind(D,matrix(rep(x,nrow(D)),ncol=length(x),byrow=TRUE))
  return(mat)
}

# Densité a priori
dprior <- function(theta, Dtheta){
  R <- prod(1 /  (Dtheta[,2] - Dtheta[,1]))
  return(R)
}

# Tirage dans la priore
rprior <- function(N, Dtheta){
  R <- matrix(runif(N*nrow(Dtheta), rep(Dtheta[,1], each = N), 
              rep(Dtheta[,2], each = N)), nrow = N, byrow = FALSE)
  return(R)
}

# Log de la densite a posteriori
Logpost <- function(theta,Yobs, Xobs, model, sigeps, Dtheta){
  if(all(theta >= Dtheta[, 1] & theta <= Dtheta[, 2])){
  n = nrow(Xobs)
  Ypred <- predict(model, newdata = data.frame(x = Dbind(Xobs,theta)),
                   cov.compute=TRUE , type = "SK", checkNames = FALSE)
  DD = Yobs - Ypred$mean
  V = Ypred$cov + diag(sigeps^2, n)
  
  Cte = -0.5*n*log(2*pi)- 0.5*log(det(V)+1e-16)+log(dprior(theta, Dtheta)) 
  R = Cte - 0.5*t(DD)%*%solve(V)%*%DD
  }else{
    R = -Inf
  }
  return(R)
}

# Somme des ecarts aux carrees (avec le code)
SS <- function(theta, X, Y){
  ff <- apply(X,1,function(x){code(x,theta)})
  o = sum((Y - ff)^2)
  return(o)
}

# Somme des ecarts aux carrees (avec l'emulateur)
SSk <- function(theta, model, X, Y){
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                   type = "SK")
  o = sum((Y - Ypred$mean)^2)
  return(o)
}
# Somme des ecarts aux carrees esperee (calcule par MC)
ESSk <-function(theta, model, X, Y, sigeps, L=1e3){
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                   cov.compute=TRUE, type = "SK", checkNames="FALSE")
  DD <- mvtnorm::rmvnorm(L, mean = Y - Ypred$mean, sigma = Ypred$cov  + diag(sigeps^2, nrow(Xobs)))
  E <- mean(apply(DD, 1, function(x){ return(norm(x, type="2")^2)} ))
  return(E)
}

# Criteres de selection des xk 
CCx <- function(theta, model, X, Dtheta, type="VAR"){
  
  YD <- predict(model, newdata = data.frame(x = Dbind(X,theta)),
                cov.compute=FALSE, type = "SK", checkNames = FALSE)
  
  if (type=="VAR"){
    xstar <- X[which.max(YD$sd),]
  }
  
  if (type=="VARprior"){
    C1 <- (YD$sd/max(YD$sd))^2
    
    if (length(theta)<=4){K <- 1000}else{K <- 1e4}
    
    TT <- rprior(K, Dtheta)
    mm <- apply(array(TT), 1, function(tt){
      Y <- predict(model, newdata = data.frame(x = Dbind(X, tt)), type = "SK",checkNames = FALSE)
      return(Y$mean)
    })
    
    C2 <- apply(mm, 1, var)
    C2 <- C2/max(C2)
    C <- C1*C2 
    xstar <- X[which.max(C), ]
  }
  
  return(matrix(xstar, nrow=1))
}


# Expected Improvement
EIk <- function(theta, mk, X, Y, model, Lmc=1e3){
  
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X,theta)),
                   cov.compute=TRUE, type = "SK", checkNames= FALSE)
  
  mu <- Ypred$mean
  Covmat <- Ypred$cov 
  YY <- mvtnorm::rmvnorm(Lmc, mean = mu, sigma = Covmat)
  SS <- apply(YY, 1, function(V){return(norm(Y - V, type="2")^2)})
  mm = mean(SS)
  ss = sd(SS)
  Ek = (mk - mm)*pnorm((mk-mm)/ss, mean=0,sd=1) + ss*dnorm((mk-mm)/ss, mean=0,sd=1)
  return(Ek)
}

# Optimisation gloutonne 
optimf <- function(fn, Domain, L=1e3, tomax=TRUE){
  
  XX = LHSD(L, Domain)
  fXX = apply(XX, 1, fn)
  
  if(tomax){
    R <- XX[which.max(fXX),]
  }
  else{
    R <- XX[which.min(fXX),]
  }
  return(R)
}

# Entropy  de classification 
Entropy <- function(theta, model, Yobs, Xobs, sigeps, L=Lmc){
  Yp <- predict(model, newdata = data.frame(x = Dbind(Xobs, theta)),
                cov.compute=TRUE, type = "SK", checkNames = FALSE)
  CC = Yp$cov+ diag(sigeps, length(Yobs))
  V = mvtnorm::rmvnorm(L, mean = Yobs - Yp$mean, sigma = CC)
  I = apply(V,1,function(U){t(U)%*%solve(CC)%*%U})
  
  pm = length(which(I<=qchisq(.99, df=length(Yobs))))/L
  
  if(pm==0 || pm==1){
    Em = 0
  }else{
    Em = - pm*log(pm) - (1-pm)*log(1-pm)
  }
  return(Em)
}

# Sommes ecarts aux carrees ponderes
WSS <- function(theta, model, Y, X, sigeps){
  
  Yp <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                type = "SK", cov.compute = TRUE, checkNames="FALSE")
  V = Y - Yp$mean
  CC = Yp$cov + diag(sigeps^2, length(Y))
  o = t(V)%*%solve(CC)%*%V
  return(o)
}

# Criteres SUR
## Mesure d'incertitude Um
Um <- function(fDadd, Dadd, model, Yobs, Xobs, sigeps, SURtype, Dx, 
               Dtheta, L){
  
  modelsim <- update(model, newX=Dadd, newy = fDadd, cov.reestim = TRUE,
                     trend.reestim = TRUE, nugget.reestim = TRUE)
  
  TT = LHSD(L, Dtheta)
  WW = apply(TT, 1, function(theta){Logpost(theta, Yobs, Xobs, model=modelsim, sigeps, Dtheta)})
  WW = exp(WW)
  WW = WW/sum(WW)
  
  if(SURtype=="SUR1"){R <- sum(diag(var(WW*TT)))}
  
  if(SURtype=="SUR2"){
    tmap = sum(WW*TT)
    Yp <- predict(modelsim, newdata = data.frame(x = Dbind(Xobs, tmap)), 
                  type = "SK", cov.compute = TRUE, light.return=TRUE, checkNames=FALSE)
    R <- t(Yobs - Yp$mean)%*%solve(Yp$cov + diag(sigeps^2, length(Yobs)))%*%(Yobs - Yp$mean)
  }
  
  if(SURtype=="SUR3"){
    tmap = sum(WW*TT)
    Yp <- predict(modelsim, newdata = data.frame(x = Dbind(Xobs, tmap)), 
                  type = "SK", cov.compute = TRUE, light.return=TRUE,
                  checkNames=FALSE)
    R <- sum(diag(Yp$cov + diag(sigeps^2, length(Yobs))))
  }
  
  return(R)
}
Jm <- function(xt, model, Yobs, Xobs, sigeps, SURtype, Dx, Dtheta, Lmc){
  
  xt = matrix(xt, nrow=1)

  Yp <- predict(model, newdata = data.frame(x = xt),
                cov.compute=FALSE, type = "SK", checkNames = FALSE)
  Ysim <- stats::rnorm(Lmc, mean=Yp$mean, sd=Yp$sd)
  
  RR = apply(array(Ysim), 1, function(U){ Um(fDadd=U, Dadd=xt, model, Yobs, Xobs, sigeps, SURtype, Dx,Dtheta, L=Lmc)})

  R = mean(na.omit(RR))
  
  return(mean(R))
}




