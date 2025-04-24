SS <- function(theta, X, Y){
  ff <- apply(X,1,function(x){code(x,theta)})
  o = sum((Y - ff)^2)
  return(o)
}

SSk <- function(theta, model, X, Y){
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                   type = "SK")
  o = sum((Y - Ypred$mean)^2)
  return(o)
}

ESSk <-function(theta, model, X, Y, sigY, L=1e3){
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X, theta)), cov.compute=TRUE, type = "SK")
  DD <- mvtnorm::rmvnorm(L, mean = Y - Ypred$mean, sigma = Ypred$cov  + diag(sigY^2, nrow(Xobs)))
  E <- mean(apply(DD, 1, function(x){ return(norm(x, type="2")^2)} ))
  return(E)
}

CCx <- function(theta, model, X, Dth, type="VAR", type_prior, param_prior){
  
  YD <- predict(model, newdata = data.frame(x = Dbind(X,theta)),
                cov.compute=FALSE, type = "SK", checkNames = FALSE)
  
  if (type=="VAR"){
    xstar <- X[which.max(YD$sd),]
  }
  
  if (type=="VARprior"){
    C1 <- (YD$sd/max(YD$sd))^2
    
    if (length(theta)<=4){K <- 1000}else{K <- 1e4}
    
    TT <- rprior(K, Dth, type_prior, param_prior)
    mm <- apply(array(TT),1, function(tt){
      Y <- predict(model, newdata = data.frame(x = Dbind(X, tt)),
                   type = "SK",checkNames = FALSE)
      return(Y$mean)
    })
    
    C2 <- apply(mm, 1, var)
    C2 <- C2/max(C2)
    C <- C1*C2 
    xstar <- X[which.max(C), ]
  }
  
  return(matrix(xstar, nrow=1))
}


WSS <- function(theta, model, Y, X, sigY){
  
  Yp <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                type = "SK", cov.compute = TRUE)
  R = sum(diag((Yp$cov + diag(sigY, length(Y)))))
  return(R)
}

WSS2 <- function(theta, model, Y, X, sigY){
  L=1e3
  Yp <- predict(model, newdata = data.frame(x = Dbind(X, theta)), 
                type = "SK", cov.compute = TRUE)
  YY <- mvtnorm::rmvnorm(L, mean = Y - Yp$mean, sigma = Yp$cov)
  SS <- apply(YY, 1, function(x){ return(norm(x, type="2")^2)})
  R = var(SS)
  return(R)
}

Entropy <- function(theta, model, Yobs, Xobs, sigY){
  L=1e3
  Yp <- predict(model, newdata = data.frame(x = Dbind(Xobs, theta)),
                cov.compute=TRUE, type = "SK", checkNames = FALSE)
  CC = Yp$cov+ sigY*diag(1,nrow(Yp$cov))
  V = mvtnorm::rmvnorm(L, mean = Yobs - Yp$mean, sigma = CC)
  I = apply(V,1,function(U){t(U)%*%solve(CC)%*%U})
  
  # Proba de classification 
  pm = length(which(I<=qchisq(.99, df=length(Yobs))))/L
  
  # Entropy 
  if(pm==0 || pm==1){
    Em = 0
  }else{
    Em = - pm*log(pm) - (1-pm)*log(1-pm)
  }
  return(Em)
}


optimf <- function(fn, Domain, L=1e3, tomax=TRUE,...){
  
  XX = LHSD(L,Domain)
  fXX = apply(XX, 1, fn,...)
  
  if(tomax){
    R <- XX[which.max(fXX),]
  }
  else{
    R <- XX[which.min(fXX),]
  }
  return(R)
}


Lpost <- function(theta, Yobs, Xobs, code, sigY, type_prior, param_prior, Dth){
  
  p = length(theta) 
  d = ncol(Xobs)
  YY <- apply(Dbind(Xobs,theta), 1, function(x){code(x[1:d],x[(d+1):(d+p)])})
  D = Yobs - YY
  V = diag(sigY, length(Yobs))
  pp = dprior(theta, Dth, type_prior, param_prior)
  R = -0.5*length(Yobs)*log(2*pi) - 0.5*log(det(V)) - 0.5*t(D)%*%solve(V)%*%D #+ log(pp)
  
  return(R)
}

Logpost <- function(theta, Yobs, Xobs, model, sigY, type_prior, param_prior, Dth){
  
  Ypred <- predict(model, newdata = data.frame(x = Dbind(Xobs,theta)),
                   cov.compute=TRUE , type = "SK", checkNames = FALSE)
  D = Yobs - Ypred$mean
  V = Ypred$cov + diag(sigY, nrow(Xobs))
  pp = dprior(theta, Dth, type_prior, param_prior)
  R = -0.5*length(Yobs)*log(2*pi) - 0.5*log(det(V)) - 0.5*t(D)%*%solve(V)%*%D #+ log(pp)
  
  return(R)
}


EIk <- function(theta, mk, X, Y, model, sigY, Lmc=1e4){
  
  Ypred <- predict(model, newdata = data.frame(x = Dbind(X,theta)),
                   cov.compute=TRUE, type = "SK", checkNames= FALSE)
  
  mu <- Ypred$mean
  Covmat <- Ypred$cov 
  YY <- mvtnorm::rmvnorm(Lmc, mean = mu, sigma = Covmat)
  SS <- apply(YY, 1, function(V){
    return(norm(Y - V, type="2")^2)
  })
  mm = mean(SS)
  ss = sd(SS)
  Ek = (mk - mm)*pnorm((mk-mm)/ss, mean=0,sd=1) + ss*dnorm((mk-mm)/ss, mean=0,sd=1)
  return(Ek)
}



CSUR <- function(xt, model, Yobs, Xobs, sigY, SURtype, Dx, Dth, Lmc, type_prior, param_prior){
  
  xt = matrix(xt, ncol=nrow(rbind(Dx, Dth)))
  
  if(nrow(rbind(Dx, Dth))<=4) N=30 else N=1e4
  if(nrow(Dth)<=2) L=15 else L=1e4
  
  Yp <- predict(model, newdata = data.frame(x = xt),
                cov.compute=FALSE, type = "SK", checkNames = FALSE)
  Ysim <- stats::rnorm(N, mean=Yp$mean, sd=Yp$sd)
  
  RR = apply(array(Ysim), 1, JJ, xt, model, Yobs, Xobs, sigY, SURtype, Dx, 
             Dth, L,type_prior, param_prior)
  R = mean(na.omit(RR))
  
  return(mean(R))
}

JJ <- function(Ys, xt, model, Yobs, Xobs, sigY, SURtype, Dx, 
               Dth, L, type_prior, param_prior){
  
  modelsim <- update(model, newX= xt, newy = Ys, cov.reestim = TRUE,
                     trend.reestim = TRUE, nugget.reestim = TRUE)
  
  TT = LHSD(L, Dth)
  WW = apply(TT, 1, function(theta){Logpost(theta, Yobs, Xobs, model, sigY, type_prior, param_prior, Dth)})
  WW = exp(WW-max(WW))
  WW = WW/sum(WW)
  
  if(SURtype=="SUR1"){R <- sum(diag(var(WW*TT)))}
  
  if(SURtype=="SUR2"){
    tmap = sum(WW*TT)
    Yp <- predict(modelsim, newdata = data.frame(x = Dbind(Xobs, tmap)), 
                  type = "SK", cov.compute = TRUE, light.return=TRUE,
                  checkNames=FALSE)
    R <- t(Yobs - Yp$mean)%*%solve(Yp$cov + diag(sigY, length(Yobs)))%*%(Yobs - Yp$mean)
  }
  
  if(SURtype=="SUR3"){
    tmap = sum(WW*TT)
    Yp <- predict(modelsim, newdata = data.frame(x = Dbind(Xobs, tmap)), 
                  type = "SK", cov.compute = TRUE, light.return=TRUE,
                  checkNames=FALSE)
    R <- sum(diag(Yp$cov + diag(sigY, length(Yobs))))
  }
  
  return(R)
}


optimff <- function(fn, X, Domain, L=1e3, tomax=TRUE,...){
  
  TT = LHSD(L, Domain)
  XX = c()
  for(i in 1:L){
    XX = rbind(XX, Dbind(X, TT[i]))
  }
  
  fXX = apply(XX, 1, fn,...)
  
  if(tomax){
    R <- XX[which.max(fXX),]
  }
  else{
    R <- XX[which.min(fXX),]
  }
  
  return(R)
}
