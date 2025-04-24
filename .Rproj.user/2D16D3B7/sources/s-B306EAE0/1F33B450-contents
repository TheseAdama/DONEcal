SURdone <- function(M, code, model, Xobs, Yobs, Dx, Dtheta, sigeps, 
                    SURtype="SUR1", rprior, dprior, Lmc=1e3, Lopt=1e4, Lmcmc=1e4){
  D = model@X
  fD = model@y
  
  for(k in 1:M){
    Xnew = optimf(fn=function(xt){Jm(xt, model, Yobs, Xobs, 
                                     sigeps, SURtype, Dx, Dtheta, Lmc)}, 
                  Domain=rbind(Dx,Dtheta), L=Lopt, tomax=FALSE)
    
    
    Ynew = code(Xnew[1:ncol(Xobs)], Xnew[(ncol(Xobs)+1):ncol(D)])

    model <- update(model, newX=matrix(Xnew, nrow=1), newy = Ynew, cov.reestim = TRUE, 
                    trend.reestim=TRUE, nugget.reestim = TRUE)
    
    D = rbind(D, Xnew)
    fD = rbind(fD, Ynew)
    
    rate <- (k/M)*100
    cat(sprintf("Progression : %.2f%%\r", rate))
    flush.console()
  }
  Chain <- MetropInGibbs(logfd=function(theta){Logpost(theta, Yobs, Xobs, model, sigeps, Dtheta)}, 
                         maxiter=Lmcmc, xinit=rowMeans(Dtheta), sig=0.25, Domaine=Dtheta)
  thetamap = colMeans(Chain)
  
  return(list(D = D, fD=fD, model = model, chain = Chain, thetamap = thetamap))
}