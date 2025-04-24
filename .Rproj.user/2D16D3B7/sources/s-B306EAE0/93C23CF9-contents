ENTdone <- function(M, code, model, Xobs, Yobs, Dx, Dtheta, sigeps, 
                  ctype="VAR", rprior, dprior, Lopt=1e4, Lmc=1e4, Lmcmc=1e4){
  D = model@X
  fD = model@y

  for(k in 1:M){
    thetak = optimf(fn=function(t){Entropy(t, model, Yobs, Xobs, sigeps, L=Lmc)}, Domain=Dtheta, L=Lopt, tomax=TRUE)
    thetak = matrix(thetak, nrow=1)

    xk = CCx(thetak, model, Xobs, Dtheta, type=ctype)
    xk = matrix(xk, nrow=1)
  
    Dnew = cbind(xk, thetak)
    fDnew = code(xk, thetak)
    model <- update(model, newX = Dnew, newy = fDnew, cov.reestim = TRUE, 
                    trend.reestim = TRUE, nugget.reestim = TRUE)
    
    D = rbind(D, Dnew)
    fD = rbind(fD, fDnew)
    
    rate <- (k/M)*100
    cat(sprintf("Progression : %.2f%%\r", rate))
    flush.console()
  } 

  Chain <- MetropInGibbs(logfd=function(theta){Logpost(theta, Yobs, Xobs, model, sigeps, Dtheta)}, 
                         maxiter=Lmcmc, xinit=rowMeans(Dtheta), sig=0.25, Domaine=Dtheta)
  
  thetamap = colMeans(Chain)
  
  return(list(D = D, fD=fD, model = model, chain = Chain, thetamap = thetamap))
}
