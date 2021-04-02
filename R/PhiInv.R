PhiInv <- function(onec, M1, M0,
                   all_Z_D, all_Z_C,
                   rho0 = 0.95) {
    ### solve for beta
    thisX <- as.matrix(cbind(1,all_Z_D))
    nvar <- ncol(thisX)
    onemorerow <- onec[1:nvar]/(rho0-1)*length(M1)
    newX <- rbind(thisX, onemorerow)
    newmarker <- c(M1, 9999)
    rqfit <- rq(newmarker ~ newX - 1, tau = 1-rho0)
    newbeta <- rqfit$coefficients

    ### solve for phi
    ctrlZ <- as.matrix(cbind(1, all_Z_C))
    ctrl_threshold_q <- ctrlZ %*% newbeta
    newphi <- (sum(M0 <= ctrl_threshold_q))/length(M0) - onec[nvar+1]

    return(c(newbeta, newphi))
}
