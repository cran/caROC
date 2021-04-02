CalculateSD_L <- function(M1, M0,
                          all_Z_D,
                          all_Z_C,
                          case_threshold_q,
                          ctrl_threshold_q,
                          pEst,
                          phi_Z, 
                          rho0 = 0.95) {
    
    ### accounting for covariates
    caseZ <- as.matrix(cbind(1, all_Z_D))
    ctrlZ <- as.matrix(cbind(1, all_Z_C))
    nvar <- ncol(caseZ)
    n1 <- nrow(caseZ)
    n0 <- nrow(ctrlZ)
    
    part1 <- matrix(0, nvar, nvar)
    for(ii in 1:n1) {
        oneZ <- caseZ[ii, ]
        tt1 <- as.numeric(M1[ii] > case_threshold_q[ii]) - rho0
        part1 <- part1 + oneZ %*% t(oneZ) * (tt1)^2
    }
    part2 <- matrix(0, nvar, nvar)
    for(jj in 1:n0) {
        oneZ <- ctrlZ[jj, ]
        phi_oneC <- exp(sum(oneZ*pEst[nvar+1:nvar]))/(1+exp(sum(oneZ*pEst[nvar+1:nvar])))
        tt2 <- as.numeric(M0[jj] <= ctrl_threshold_q[jj]) - phi_oneC
        part2 <- part2 + oneZ %*% t(oneZ) * (tt2)^2
    }
    
    Sigma <- matrix(0, nvar*2, nvar*2)
    Sigma[1:nvar, 1:nvar] <- part1/n1^2
    Sigma[nvar+1:nvar, nvar+1:nvar] <- part2/n0^2
    Sigma.cho <- GetEDecom(Sigma)
    
    storeInv <- matrix(0, nvar*2, nvar*2)
    for(i in 1:(nvar*2)) {
        storeInv[,i] <- PhiInv_L(Sigma.cho[,i],
                                 M1, M0,
                                 all_Z_D, all_Z_C,
                                 pEst,
                                 rho0 = rho0)
    }
    D <- storeInv - pEst
    allsd <- sqrt(diag(D%*%t(D)))
    
    PhiVar <- GetPhiVar(pEst, D, phi_Z)
    allsd <- c(allsd, sqrt(PhiVar))
    
    return(allsd)
}