CalculateSD <- function(M1, M0,
                        all_Z_D,
                        all_Z_C,
                        case_threshold_q,
                        ctrl_threshold_q,
                        pEst,
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
    part2 <- 0
    for(jj in 1:n0) {
        tt2 <- as.numeric(M0[jj] <= ctrl_threshold_q[jj]) - pEst[nvar+1]
        part2 <- part2 + tt2^2
    }

    Sigma <- matrix(0, nvar+1, nvar+1)
    Sigma[1:nvar, 1:nvar] <- part1/n1^2
    Sigma[nvar+1, nvar+1] <- part2/n0^2
    if(0 %in% diag(Sigma)) {
        Sigma[diag(Sigma) == 0] <- 0.001
    }
    Sigma.cho <- t(chol(Sigma))

    storeInv <- matrix(0, nvar+1, nvar+1)
    for(i in 1:(nvar+1)) {
        storeInv[,i] <- PhiInv(Sigma.cho[,i],
                               M1, M0,
                               all_Z_D, all_Z_C,
                               rho0 = rho0)
    }
    D <- storeInv - pEst
    allsd <- sqrt(diag(D%*%t(D)))

    return(allsd)
}
