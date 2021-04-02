ParseuserFormula_L <- function(userFormula) {
    formula2 <- gsub(" ", "", userFormula)
    tt1 <- unlist(strsplit(formula2, split = "~"))
    allvar <- unlist(strsplit(tt1[2], split = "\\+"))
    rightside <- paste0(allvar, collapse = " + ")
    return(list(response = tt1[1],
                allvar = allvar,
                rightside = rightside))
}

GetEDecom <-   function(mat) {
    eout <- eigen(mat)
    L <- eout$vectors %*% diag(sqrt(eout$values))
    return(L)
}

mylogit <- function(p) {
    return(log(p/(1-p)))
}

GetPhiVar <- function(pEst, D, phi_Z) {
    nvar <- length(phi_Z)

    egamma <- pEst[nvar+1:nvar]
    gammaVar <- (D%*%t(D))[nvar+1:nvar, nvar+1:nvar]
    pptmp <- exp(phi_Z%*%egamma)/((1+exp(phi_Z%*%egamma))^2)
    phiout <- matrix(phi_Z, nrow=1) %*% gammaVar %*% matrix(t(phi_Z), ncol=1) * pptmp^2

    return(phiout)
}

PhiInv_L <- function(onec, M1, M0,
                     all_Z_D, all_Z_C,
                     pEst,
                     rho0 = 0.95) {
    n0 = nrow(all_Z_C)
    n1 = nrow(all_Z_D)
    nvar = ncol(all_Z_C)+1

    ### solve for beta
    thisX <- as.matrix(cbind(1,all_Z_D))
    nvar <- ncol(thisX)
    onemorerow <- onec[1:nvar]/(rho0-1)*length(M1)
    newX <- rbind(thisX, onemorerow)
    newmarker <- c(M1, 9999)
    rqfit <- quantreg::rq(newmarker ~ newX - 1, tau = 1-rho0)
    newbeta <- rqfit$coefficients

    ### solve for gamma
    ctrlZ <- as.matrix(cbind(1, all_Z_C))
    ctrl_Y <- as.numeric(M0 <= ctrlZ %*% newbeta)
    C_data_new <- as.data.frame(cbind(ctrl_Y, all_Z_C))
    mylogit <- glm(ctrl_Y ~ ., data = C_data_new, family = "binomial")
    start_gamma <- mylogit$coefficients

    if(all(onec[nvar+1:nvar] == 0)) {
        newgamma <- start_gamma
    } else {
        ctrlZ <- as.matrix(cbind(1, all_Z_C))
        ctrl_threshold_q <- ctrlZ %*% newbeta
        Fgammac <- function(tmpgamma, ggc) {
            Fout <- rep(0, nvar)
            for(jj in 1:n0) {
                oneZ <- ctrlZ[jj,]
                phi_oneC <- exp(ctrlZ[jj,]%*%tmpgamma)/(1+exp(ctrlZ[jj,]%*%tmpgamma))
                tt2 <- as.numeric(M0[jj] <= ctrl_threshold_q[jj]) - phi_oneC
                Fout <- Fout + oneZ * c(tt2)
            }
            Fout <- Fout/n0 - ggc
            return(Fout)
        }
        getLogLike <- function(tmpgamma) {
            llout <- 0
            for(jj in 1:n0) {
                oneZ <- ctrlZ[jj,]
                oneidx <- as.numeric(M0[jj] <= ctrl_threshold_q[jj])
                tt4 <- -log(1+exp(oneZ %*% tmpgamma)) + oneidx*oneZ%*%tmpgamma
                llout <- llout + tt4
            }
            return(tt4)
        }
        Jgamma <- function(tmpgamma) {
            Jout <- matrix(0, nvar, nvar)
            for(jj in 1:n0) {
                oneZ <- ctrlZ[jj,]
                phi_oneC2 <- exp(ctrlZ[jj,]%*%tmpgamma)/(1+exp(ctrlZ[jj,]%*%tmpgamma))^2
                Jout <- Jout + oneZ %*% t(oneZ) * c(phi_oneC2)
            }
            tmpout <- -Jout/n0
            return(tmpout)
        }
        myNR <- function(start_gamma, diff_thres, ssc) {
            old_gamma <- start_gamma
            thisdiff <- 100
            old_ll <- getLogLike(old_gamma)
            ncount <- 1
            while(thisdiff > diff_thres) {
                Fout <- Fgammac(old_gamma, ssc)
                Jout <- Jgamma(old_gamma)
                ssout <- solve(Jout)
                delta <- solve(Jout) %*% Fout
                thisgamma <- old_gamma - delta
                new_ll <- getLogLike(thisgamma)
                while((new_ll <= old_ll) & (mean(abs(delta)) > 1e-3)) {
                    delta <- delta/2
                    thisgamma <- old_gamma - delta
                    new_ll <- getLogLike(thisgamma)
                }
                thisdiff <- mean(abs(thisgamma - old_gamma))
                old_ll <- new_ll
                old_gamma <- thisgamma
                ncount <- ncount + 1
                if(is.na(thisdiff)) {
                    thisgamma <- start_gamma
                    warning("Inference procedure does not converge: possibly sample size is too small.")
                    break
                }
            }
            return(thisgamma)
        }
        newgamma <- myNR(start_gamma, diff_thres = 1e-3, onec[nvar+1:nvar])
    }

    return(c(newbeta, newgamma))
}
