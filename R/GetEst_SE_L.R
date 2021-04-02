GetEst_SE_L <- function(diseaseData,
                        controlData,
                        userFormula,
                        fixSS,
                        target_covariates,
                        whichSE = "sample",
                        nbootstrap = 100,
                        verbose = TRUE) {

    outform <- ParseuserFormula_L(userFormula = userFormula)
    M0 <- eval(parse(text = paste("controlData$",outform[[1]])))
    M1 <- eval(parse(text = paste("diseaseData$",outform[[1]])))

    rqfit <- NULL
    expr1 <- paste0("rqfit <- quantreg::rq(", userFormula, ", tau = ", 1-fixSS, ", data = diseaseData)")
    eval(parse(text = expr1))
    est_b <- as.numeric(rqfit$coefficients)
    if (length(outform[[2]]) == 1) {
        expr2 <- paste0("controlData$", outform[[2]])
        expr3 <- paste0("diseaseData$", outform[[2]])
    } else if (length(outform[[2]]) > 1) {
        expr2 <- paste0("controlData[, ", outform[2], "]")
        expr3 <- paste0("diseaseData[, ", outform[2], "]")
    }
    all_Z_C <- eval(parse(text = expr2))
    all_Z_D <- eval(parse(text = expr3))

    case_threshold_q <- as.matrix(cbind(1, all_Z_D)) %*% est_b
    ctrl_threshold_q <- as.matrix(cbind(1, all_Z_C)) %*% est_b

    ctrl_Y <- as.numeric(M0 <= ctrl_threshold_q)
    C_data_new <- as.data.frame(cbind(ctrl_Y, all_Z_C))
    expr2 <- paste0("mylogit <- glm( ctrl_Y ~", outform$rightside, ", data = C_data_new, family = 'binomial')")
    eval(parse(text = expr2))
    est_gamma <- as.numeric(mylogit$coefficients)

    pEst <- c(est_b, est_gamma)

    if(!is.matrix(target_covariates)) {
        nZ <- 1
        pSD <- rep(NA, length(pEst)+nZ)

        phi_vec <- exp(target_covariates%*%est_gamma)/(1+exp(target_covariates%*%est_gamma))
        if (whichSE == "sample") {
            pSD <- CalculateSD_L(M1, M0, all_Z_D, all_Z_C,
                                 case_threshold_q, ctrl_threshold_q,
                                 pEst, target_covariates, rho0 = fixSS)
        } else if (whichSE == "bootstrap") {
            pSD <- GetBootstrapSD_L(diseaseData,
                                    controlData,
                                    userFormula,
                                    target_covariates,
                                    fixSS,
                                    nbootstrap,
                                    verbose)
        }


    } else {
        nZ <- nrow(target_covariates)
        pSD <- rep(NA, length(pEst)+nZ)

        phi_vec = phi_SE_vec = rep(NA, nrow(target_covariates))
        for(kk in 1:nrow(target_covariates)) {
            phi_vec[kk] <- exp(target_covariates[kk,]%*%est_gamma)/(1+exp(target_covariates[kk,]%*%est_gamma))
            phi_Z <- target_covariates[kk,]
            if (whichSE == "sample") {
                SEout <- CalculateSD_L(M1, M0, all_Z_D, all_Z_C,
                                       case_threshold_q, ctrl_threshold_q,
                                       pEst, phi_Z, rho0 = fixSS)
            } else if (whichSE == "bootstrap") {
                SEout <- GetBootstrapSD_L(diseaseData,
                                        controlData,
                                        userFormula,
                                        target_covariates,
                                        fixSS,
                                        nbootstrap,
                                        verbose)
            }

            phi_SE_vec[kk] <- SEout[length(SEout)]
            if(kk == 1) {
                pEst_SE <- SEout[1:(length(est_b)+length(est_gamma))]
            }
        }
        pSD <- c(pEst_SE, phi_SE_vec)
    }

    allnames <- c(paste0(c(rep("beta_", length(est_b)),rep("gamma_", length(est_gamma))),
                         rep(c("Intercept", outform$allvar), 2)),
                  paste0(rep("AdjSpec", nZ), 1:nZ))
    pEst <- c(pEst, phi_vec)
    names(pEst) <- allnames
    names(pSD) <- allnames

    return(list(Estimates = pEst,
                SE = pSD,
                M1 = M1,
                M0 = M0,
                all_Z_C = all_Z_C,
                all_Z_D = all_Z_D))
}
