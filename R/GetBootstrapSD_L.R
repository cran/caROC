GetBootstrapSD_L <- function(diseaseData,
                           controlData,
                           userFormula,
                           target_covariates,
                           fixSS,
                           nbootstrap = 100,
                           verbose = TRUE) {

    outform <- ParseuserFormula_L(userFormula = userFormula)

    if(!is.matrix(target_covariates)) {
        neZ <- 1
    } else {
        neZ <- nrow(target_covariates)
    }
    allrec_bt <- matrix(NA, nbootstrap, 2+length(outform$allvar)*2+neZ)

    if(verbose) {
        message("Estimating standard error using bootstrap...")
    }

    for (nb in 1:nbootstrap) {
        oneindx1 <- sample(1:nrow(controlData), nrow(controlData), replace = TRUE)
        oneindx2 <- sample(1:nrow(diseaseData), nrow(diseaseData), replace = TRUE)

        controlData_bt <- controlData[oneindx1, ]
        diseaseData_bt <- diseaseData[oneindx2, ]

        if (length(outform[[2]]) == 1) {
            expr2 <- paste0("controlData_bt$", outform[[2]])
            expr3 <- paste0("diseaseData_bt$", outform[[2]])
        } else if (length(outform[[2]]) > 1) {
            expr2 <- paste0("controlData_bt[, ", outform[2], "]")
            expr3 <- paste0("diseaseData_bt[, ", outform[2], "]")
        }
        Z_C_bt <- as.matrix(eval(parse(text = expr2)))
        Z_D_bt <- as.matrix(eval(parse(text = expr3)))

        M0_bt <- eval(parse(text = paste("controlData_bt$",outform[[1]])))
        M1_bt <- eval(parse(text = paste("diseaseData_bt$",outform[[1]])))

        rqfit <- NULL
        expr1 <- paste0("rqfit <- quantreg::rq(", userFormula, ", tau = ", 1-fixSS, ", data = diseaseData_bt)")
        eval(parse(text = expr1))
        est_b_bt <- as.numeric(rqfit$coefficients)
        ctrl_threshold_bt <- cbind(1, Z_C_bt) %*% est_b_bt
        case_threshold_bt <- cbind(1, Z_D_bt) %*% est_b_bt

        ctrl_Y <- as.numeric(M0_bt <= ctrl_threshold_bt)
        C_data_newbt <- as.data.frame(cbind(ctrl_Y, Z_C_bt))
        expr2 <- paste0("mylogit <- glm(ctrl_Y~., data = C_data_newbt, family = 'binomial')")
        eval(parse(text = expr2))
        est_gamma_bt <- as.numeric(mylogit$coefficients)
        pEst_bt <- c(est_b_bt, est_gamma_bt)

        allrec_bt[nb, 1:(length(outform$allvar)*2+2)] <- pEst_bt

        if(!is.matrix(target_covariates)) {
            phi_vec <- NA
            phi_vec <- exp(target_covariates%*%est_gamma_bt)/(1+exp(target_covariates%*%est_gamma_bt))
        } else {
            phi_vec <- rep(NA, nrow(target_covariates))
            for(kk in 1:nrow(target_covariates)) {
                phi_vec[kk] <- exp(target_covariates[kk,]%*%est_gamma_bt)/(1+exp(target_covariates[kk,]%*%est_gamma_bt))
            }
        }
        allrec_bt[nb, (length(outform$allvar)*2+2)+1:neZ] <- phi_vec
    }

    return(apply(allrec_bt, 2, sd))
}
