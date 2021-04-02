AdjSpec_AllPoints_none_L <- function(diseaseData,
                              controlData,
                              userFormula,
                              target_covariates) {

    outform <- ParseuserFormula_L(userFormula = userFormula)
    M0 <- eval(parse(text = paste("controlData$",outform[[1]])))
    M1 <- eval(parse(text = paste("diseaseData$",outform[[1]])))
    if (length(outform[[2]]) == 1) {
        expr2 <- paste0("controlData$", outform[[2]])
        expr3 <- paste0("diseaseData$", outform[[2]])
    } else if (length(outform[[2]]) > 1) {
        expr2 <- paste0("controlData[, ", outform[2], "]")
        expr3 <- paste0("diseaseData[, ", outform[2], "]")
    }
    all_Z_C <- eval(parse(text = expr2))
    all_Z_D <- eval(parse(text = expr3))

    rqfit <- NULL
    expr1 <- paste0("rqfit <- quantreg::rq(", userFormula, ", tau = -1, data = diseaseData)")
    eval(parse(text = expr1))

    all_tau <- rqfit$sol[1, ]
    sens_vec <- 1 - all_tau
    est_b <- rqfit$sol[4:(nrow(rqfit$sol)),]

    est_gamma <- matrix(NA, 3, length(all_tau))
    for(i in 1:dim(est_b)[2]) {
        control_threshold_q <- as.matrix(cbind(1, all_Z_C)) %*% est_b[,i]
        ctrl_Y <- as.numeric(M0 <= control_threshold_q)
        C_data_new <- as.data.frame(cbind(ctrl_Y, all_Z_C))
        expr2 <- paste0("mylogit <- glm( ctrl_Y ~", outform$rightside, ", data = C_data_new, family = 'binomial')")
        eval(parse(text = expr2))
        est_gamma[,i] <- as.numeric(mylogit$coefficients)
    }

    if(!is.matrix(target_covariates)) {
        eval_phi <- matrix(NA, 1, length(all_tau))

        for(i in 1:length(all_tau)) {
            timeout <- sum(as.numeric(target_covariates)*as.numeric(est_gamma[,i]))
            eval_phi[1,i] <- exp(timeout)/(1+exp(timeout))
        }
    } else if(is.matrix(target_covariates)) {
        eval_phi <- matrix(NA, nrow(target_covariates), length(all_tau))

        for(jj in 1:nrow(target_covariates)) {
            for(i in 1:length(all_tau)) {
                timeout <- sum(as.numeric(target_covariates[jj, ])*as.numeric(est_gamma[,i]))
                eval_phi[jj,i] <- exp(timeout)/(1+exp(timeout))
            }
        }
        rownames(eval_phi) <- apply(target_covariates, 1, paste, collapse = "_")
    }

    return(list(sensitivity = sens_vec,
           specificity = eval_phi,
           mono_adj = "none"))

}
