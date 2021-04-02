GetBootstrapSD <- function(diseaseData,
                           controlData,
                           userFormula,
                           outform,
                           fixSens,
                           Nbootstrap = 100,
                           verbose = TRUE) {

    allrec_bt <- matrix(NA, Nbootstrap, length(outform$allvar)+2)

    if(verbose) {
        message("Estimating standard error using bootstrap...")
    }

    for (nb in 1:Nbootstrap) {
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
        expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", 1-fixSens, ", data = diseaseData_bt)")
        eval(parse(text = expr1))
        est_b_bt <- as.numeric(rqfit$coefficients)
        ctrl_threshold_bt <- cbind(1, Z_C_bt) %*% est_b_bt
        case_threshold_bt <- cbind(1, Z_D_bt) %*% est_b_bt
        onephi <- sum(M0_bt <= ctrl_threshold_bt)/length(M0_bt)
        allrec_bt[nb, 1:(length(outform$allvar)+1)] <- est_b_bt
        allrec_bt[nb, length(outform$allvar)+2] <- onephi
    }

    return(apply(allrec_bt, 2, sd))
}
