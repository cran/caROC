AdjSpec_AllPoints_curve <- function(diseaseData,
                                   controlData,
                                   userFormula) {

    outform <- ParseuserFormula(userFormula = userFormula)
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
    expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = -1, data = diseaseData)")
    eval(parse(text = expr1))

    all_tau <- rqfit$sol[1, ]
    sens_vec <- 1 - all_tau
    est_b <- rqfit$sol[4:(nrow(rqfit$sol)),]
    control_threshold_q <- as.matrix(cbind(1, all_Z_C)) %*% est_b
    spec_adj <- rep(0, length(all_tau))
    spec_adj <- colSums(M0 <= control_threshold_q)/length(M0)

    tmpout <- MonoRespect(tau = all_tau, orig_measure = spec_adj, startTau = 0.5)
    spec_curve <- tmpout$new_meas
    sens_curve <- 1-tmpout$tau

    return(list(sensitivity = sens_curve,
                specificity = spec_curve,
                mono_adj = "ROC"))

}
