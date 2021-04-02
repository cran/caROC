GetEst_SE <- function(diseaseData,
                     controlData,
                     userFormula,
                     fixSS,
                     whichSE = "sample",
                     nbootstrap = 100,
                     verbose = TRUE) {

    outform <- ParseuserFormula(userFormula = userFormula)
    M0 <- eval(parse(text = paste("controlData$",outform[[1]])))
    M1 <- eval(parse(text = paste("diseaseData$",outform[[1]])))

    rqfit <- NULL
    expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", 1-fixSS, ", data = diseaseData)")
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
    phi_q <- sum(M0 <= ctrl_threshold_q)/length(M0)

    allnames <- c("Intercept", outform$allvar, "AdjSpec")
    pEst <- c(est_b, phi_q)
    names(pEst) <- allnames

    if (whichSE == "sample") {
        pSD <- CalculateSD(M1, M0, all_Z_D, all_Z_C,
                           case_threshold_q, ctrl_threshold_q,
                           pEst, rho0 = fixSS)
    } else if (whichSE == "bootstrap") {
        pSD <- GetBootstrapSD(diseaseData,
                              controlData,
                              userFormula,
                              outform,
                              fixSS,
                              nbootstrap,
                              verbose)
    } else {
        stop("whichSE can only be sample or bootstrap!")
    }
    names(pSD) <- allnames

    return(list(Estimates = pEst,
           SE = pSD,
           M1 = M1,
           M0 = M0,
           all_Z_C = all_Z_C,
           all_Z_D = all_Z_D))
}
