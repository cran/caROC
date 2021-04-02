AdjSpec_curve <- function(diseaseData,
                         controlData,
                         userFormula,
                         fixSS,
                         whichSE = "sample",
                         nbootstrap = 100,
                         CI_alpha = 0.95,
                         logit_CI = TRUE,
                         verbose = TRUE) {

    rec_pEst = rec_pEst_curve = rep(0, length(fixSS))
    rec_pSE <- rep(0, length(fixSS))
    for(i in 1:length(fixSS)) {
        EstSE <- GetEst_SE(diseaseData = diseaseData,
                           controlData = controlData,
                           userFormula = userFormula,
                           fixSS = fixSS[i],
                           whichSE = whichSE,
                           nbootstrap = 100,
                           verbose = verbose)
        rec_pEst[i] <- EstSE$Estimates[length(EstSE$Estimates)]
        rec_pSE[i] <- EstSE$SE[length(EstSE$SE)]

        rqfit <- rq(userFormula, tau = -1, data = diseaseData)
        all_Z <- rbind(EstSE$all_Z_D, EstSE$all_Z_C)

        estparam2 <- rbind(rqfit$sol[c(4:nrow(rqfit$sol)), ])
        ctrl_threshold_q_all <- as.matrix(cbind(1, EstSE$all_Z_C)) %*% estparam2
        rho_q_all <- colSums(EstSE$M0 <= ctrl_threshold_q_all)/length(EstSE$M0)
        out_naive <- MonoRespect(tau = rqfit$sol[1,], orig_measure = rho_q_all, startTau = 0.5,
                                 interestTau = 1-fixSS[i])
        rec_pEst_curve[i] <- out_naive$adjrho
    }

    AllCI <- GetAllCI(pEst = rec_pEst_curve,
                      pSE = rec_pSE,
                      alpha = CI_alpha,
                      logit_CI = logit_CI)
    colnames(AllCI) <- fixSS

    return(list(Estimate = rec_pEst_curve,
                SE = rec_pSE,
                ConfidenceInterval = AllCI))
}
