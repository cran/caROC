AdjSpec_none <- function(diseaseData,
                         controlData,
                         userFormula,
                         fixSS,
                         whichSE = "sample",
                         nbootstrap = 100,
                         CI_alpha = 0.95,
                         logit_CI = TRUE,
                         verbose = TRUE) {

    rec_pEst <- rep(0, length(fixSS))
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
    }

    AllCI <- GetAllCI(pEst = rec_pEst,
             pSE = rec_pSE,
             alpha = CI_alpha,
             logit_CI = logit_CI)
    colnames(AllCI) <- fixSS

    return(list(Estimate = rec_pEst,
                SE = rec_pSE,
                ConfidenceInterval = AllCI))
}
