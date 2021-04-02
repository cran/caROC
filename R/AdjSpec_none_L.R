AdjSpec_none_L <- function(diseaseData,
                         controlData,
                         userFormula,
                         target_covariates,
                         fixSS,
                         whichSE = "sample",
                         nbootstrap = 100,
                         CI_alpha = 0.95,
                         logit_CI = TRUE,
                         verbose = TRUE) {
    
    if(!is.matrix(target_covariates)) {
        rec_pEst <- rep(0, length(fixSS))
        rec_pSE <- rep(0, length(fixSS))
        for(i in 1:length(fixSS)) {
            EstSE <- GetEst_SE_L(diseaseData = diseaseData,
                                 controlData = controlData,
                                 userFormula = userFormula,
                                 target_covariates = target_covariates,
                                 fixSS = fixSS[i],
                                 whichSE = whichSE,
                                 nbootstrap = nbootstrap,
                                 verbose = verbose)
            rec_pEst[i] <- EstSE$Estimates[length(EstSE$Estimates)]
            rec_pSE[i] <- EstSE$SE[length(EstSE$SE)]
        }
        
        AllCI <- GetAllCI(pEst = rec_pEst,
                          pSE = rec_pSE,
                          alpha = CI_alpha,
                          logit_CI = logit_CI)
        colnames(AllCI) <- fixSS
    } else if (is.matrix(target_covariates)) {
        rec_pEst = rec_pSE = AllCI = list()
        for (jj in 1:nrow(target_covariates)) {
            one_pEst = one_pSE = rep(0, length(fixSS))
            for(i in 1:length(fixSS)) {
                EstSE <- GetEst_SE_L(diseaseData = diseaseData,
                                     controlData = controlData,
                                     userFormula = userFormula,
                                     target_covariates = target_covariates,
                                     fixSS = fixSS[i],
                                     whichSE = whichSE,
                                     nbootstrap = nbootstrap,
                                     verbose = verbose)
                one_pEst[i] <- EstSE$Estimates[length(EstSE$Estimates)]
                one_pSE[i] <- EstSE$SE[length(EstSE$SE)]
            }
            OneCI <- GetAllCI(pEst = one_pEst,
                              pSE = one_pSE,
                              alpha = CI_alpha,
                              logit_CI = logit_CI)
            colnames(OneCI) <- fixSS
            
            rec_pEst[[jj]] <- one_pEst
            rec_pSE[[jj]] <- one_pSE
            AllCI[[jj]] <- OneCI
        }
        
        names(rec_pSE) = names(rec_pEst) = names(AllCI) = apply(target_covariates, 1, paste, collapse = "_")
    }
    
    
    return(list(Estimate = rec_pEst,
                SE = rec_pSE,
                ConfidenceInterval = AllCI))
}
