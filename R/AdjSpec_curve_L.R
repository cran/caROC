AdjSpec_curve_L <- function(diseaseData,
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
        rec_pEst = rec_pEst_curve = rep(0, length(fixSS))
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

            rqfit <- quantreg::rq(userFormula, tau = -1, data = diseaseData)
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
    } else if (is.matrix(target_covariates)) {
        rec_pEst = rec_pEst_curve =rec_pSE = AllCI = list()
        for (jj in 1:nrow(target_covariates)) {
            one_pEst = one_pEst_curve = one_pSE = rep(0, length(fixSS))
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

                rqfit <- quantreg::rq(userFormula, tau = -1, data = diseaseData)
                all_Z <- rbind(EstSE$all_Z_D, EstSE$all_Z_C)

                estparam2 <- rbind(rqfit$sol[c(4:nrow(rqfit$sol)), ])
                ctrl_threshold_q_all <- as.matrix(cbind(1, EstSE$all_Z_C)) %*% estparam2
                rho_q_all <- colSums(EstSE$M0 <= ctrl_threshold_q_all)/length(EstSE$M0)
                out_naive <- MonoRespect(tau = rqfit$sol[1,], orig_measure = rho_q_all, startTau = 0.5,
                                         interestTau = 1-fixSS[i])
                one_pEst_curve[i] <- out_naive$adjrho
            }
            OneCI <- GetAllCI(pEst = one_pEst_curve,
                              pSE = one_pSE,
                              alpha = CI_alpha,
                              logit_CI = logit_CI)
            colnames(OneCI) <- fixSS

            rec_pEst_curve[[jj]] <- one_pEst_curve
            rec_pSE[[jj]] <- one_pSE
            AllCI[[jj]] <- OneCI
        }

        names(rec_pSE) = names(rec_pEst_curve) = names(AllCI) = apply(target_covariates, 1, paste, collapse = "_")
    }


    return(list(Estimate = rec_pEst_curve,
                SE = rec_pSE,
                ConfidenceInterval = AllCI))
}
