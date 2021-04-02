ssGetSensitivity <- function(diseaseData,
                           controlData,
                           userFormula,
                           target_covariates,
                           control_specificity,
                           whichSE,
                           mono_resp_method,
                           nbootstrap,
                           CI_alpha,
                           logit_CI,
                           verbose) {
    
    if(verbose) {
        message("Calculate covariate-adjusted sensitivity at controlled specificity.")
    }
    
    revData <- reverseData_L(diseaseData = diseaseData,
                           controlData = controlData,
                           userFormula = userFormula)
    
    if (mono_resp_method == "none") {
        myROC <- AdjSpec_none_L(diseaseData = revData$diseaseData,
                              controlData = revData$controlData,
                              userFormula = userFormula,
                              target_covariates = target_covariates,
                              fixSS = control_specificity,
                              whichSE = whichSE,
                              nbootstrap = nbootstrap,
                              CI_alpha = CI_alpha,
                              logit_CI = logit_CI,
                              verbose = verbose)
    } else if (mono_resp_method == "ROC") {
        myROC <- AdjSpec_curve_L(diseaseData = revData$diseaseData,
                               controlData = revData$controlData,
                               userFormula = userFormula,
                               target_covariates = target_covariates,
                               fixSS = control_specificity,
                               whichSE = whichSE,
                               nbootstrap = nbootstrap,
                               CI_alpha = CI_alpha,
                               logit_CI = logit_CI,
                               verbose = verbose)
    } else {
        stop("mono_resp_method need to be among the following: ROC/none!")
    }
    
    return(myROC)
}
