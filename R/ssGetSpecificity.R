ssGetSpecificity <- function(diseaseData,
                           controlData,
                           userFormula,
                           target_covariates = target_covariates,
                           control_sensitivity,
                           whichSE,
                           mono_resp_method,
                           nbootstrap,
                           CI_alpha,
                           logit_CI,
                           verbose) {
    
    if(verbose) {
        message("Calculate covariate-adjusted specificity at controlled sensitivity.")
    }
    
    if (mono_resp_method == "none") {
        myROC <- AdjSpec_none_L(diseaseData = diseaseData,
                              controlData = controlData,
                              userFormula = userFormula,
                              target_covariates = target_covariates,
                              fixSS = control_sensitivity,
                              whichSE = whichSE,
                              nbootstrap = nbootstrap,
                              CI_alpha = CI_alpha,
                              logit_CI = logit_CI,
                              verbose = verbose)
    } else if (mono_resp_method == "ROC") {
        myROC <- AdjSpec_curve_L(diseaseData = diseaseData,
                               controlData = controlData,
                               userFormula = userFormula,
                               target_covariates = target_covariates,
                               fixSS = control_sensitivity,
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
