GetSensitivity <- function(diseaseData,
                           controlData,
                           userFormula,
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

    revData <- reverseData(diseaseData = diseaseData,
                           controlData = controlData,
                           userFormula = userFormula)

    if (mono_resp_method == "none") {
        myROC <- AdjSpec_none(diseaseData = revData$diseaseData,
                              controlData = revData$controlData,
                              userFormula = userFormula,
                              fixSS = control_specificity,
                              whichSE = whichSE,
                              nbootstrap = nbootstrap,
                              CI_alpha = CI_alpha,
                              logit_CI = logit_CI,
                              verbose = verbose)
    } else if (mono_resp_method == "ROC") {
        myROC <- AdjSpec_curve(diseaseData = revData$diseaseData,
                               controlData = revData$controlData,
                               userFormula = userFormula,
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
