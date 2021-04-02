sscaROC <- function(diseaseData,
                  controlData,
                  userFormula,
                  target_covariates,
                  control_sensitivity = NULL,
                  control_specificity = NULL,
                  mono_resp_method = "ROC",
                  whichSE = "sample",
                  global_ROC_controlled_by = "sensitivity",
                  nbootstrap = 100,
                  CI_alpha = 0.95,
                  logit_CI = TRUE,
                  verbose = TRUE) {

    if (!is.null(control_sensitivity) & !is.null(control_specificity)) {

        stop("Sensitivity and specificity cannot be contolled at the same!")

    } else if (is.null(control_sensitivity) & is.null(control_specificity)) {

        myROC <- ssGetROC(diseaseData = diseaseData,
                    controlData = controlData,
                    userFormula = userFormula,
                    target_covariates = target_covariates,
                    mono_resp_method = mono_resp_method,
                    global_ROC_controlled_by = global_ROC_controlled_by,
                    verbose = verbose)

    } else if (!is.null(control_sensitivity)) {

        if(!is.numeric(control_sensitivity)) {
            control_sensitivity <- as.numeric(control_sensitivity)
        }
        if(sum(control_sensitivity>1)>0 | sum(control_sensitivity<0)>0) {
            stop("control_sensitivity should be between 0 and 1!")
        }

        myROC <- ssGetSpecificity(diseaseData = diseaseData,
                                controlData = controlData,
                                userFormula = userFormula,
                                target_covariates = target_covariates,
                                control_sensitivity = control_sensitivity,
                                whichSE = whichSE,
                                mono_resp_method = mono_resp_method,
                                nbootstrap = nbootstrap,
                                CI_alpha = CI_alpha,
                                logit_CI = logit_CI,
                                verbose = verbose)

    } else if (!is.null(control_specificity)) {

        if(!is.numeric(control_specificity)) {
            control_sensitivity <- as.numeric(control_sensitivity)
        }
        if(sum(control_specificity>1)>0 | sum(control_specificity<0)>0) {
            stop("control_specificity should be between 0 and 1!")
        }

        myROC <- ssGetSensitivity(diseaseData = diseaseData,
                                controlData = controlData,
                                userFormula = userFormula,
                                target_covariates = target_covariates,
                                control_specificity = control_specificity,
                                whichSE = whichSE,
                                mono_resp_method = mono_resp_method,
                                nbootstrap = nbootstrap,
                                CI_alpha = CI_alpha,
                                logit_CI = logit_CI,
                                verbose = verbose)
    }

    return(myROC)
}
