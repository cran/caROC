ssGetROC <- function(diseaseData,
                   controlData,
                   userFormula,
                   target_covariates,
                   mono_resp_method,
                   global_ROC_controlled_by,
                   verbose) {
    
    if(verbose) {
        message("Calculate sensitivity and specificity for all thresholds along ROC.")
    }
    
    if (global_ROC_controlled_by == "sensitivity") {
        if(verbose) {
            message("Global ROC by controlling sensitivity.")
        }
        
        if (mono_resp_method == "none") {
            myROC <- AdjSpec_AllPoints_none_L(diseaseData = diseaseData,
                                            controlData = controlData,
                                            userFormula = userFormula,
                                            target_covariates = target_covariates)
        } else if (mono_resp_method == "ROC") {
            myROC <- AdjSpec_AllPoints_curve_L(diseaseData = diseaseData,
                                             controlData = controlData,
                                             userFormula = userFormula,
                                             target_covariates = target_covariates)
        } else {
            stop("mono_resp_method need to be among the following: ROC/none!")
        }
        
        
    } else if (global_ROC_controlled_by == "specificity") {
        if(verbose) {
            message("Global ROC by controlling specificity.")
        }
        
        revData <- reverseData_L(diseaseData = diseaseData,
                               controlData = controlData,
                               userFormula = userFormula)
        
        if (mono_resp_method == "none") {
            tmpROC <- AdjSpec_AllPoints_none_L(diseaseData = revData$diseaseData,
                                             controlData = revData$controlData,
                                             userFormula = userFormula,
                                             target_covariates = target_covariates)
        } else if (mono_resp_method == "ROC") {
            tmpROC <- AdjSpec_AllPoints_curve_L(diseaseData = revData$diseaseData,
                                              controlData = revData$controlData,
                                              userFormula = userFormula,
                                              target_covariates = target_covariates)
        } else {
            stop("mono_resp_method need to be among the following: ROC/none!")
        }
        
        myROC <- list(sensitivity = tmpROC$specificity,
                      specificity = tmpROC$sensitivity,
                      mono_adj = tmpROC$mono_adj)
        
    } else {
        stop("Global ROC can only be controlled by sensitivity or specificity!")
    }
    
    return(myROC)
}
