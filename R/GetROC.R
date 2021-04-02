GetROC <- function(diseaseData,
                   controlData,
                   userFormula,
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
            myROC <- AdjSpec_AllPoints_none(diseaseData = diseaseData,
                                            controlData = controlData,
                                            userFormula = userFormula)
        } else if (mono_resp_method == "ROC") {
            myROC <- AdjSpec_AllPoints_curve(diseaseData = diseaseData,
                                             controlData = controlData,
                                             userFormula = userFormula)
        } else {
            stop("mono_resp_method need to be among the following: ROC/none!")
        }

        if (min(myROC$sensitivity)!=0 | max(myROC$specificity)!=1) {
            myROC <- list(sensitivity = c(myROC$sensitivity, 0),
                          specificity = c(myROC$specificity, 1),
                          mono_adj = myROC$mono_adj)
        }
        if (max(myROC$sensitivity)!=1 | min(myROC$specificity)!=0) {
            myROC <- list(sensitivity = c(1, myROC$sensitivity),
                          specificity = c(0, myROC$specificity),
                          mono_adj = myROC$mono_adj)
        }

    } else if (global_ROC_controlled_by == "specificity") {
        if(verbose) {
            message("Global ROC by controlling specificity.")
        }

        revData <- reverseData(diseaseData = diseaseData,
                               controlData = controlData,
                               userFormula = userFormula)

        if (mono_resp_method == "none") {
            tmpROC <- AdjSpec_AllPoints_none(diseaseData = revData$diseaseData,
                                            controlData = revData$controlData,
                                            userFormula = userFormula)
        } else if (mono_resp_method == "ROC") {
            tmpROC <- AdjSpec_AllPoints_curve(diseaseData = revData$diseaseData,
                                             controlData = revData$controlData,
                                             userFormula = userFormula)
        } else {
            stop("mono_resp_method need to be among the following: ROC/none!")
        }

        myROC <- list(sensitivity = tmpROC$specificity,
                      specificity = tmpROC$sensitivity,
                      mono_adj = tmpROC$mono_adj)

        if (min(myROC$sensitivity)!=0 | max(myROC$specificity)!=1) {
            myROC <- list(sensitivity = c(0, myROC$sensitivity),
                          specificity = c(1, myROC$specificity),
                          mono_adj = myROC$mono_adj)
        }
        if (max(myROC$sensitivity)!=1 | min(myROC$specificity)!=0) {
            myROC <- list(sensitivity = c(myROC$sensitivity, 1),
                          specificity = c(myROC$specificity, 0),
                          mono_adj = myROC$mono_adj)
        }

    } else {
        stop("Global ROC can only be controlled by sensitivity or specificity!")
    }

    return(myROC)
}
