plot_sscaROC_CB <- function(myROC_CB, add = TRUE,...) {
    if (add) {
        if (myROC_CB$global_ROC_controlled_by == "sensitivity") {
            points(1-myROC_CB$Specificity_upper, myROC_CB$Sensitivity, type = "s", ...)
            points(1-myROC_CB$Specificity_lower, myROC_CB$Sensitivity, type = "s", ...)
        } else if (myROC_CB$global_ROC_controlled_by == "specificity") {
            points(1-myROC_CB$Specificity, myROC_CB$Sensitivity_upper, type = "s", ...)
            points(1-myROC_CB$Specificity, myROC_CB$Sensitivity_lower, type = "s", ...)
        }
    } else {
        if (myROC_CB$global_ROC_controlled_by == "sensitivity") {
            plot(1-myROC_CB$Specificity_upper, myROC_CB$Sensitivity, type = "s", ...)
            points(1-myROC_CB$Specificity_lower, myROC_CB$Sensitivity, type = "s", ...)
        } else if (myROC_CB$global_ROC_controlled_by == "specificity") {
            plot(1-myROC_CB$Specificity, myROC_CB$Sensitivity_upper, type = "s", ...)
            points(1-myROC_CB$Specificity, myROC_CB$Sensitivity_lower, type = "s", ...)
        }
    }
}
