plot_caROC <- function(myROC,...){

    message("Plotting the ROC curve...")

    if (myROC$mono_adj == "ROC") {
        linetype = "l"
    } else {
        linetype = "s"
    }

    plot(1-myROC$specificity, myROC$sensitivity, xlab = "1 - Specificity",
         ylab = "Sensitivity", type = linetype, ...)
}

