plot_sscaROC <- function(myROC,...){

    message("Plotting the ROC curve...")

    if (myROC$mono_adj == "ROC") {
        if(!is.list(myROC$sensitivity) & !is.list(myROC$specificity)) {
            plot(1-myROC$specificity, myROC$sensitivity, xlab = "1 - Specificity",
                 ylab = "Sensitivity", type = "l", ...)
        } else {
            for(j in 1:length(myROC$sensitivity)) {
                newcol <- RColorBrewer::brewer.pal(max(length(myROC$specificity), 3), "Set1")
                if(j == 1) {
                    plot(1-myROC$specificity[[j]], myROC$sensitivity[[j]], xlab = "1 - Specificity",
                         ylab = "Sensitivity", col = newcol[j], type = "l", lty = j, ...)
                } else {
                    points(1-myROC$specificity[[j]], myROC$sensitivity[[j]], xlab = "1 - Specificity",
                         ylab = "Sensitivity", col = newcol[j], type = "l", lty = j, ...)
                }
                if(length(myROC$sensitivity)>1) {
                    legend("bottomright", legend = names(myROC$sensitivity), col = newcol, lty = 1:length(myROC$sensitivity), ...)
                }
            }
        }
    } else {
        if(!is.matrix(myROC$sensitivity) & !is.matrix(myROC$specificity)) {
            plot(1-myROC$specificity, myROC$sensitivity, xlab = "1 - Specificity",
                 ylab = "Sensitivity", type = "l", ...)
        } else if(is.matrix(myROC$sensitivity)) {
            for(i in 1:nrow(myROC$sensitivity)) {
                newcol <- RColorBrewer::brewer.pal(max(nrow(myROC$sensitivity), 3), "Set1")
                if(i == 1) {
                    plot(1-myROC$specificity, myROC$sensitivity[i,], xlab = "1 - Specificity",
                         ylab = "Sensitivity", col = newcol[i], type = "l", lty = i,...)
                } else {
                    points(1-myROC$specificity, myROC$sensitivity[i,], xlab = "1 - Specificity",
                         ylab = "Sensitivity", col = newcol[i], type = "l", lty = i,...)
                }
                if(nrow(myROC$sensitivity) > 1) {
                    legend("bottomright", legend = rownames(myROC$sensitivity), col = newcol, lty = 1:nrow(myROC$sensitivity), ...)
                }
            }
        } else if(is.matrix(myROC$specificity)) {
            for(i in 1:nrow(myROC$specificity)) {
                newcol <- RColorBrewer::brewer.pal(max(nrow(myROC$specificity), 3), "Set1")
                if(i == 1) {
                    plot(1-myROC$specificity[i,], myROC$sensitivity, xlab = "1 - Specificity",
                         ylab = "Sensitivity", col = newcol[i], type = "l", lty = i, ...)
                } else {
                    points(1-myROC$specificity[i,], myROC$sensitivity, xlab = "1 - Specificity",
                           ylab = "Sensitivity", col = newcol[i], type = "l", lty = i, ...)
                }
                if(nrow(myROC$specificity) > 1) {
                    legend("bottomright", legend = rownames(myROC$specificity), col = newcol, lty = 1:nrow(myROC$specificity), ...)
                }
            }
        }
    }
}

