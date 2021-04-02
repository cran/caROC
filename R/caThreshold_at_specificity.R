caThreshold_at_specificity <- function(controlData,
                           userFormula,
                           control_specificity,
                           new_covariates) {

    if(!is.numeric(control_specificity)) {
        control_specificity <- as.numeric(control_specificity)
    }
    if(sum(control_specificity>1)>0 | sum(control_specificity<0)>0) {
        stop("control_specificity should be between 0 and 1!")
    }

    rqfit <- NULL
    if (length(control_specificity) == 1) {
        expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", control_specificity, ", data = controlData)")
        eval(parse(text = expr1))
        new_thres <- predict.rq(rqfit, newdata = new_covariates)

    } else {
        new_thres <- c()
        for(j in 1:length(control_specificity)) {
            expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", control_specificity[j], ", data = controlData)")
            eval(parse(text = expr1))
            tmp <- predict.rq(rqfit, newdata = new_covariates)
            new_thres <- cbind(new_thres, tmp)
        }
        colnames(new_thres) = paste0("control_spec=", control_specificity)
    }

    return(new_thres)

}
