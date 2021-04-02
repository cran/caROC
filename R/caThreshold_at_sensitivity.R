caThreshold_at_sensitivity <- function(diseaseData,
                           userFormula,
                           control_sensitivity,
                           new_covariates) {

    if(!is.numeric(control_sensitivity)) {
        control_sensitivity <- as.numeric(control_sensitivity)
    }
    if(sum(control_sensitivity>1)>0 | sum(control_sensitivity<0)>0) {
        stop("control_sensitivity should be between 0 and 1!")
    }

    rqfit <- NULL
    if (length(control_sensitivity) == 1) {
        expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", 1-control_sensitivity, ", data = diseaseData)")
        eval(parse(text = expr1))
        new_thres <- predict.rq(rqfit, newdata = new_covariates)

    } else {
        new_thres <- c()
        for(j in 1:length(control_sensitivity)) {
            expr1 <- paste0("rqfit <- rq(", userFormula, ", tau = ", 1-control_sensitivity[j], ", data = diseaseData)")
            eval(parse(text = expr1))
            tmp <- predict.rq(rqfit, newdata = new_covariates)
            new_thres <- cbind(new_thres, tmp)
        }
        colnames(new_thres) = paste0("control_sens=", control_sensitivity)
    }

    return(new_thres)

}
