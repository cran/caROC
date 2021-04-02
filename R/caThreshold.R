caThreshold <- function(userFormula,
                        new_covariates,
                        diseaseData = NULL,
                        controlData = NULL,
                        control_sensitivity = NULL,
                        control_specificity = NULL) {

    if (is.null(control_sensitivity) & is.null(control_specificity)) {
        stop("Give a sensitivity/specificity to control at!")
    } else if (!is.null(control_sensitivity) & !is.null(control_specificity)) {
        stop("Sensitivity and specificity cannot be controlled at the same time!")
    } else if (!is.null(control_sensitivity)) {

        if(is.null(diseaseData)) {
            stop("I need data from patients (diseaseData)!")
        } else {
            cathres <- caThreshold_at_sensitivity(diseaseData = diseaseData,
                                       userFormula = userFormula,
                                       control_sensitivity = control_sensitivity,
                                       new_covariates = new_covariates)
        }

    } else if (!is.null(control_specificity)) {

        if(is.null(controlData)) {
            stop("I need data from normal controls (controlData)!")
        } else {
            cathres <- caThreshold_at_specificity(controlData = controlData,
                                       userFormula = userFormula,
                                       control_specificity = control_specificity,
                                       new_covariates = new_covariates)
        }
    }

    return(cathres)
}
