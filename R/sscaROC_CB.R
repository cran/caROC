sscaROC_CB <- function(diseaseData,
                     controlData,
                     userFormula,
                     mono_resp_method = "none",
                     target_covariates,
                     global_ROC_controlled_by = "sensitivity",
                     CB_alpha = 0.95,
                     logit_CB = FALSE,
                     nbootstrap = 100,
                     nbin = 100,
                     verbose = FALSE) {

    if (global_ROC_controlled_by == "sensitivity") {

        message("Global ROC by controlling sensitivity.")

    } else if (global_ROC_controlled_by == "specificity") {

        message("Global ROC by controlling specificity.")

        revData <- reverseData_L(diseaseData = diseaseData,
                               controlData = controlData,
                               userFormula = userFormula)
        diseaseData <- revData$diseaseData
        controlData <- revData$controlData

    } else {
        stop("Global ROC can only be controlled by sensitivity or specificity!")
    }

    origROC <- sscaROC(diseaseData = diseaseData,
                     controlData = controlData,
                     userFormula = userFormula,
                     target_covariates = target_covariates,
                     mono_resp_method = mono_resp_method,
                     verbose = verbose)

    sens_vec <- seq(0, 1, 1/nbin)
    supmat <- matrix(0, nbootstrap, length(sens_vec))
    pb = txtProgressBar(min = 0, max = nbootstrap, style = 3)
    for(nidx in 1:nbootstrap) {
        setTxtProgressBar(pb, nidx)

        oneindx1 <- sample(1:nrow(controlData), nrow(controlData), replace = TRUE)
        oneindx2 <- sample(1:nrow(diseaseData), nrow(diseaseData), replace = TRUE)

        controlData_bt <- controlData[oneindx1, ]
        diseaseData_bt <- diseaseData[oneindx2, ]

        btROC <- sscaROC(diseaseData = diseaseData_bt,
                       controlData = controlData_bt,
                       userFormula = userFormula,
                       target_covariates = target_covariates,
                       mono_resp_method = mono_resp_method,
                       verbose = verbose)

        outres <- GetAbsSup(origROC, btROC, sens_vec, logit_CB)
        supmat[nidx,] <- outres$alldiff

    }
    close(pb)
    finalorigY <- outres$finalorigY


    CB_upper = CB_lower = rep(NA, length(sens_vec))
    if (logit_CB) {
        deltaterm <- apply(supmat, 2, function(x) quantile(x, CB_alpha))
        deltaterm[deltaterm>2] <- max(deltaterm[deltaterm != Inf])
        CB_upper <- revlogit(finalorigY + deltaterm)
        CB_lower <- revlogit(finalorigY - deltaterm)
    } else {
        deltaterm <- apply(supmat, 2, function(x) quantile(x, CB_alpha))
        CB_upper <- finalorigY + deltaterm
        CB_lower <- finalorigY - deltaterm
        CB_upper[CB_upper > 1] = 1
        CB_upper[CB_upper < 0] = 0
        CB_lower[CB_lower > 1] = 1
        CB_lower[CB_lower < 0] = 0

    }

    if (global_ROC_controlled_by == "sensitivity") {
        return(list(Sensitivity = sens_vec,
                    Specificity_upper = CB_upper,
                    Specificity_lower = CB_lower,
                    global_ROC_controlled_by = global_ROC_controlled_by))
    } else {
        return(list(Specificity = sens_vec,
                    Sensitivity_upper = CB_upper,
                    Sensitivity_lower = CB_lower,
                    global_ROC_controlled_by = global_ROC_controlled_by))
    }

}
