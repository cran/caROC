GetAllCI <- function(pEst, pSE,
                     alpha = 0.95,
                     logit_CI = TRUE) {
    
    sai <- qnorm(1 - (1-alpha)/2)
    allCI <- matrix(0, 2, length(pEst))
    if (logit_CI) {
        allCI[1, ] <- revlogit(logit(pEst) - sai * pSE/pEst/(1-pEst))
        allCI[2, ] <- revlogit(logit(pEst) + sai * pSE/pEst/(1-pEst))
        rownames(allCI) <- c("LogitCI_Lower", "LogitCI_Upper")
    } else {
        allCI[1, ] <- pEst - sai * pSE
        allCI[2, ] <- pEst + sai * pSE
        rownames(allCI) <- c("CI_Lower", "CI_Upper")
    }
    
    return(allCI)
}

logit <- function(p) {
    return(log(p/(1-p)))
}

revlogit <- function(t) {
    if (length(t)>1) {
        if(any(t>300) | any(is.na(t))) {
            t[t>300] <- 300
        }
    } else {
        if(t>300 | any(is.na(t))) {
            t <- 300
        }
    }
    
    return(exp(t)/(1+exp(t)))
}
