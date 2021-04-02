MonoRespect <- function(tau, orig_measure, startTau = 0.5,
                        interestTau = NULL) {

    if (is.null(startTau)) {
        startTau = tau[round(median(1:length(tau)))]
    }
    leftRec <- data.frame(tau = rep(NA, length(tau)),
                          new_meas = rep(NA, length(tau)))
    rightRec <- data.frame(tau = rep(NA, length(tau)),
                           new_meas = rep(NA, length(tau)))

    myi = 1
    leftRec$tau[myi] = startTau
    leftRec$new_meas[myi] = orig_measure[which.min(abs(startTau-tau))]
    new_idx = 1
    while (!is.na(new_idx)) {
        tmp <- which(orig_measure <= leftRec$new_meas[myi] & tau < leftRec$tau[myi])
        if (length(tmp) >= 1) {
            new_idx <- max(tmp)
            myi <- myi+1
            leftRec$tau[myi] <- tau[new_idx]
            leftRec$new_meas[myi] <- orig_measure[new_idx]
        } else {
            new_idx = NA
        }
    }
    leftRec <- na.omit(leftRec)

    myi = 1
    rightRec$tau[myi] = startTau
    rightRec$new_meas[myi] = orig_measure[which.min(abs(startTau-tau))]
    new_idx = 1
    while (!is.na(new_idx)) {
        tmp <- which(orig_measure >= rightRec$new_meas[myi] & tau > rightRec$tau[myi])
        if (length(tmp) >= 1) {
            new_idx <- min(tmp)
            myi <- myi+1
            rightRec$tau[myi] <- tau[new_idx]
            rightRec$new_meas[myi] <- orig_measure[new_idx]
        } else {
            new_idx = NA
        }
    }
    rightRec <- na.omit(rightRec)

    rev_leftRec <- cbind(rev(leftRec$tau), rev(leftRec$new_meas))
    colnames(rev_leftRec) <- colnames(leftRec)
    monoRes <- rbind(rev_leftRec, rightRec[-1,])

    if (!is.null(interestTau)) {
        if (interestTau %in% monoRes$tau) {
            indx <- which(monoRes$tau == interestTau)
            thisrho <- monoRes$new_meas[indx]
        } else {
            if (length(monoRes$tau[monoRes$tau < interestTau])>0) {
                indx1 <- which(monoRes$tau == max(monoRes$tau[monoRes$tau < interestTau]))
            } else {
                indx1 <- 1
            }
            if(length(monoRes$tau[monoRes$tau > interestTau])>0) {
                indx2 <- which(monoRes$tau == min(monoRes$tau[monoRes$tau > interestTau]))
            } else {
                indx2 <- length(monoRes$tau)
            }

            tau_left <- monoRes$tau[indx1]
            tau_right <- monoRes$tau[indx2]
            rho_left <- monoRes$new_meas[indx1]
            rho_right <- monoRes$new_meas[indx2]
            denomi <- tau_right - tau_left
            mydenomi <- ifelse(denomi == 0, 1, denomi)
            thisrho <- rho_left + (rho_right - rho_left)/mydenomi * (interestTau - tau_left)
        }
        return(list(monoRes = monoRes,
                    adjrho = thisrho))
    } else {
        return(monoRes)
    }
}
