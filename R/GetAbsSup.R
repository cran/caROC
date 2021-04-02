GetAbsSup <- function(origROC, btROC, sens_vec,
                      logit_CB = FALSE) {

    tmpX <- btROC$sensitivity
    tmpY <- btROC$specificity
    if(tmpX[1]>tmpX[length(tmpX)]) {
        btX <- rev(tmpX)
        btY <- rev(tmpY)
    } else {
        btX <- tmpX
        btY <- tmpY
    }

    tmpX <- origROC$sensitivity
    tmpY <- origROC$specificity
    if(tmpX[1]>tmpX[length(tmpX)]) {
        origX <- rev(tmpX)
        origY <- rev(tmpY)
    } else {
        origX <- tmpX
        origY <- tmpY
    }

    if(logit_CB) {
        # btY[btY == 1] <- 0.9999
        # origY[origY == 1] <- 0.9999
        btY <- logit(btY)
        origY <- logit(origY)
    }

    longX <- sort(c(btX, origX))
    longbtY = longorigY = rep(NA, length(longX))
    for(j in 1:length(longX)) {
        ssval <- longX[j]
        longbtY[j] <- btY[findInterval(ssval, c(-Inf, btX))-1]
        longorigY[j] <- origY[findInterval(ssval, c(-Inf, origX))-1]
    }
    longabsdiff <- abs(longbtY - longorigY)
    longX2 <- longX[!is.na(longabsdiff)]
    longabsdiff2 <- na.omit(longabsdiff)

    alldiff = finalorigY = rep(NA, length(sens_vec))
    vv_rec <- rep(1,length(sens_vec)+1)
    for(i in 1:length(sens_vec)) {
        ssval <- sens_vec[i]
        if(findInterval(ssval, c(-Inf, longX2)) == 1) {
            vv_rec[i+1] <- 1
        } else {
            vv_rec[i+1] <- findInterval(ssval, c(-Inf, longX2)) - 1
        }
        alldiff[i] <- max(longabsdiff2[vv_rec[i]:vv_rec[i+1]])
        finalorigY[i] <- longorigY[vv_rec[i+1]]
    }
    return(list(alldiff = alldiff,
                finalorigY = finalorigY))
}
