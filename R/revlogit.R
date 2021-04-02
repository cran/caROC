revlogit <- function(t) {
    if (length(t)>1) {
        if(any(t>300)) {
            t[t>300] <- 300
        }
    } else {
        if(t>300) {
            t <- 300
        }
    }

    return(exp(t)/(1+exp(t)))
}
