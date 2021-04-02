ParseuserFormula <- function(userFormula) {
    formula2 <- gsub(" ", "", userFormula)
    tt1 <- unlist(strsplit(formula2, split = "~"))
    allvar <- unlist(strsplit(tt1[2], split = "\\+"))
    return(list(response = tt1[1],
                allvar = allvar))
}
