reverseData_L <- function(diseaseData,
                        controlData,
                        userFormula) {

    outform <- ParseuserFormula_L(userFormula = userFormula)
    revdiseaseData <- controlData
    revcontrolData <- diseaseData

    expr1 <- parse(text = paste0("revcontrolData$",outform[[1]], "<- -diseaseData$",outform[[1]]))
    expr2 <- parse(text = paste0("revdiseaseData$",outform[[1]], "<- -controlData$",outform[[1]]))
    eval(expr1)
    eval(expr2)

    return(list(diseaseData = revdiseaseData,
                controlData = revcontrolData))
}
