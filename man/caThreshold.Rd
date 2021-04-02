\name{caThreshold}
\alias{caThreshold}
\title{
Calculate covariate-adjusted threshold.
}
\description{
This function is used to calculate covariate-adjusted threshold(s)
at controlled sensitivity levels or specificity levels.
}
\usage{
caThreshold(userFormula, new_covariates, diseaseData = NULL,
controlData = NULL, control_sensitivity = NULL,
control_specificity = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{userFormula}{A character string to represent the function for covariate     adjustment. For example, let Y denote biomarker, Z1 and Z2 denote two covariates. Then userFormula = "Y ~ Z1 + Z2".}
  \item{new_covariates}{A data frame containing covariates for new data. For example, if my userFormula is "Y ~ Z1 + Z2", new_covariates could be data.frame(Z1 = rnorm(100), Z2 = rnorm(100)).}
  \item{diseaseData}{Data from patients including dependent (biomarker) and independent (covariates) variables.}
  \item{controlData}{Data from controls including dependent (biomarker) and independent (covariates) variables.}
  \item{control_sensitivity}{The level(s) of sensitivity to be controlled at. Could be a scalar (e.g. 0.7) or a numeric vector (e.g. c(0.7, 0.8, 0.9)).}
  \item{control_specificity}{The level(s) of specificity to be controlled at. Could be a scalar (e.g. 0.7) or a numeric vector (e.g. c(0.7, 0.8, 0.9)).}
}
\value{
A vector of covariate-adjusted threshold for all subjects if a scalar sensitivity/specificity is given. A data matrix of covariate-adjusted thresholds for all subjects if a vector of sensitivity/specificity is given.
}
\author{
Ziyi Li <ziyi.li@emory.edu>
}
\examples{
n1 = n0 = 500

## generate data
Z_D <- rbinom(n1, size = 1, prob = 0.3)
Z_C <- rbinom(n0, size = 1, prob = 0.7)

Y_C_Z0 <- rnorm(n0, 0.1, 1)
Y_D_Z0 <- rnorm(n1, 1.1, 1)
Y_C_Z1 <- rnorm(n0, 0.2, 1)
Y_D_Z1 <- rnorm(n1, 0.9, 1)

M0 <- Y_C_Z0 * (Z_C == 0) + Y_C_Z1 * (Z_C == 1)
M1 <- Y_D_Z0 * (Z_D == 0) + Y_D_Z1 * (Z_D == 1)

diseaseData <- data.frame(M = M1, Z = Z_D)
controlData <- data.frame(M = M0, Z = Z_C)
userFormula = "M~Z"

### generate new covariates
new_covariates <- data.frame(Z = rbinom(20, size = 1, prob = 0.5))

### calculate covariate-adjusted thresholds at controlled
### sensitivity level 0.7, 0.8, 0.9
caThreshold(userFormula, new_covariates,
            diseaseData = diseaseData,
            controlData = NULL,
            control_sensitivity = c(0.7,0.8,0.9),
            control_specificity = NULL)

### calculate covariate-adjusted thresholds at controlled
### sensitivity level 0.7
caThreshold(userFormula,new_covariates,
            diseaseData = diseaseData,
            controlData = NULL,
            control_sensitivity = 0.7,
            control_specificity = NULL)

### calculate covariate-adjusted thresholds at controlled
### specificity level 0.7, 0.8, 0.9
caThreshold(userFormula,new_covariates,
            diseaseData = NULL,
            controlData = controlData,
            control_sensitivity = NULL,
            control_specificity = c(0.7,0.8,0.9))

### calculate covariate-adjusted thresholds at controlled
### specificity level 0.7
caThreshold(userFormula,new_covariates,
            diseaseData = NULL,
            controlData = controlData,
            control_sensitivity = NULL,
            control_specificity = 0.7)
}
\keyword{ threshold }% use one of  RShowDoc("KEYWORDS")
\keyword{ sensitivity }% __ONLY ONE__ keyword per line
\keyword{ specificity }
