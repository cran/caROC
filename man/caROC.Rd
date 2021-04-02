\name{caROC}
\alias{caROC}
\title{
Covariate-adjusted ROC}
\description{
Compute covariate-adjusted specificity at controlled sensitivity level, or covariate-adjusted sensitivity at controlled specificity level, or covariate-adjust receiver operating characteristic curve.
}
\usage{
caROC(diseaseData, controlData, userFormula, control_sensitivity = NULL,
control_specificity = NULL, mono_resp_method = "ROC",
whichSE = "sample", global_ROC_controlled_by = "sensitivity",
nbootstrap = 100, CI_alpha = 0.95, logit_CI = TRUE,
verbose = TRUE)
}
\arguments{
  \item{diseaseData}{
Data from patients including dependent (biomarker) and independent (covariates) variables.
}
  \item{controlData}{
Data from controls including dependent (biomarker) and independent (covariates) variables.
}
  \item{userFormula}{
  A character string to represent the function for covariate adjustment. For example, let Y denote biomarker, Z1 and Z2 denote two covariates. Then userFormula = "Y ~ Z1 + Z2".
  }
  \item{control_sensitivity}{
The level(s) of sensitivity to be controlled at. Could be a scalar (e.g. 0.7) or a numeric vector (e.g. c(0.7, 0.8, 0.9)).
}
  \item{control_specificity}{
The level(s) of specificity to be controlled at. Could be a scalar (e.g. 0.7) or a numeric vector (e.g. c(0.7, 0.8, 0.9)).
}
  \item{mono_resp_method}{
The method used to restore monotonicity of the ROC curve or computed sensitivity/specificity value. It should one from the following: "none", "ROC". "none" is not applying any monotonicity respecting method. "ROC" is to apply ROC-based monotonicity respecting approach. Default value is "ROC".
}
  \item{whichSE}{
The method used to compute standard error. It should be one from the following: "sample", "bootstrap", meaning to calculate the standard error using sample-based approach or bootstrap. Default is "sample".
}
  \item{global_ROC_controlled_by}{
Whether sensitivity/specificity is used to control when computing global ROC. It should one from the following: "sensitivity", "specificity". Default is "sensitivity".
}
  \item{nbootstrap}{
Number of boostrap iterations. Default is 100.
}
  \item{CI_alpha}{
Percentage of confidence interval. Default is 0.95.
}
  \item{logit_CI}{
Whether to apply logit-based confidence interval. Logit-transformed CI has been identified to be more robust near border area.
}
    \item{verbose}{
Whether to print out messages. Default value is true.
}
}
\value{
If control_sensitivity or control_specificity is provided, compute
covariate-adjusted specificity (sensitivity) at controlled sensitivity
(specificity) level.
\item{Estimate}{Covariate-adjusted sensitivity/specificity.}
\item{SE}{Estimated standard error.}
\item{CI}{Estimated confidence intervals.}
If both control_sensitivity and control_specificity are null,
compuate covariate-adjusted ROC curve.
\item{sensitivity}{Estimated sensitivity.}
\item{specificity}{Estimated specificity.}
\item{mono_adj}{Monotonicity adjustment method.}
}
\author{
Ziyi.li <ziyi.li@emory.edu>
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

## calculate covariate-adjusted specificity at
## controlled sensitivity levels (0.2, 0.8, 0.9)
caROC(diseaseData,controlData,userFormula,
      control_sensitivity = c(0.2,0.8, 0.9),
      control_specificity = NULL,mono_resp_method = "ROC",
      whichSE = "bootstrap",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)

## calculate covariate-adjusted sensitivity at
## controlled specificity levels (0.2, 0.8, 0.9)
caROC(diseaseData,controlData,userFormula,
      control_sensitivity = NULL,
      control_specificity = c(0.7,0.8, 0.9),mono_resp_method = "none",
      whichSE = "sample",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)

## calculate the whole covariate-adjusted ROC curve
ROC1 <- caROC(diseaseData,controlData,userFormula = "M~Z",
                 mono_resp_method = "none")
ROC2 <- caROC(diseaseData,controlData,userFormula = "M~Z",
                 mono_resp_method = "ROC")
}
\keyword{ sensitivity }% use one of  RShowDoc("KEYWORDS")
\keyword{ specificity }% __ONLY ONE__ keyword per line
\keyword{ ROC }
