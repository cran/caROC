\name{sscaROC_CB}
\alias{sscaROC_CB}
\title{
Get confidence band for covariate-adjusted ROC curve for specified sub-population.
}
\description{
Use this function to compute the confidence band for covariate-adjusted ROC curve, with or without monotonicity respecting methods for sub-population.
}
\usage{
sscaROC_CB(diseaseData, controlData, userFormula, mono_resp_method = "none",
target_covariates, global_ROC_controlled_by = "sensitivity", CB_alpha = 0.95,
logit_CB = FALSE, nbootstrap = 100, nbin = 100, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{mono_resp_method}{
The method used to restore monotonicity of the ROC curve or computed sensitivity/specificity value. It should one from the following: "none", "ROC". "none" is not applying any monotonicity respecting method. "ROC" is to apply ROC-based monotonicity respecting approach. Default value is "ROC".
}
  \item{target_covariates}{
Covariates of the interested sub-population. It could be a vector, e.g. c(1, 0.5, 0.8), or a matrix, e.g. target_covariates = matrix(c(1, 0.7, 0.9, 1, 0.8, 0.8), 2, 3, byrow = TRUE)
}
  \item{global_ROC_controlled_by}{
Whether sensitivity/specificity is used to control when computing global ROC. It should one from the following: "sensitivity", "specificity". Default is "sensitivity".
}
  \item{CB_alpha}{
Percentage of confidence band. Default is 0.95.
}
  \item{logit_CB}{
Whether to use logit-transformed (then transform back) confidence band. Default is FALSE.
  }
  \item{nbootstrap}{
Number of boostrap iterations. Default is 100.
}
  \item{nbin}{
Number of bins used for constructing confidence band. Default is 100.
}
  \item{verbose}{
Whether to print out messages during bootstrap. Default value is FALSE.
}
}
\value{
If global ROC is controlled by sensitivity, a list will be output including the following
\item{Sensitivity}{Vector of sensitivities;}
\item{Specificity_upper}{Upper confidence band for specificity estimations;}
\item{Specificity_lower}{Lower confidence band for specificity estimations;}
\item{global_ROC_controlled_by}{"sensitivity".}

If global ROC is controlled by Specificity, a list will be output including the following
\item{Specificity}{Vector of specificity;}
\item{Sensitivity_upper}{Upper confidence band for sensitivity estimations;}
\item{Sensitivity_lower}{Lower confidence band for sensitivity estimations;}
\item{global_ROC_controlled_by}{"specificity".}
}
\author{
Ziyi.li <zli16@mdanderson.org>
}
\examples{
n1 = n0 = 500

## generate data
Z_D1 <- rbinom(n1, size = 1, prob = 0.3)
Z_D2 <- rnorm(n1, 0.8, 1)
Z_C1 <- rbinom(n0, size = 1, prob = 0.7)
Z_C2 <- rnorm(n0, 0.8, 1)
Y_C_Z0 <- rnorm(n0, 0.1, 1)
Y_D_Z0 <- rnorm(n1, 1.1, 1)
Y_C_Z1 <- rnorm(n0, 0.2, 1)
Y_D_Z1 <- rnorm(n1, 0.9, 1)

M0 <- Y_C_Z0 * (Z_C1 == 0) + Y_C_Z1 * (Z_C1 == 1) + Z_C2
M1 <- Y_D_Z0 * (Z_D1 == 0) + Y_D_Z1 * (Z_D1 == 1) + 1.5 * Z_D2

diseaseData <- data.frame(M = M1, Z1 = Z_D1, Z2 = Z_D2)
controlData <- data.frame(M = M0, Z1 = Z_C1, Z2 = Z_C2)

userFormula = "M~Z1+Z2"
target_covariates = c(1, 0.7, 0.9)

# default nbootstrap is 100
# set nboostrap as 10 here to improve example speed
\donttest{
myROCband <- sscaROC_CB(diseaseData,
                       controlData,
                       userFormula,
                       mono_resp_method = "none",
                       target_covariates,
                       global_ROC_controlled_by = "sensitivity",
                       CB_alpha = 0.95,
                       logit_CB = FALSE,
                       nbootstrap = 10,
                       nbin = 100,
                       verbose = FALSE)
}
}
\keyword{ confidence band }% use one of  RShowDoc("KEYWORDS")
\keyword{ ROC }

