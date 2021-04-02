\name{caROC_CB}
\alias{caROC_CB}
\title{
Get confidence band for covariate-adjusted ROC curve.
}
\description{
Use this function to compute the confidence band for covariate-adjusted ROC curve, with or without monotonicity respecting methods.
}
\usage{
caROC_CB(diseaseData, controlData, userFormula,
mono_resp_method, global_ROC_controlled_by = "sensitivity",
CB_alpha = 0.95, logit_CB = FALSE, nbootstrap = 100,
nbin = 100, verbose = FALSE)
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

### calculate confidence band by controlling sensitivity
### using different monotonicity respecting methods

\donttest{
ROC_CB1 <- caROC_CB(diseaseData,controlData,userFormula,
                      mono_resp_method = "none",
                      CB_alpha = 0.95,
                      nbin = 100,verbose = FALSE)
ROC_CB2 <- caROC_CB(diseaseData,controlData,userFormula,
                      mono_resp_method = "ROC",
                      CB_alpha = 0.95,
                      nbin = 100,verbose = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ confidence band }% use one of  RShowDoc("KEYWORDS")
\keyword{ ROC }
