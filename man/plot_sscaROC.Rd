\name{plot_sscaROC}
\alias{plot_sscaROC}
\title{
Plot covariate-adjusted ROC for specific subpopulations.
}
\description{
Function to plot the ROC curve generated from sscaROC().
}
\usage{
plot_sscaROC(myROC, ...)
}
\arguments{
  \item{myROC}{
ROC output from sscaROC() function.
}
  \item{\dots}{
Arguments to tune generated plots.
}
}
\details{
This function can be used to plot other ROC curve, as long as the input contains two components "sensitivity" and "specificity".
}
\value{
Plot the ROC curve.
}
\author{
Ziyi Li <zli16@mdanderson.org>
}
\examples{
n1 = n0 = 1000

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

\donttest{
myROC <- sscaROC(diseaseData,
                 controlData,
                 userFormula,
                 target_covariates,
                 global_ROC_controlled_by = "sensitivity",
                 mono_resp_method = "none")
plot_sscaROC(myROC, lwd = 1.6)
}
}
\keyword{ plot }% use one of  RShowDoc("KEYWORDS")
\keyword{ ROC }% __ONLY ONE__ keyword per line

