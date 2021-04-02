\name{plot_sscaROC_CB}
\alias{plot_sscaROC_CB}
\title{
Plot confidence band of covariate-adjusted ROC in specific subpopulations.
}
\description{
A function to plot the confidence band of covariate-adjusted ROC in specific subpopulations.
}
\usage{
plot_sscaROC_CB(myROC_CB, add = TRUE, ...)
}
\arguments{
  \item{myROC_CB}{
Output from sscaROC_CB() function.
}
  \item{add}{
Whether to add confidence band to existing plot (TRUE) or draw a new one (FALSE). Default is TRUE.
}
  \item{\dots}{
Any parameters related with the plot.
}
}

\value{
No values will be return. This function is for plotting only.
}
\author{
Ziyi Li<zli16@mdanderson.org>
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

\donttest{
# example that takes more than a minute to run
myROC <- sscaROC(diseaseData,
                controlData,
                userFormula,
                target_covariates,
                global_ROC_controlled_by = "sensitivity",
                mono_resp_method = "none")

# default nbootstrap is 100
# set nboostrap as 10 here to improve example speed
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

plot_sscaROC(myROC, lwd = 1.6)
plot_sscaROC_CB(myROCband, col = "purple", lty = 2)
}
}
\keyword{ confidence band }
\keyword{ ROC }

