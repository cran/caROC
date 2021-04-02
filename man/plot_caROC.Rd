\name{plot_caROC}
\alias{plot_caROC}
\title{
Plot covariate-adjusted ROC.
}
\description{
Function to plot the ROC curve generated from caROC().
}
\usage{
plot_caROC(myROC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myROC}{
ROC output from caROC() function.
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

ROC1 <- caROC(diseaseData,controlData,userFormula,
                 mono_resp_method = "none")
ROC2 <- caROC(diseaseData,controlData,userFormula,
                 mono_resp_method = "ROC")

plot_caROC(ROC1)
plot_caROC(ROC2, col = "blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }% use one of  RShowDoc("KEYWORDS")
\keyword{ ROC }% __ONLY ONE__ keyword per line
