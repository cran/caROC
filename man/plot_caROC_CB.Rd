\name{plot_caROC_CB}
\alias{plot_caROC_CB}
\title{
Plot confidence band of covariate-adjusted ROC.
}
\description{
A function to plot the confidence band of covariate-adjusted ROC.
}
\usage{
plot_caROC_CB(myROC_CB, add = TRUE, ...)
}
\arguments{
  \item{myROC_CB}{
Output from caROC_CB() function.
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
Ziyi Li<ziyi.li@emory.edu>
}
\examples{
library(caROC)
n1 = n0 = 100

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
formula = "M~Z"

ROC_CB1 <- caROC_CB(diseaseData,controlData,formula,
                       mono_resp_method = "none",
                       CB_alpha = 0.95,
                       nbin = 100,verbose = FALSE)
### plot confidence band individually
plot_caROC_CB(ROC_CB1, add = FALSE, lty = 2, col = "blue")

### plot confidence band together with the ROC curve
ROC1 <- caROC(diseaseData,controlData,formula,
                 mono_resp_method = "none", verbose = FALSE)
plot_caROC(ROC1)
plot_caROC_CB(ROC_CB1, add = TRUE, lty = 2, col = "blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Plot }% use one of  RShowDoc("KEYWORDS")
\keyword{ ROC }% use one of  RShowDoc("KEYWORDS")
\keyword{ Confidence band }% __ONLY ONE__ keyword per line
