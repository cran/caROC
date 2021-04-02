## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  library(devtools)
#  install_github("ziyili20/caROC")

## ----quick_start, eval = FALSE, message=FALSE---------------------------------
#  library(caROC)
#  ## get specificity at controlled sensitivity levels 0.2, 0.8, 0.9
#  caROC(diseaseData,controlData,formula,
#        control_sensitivity = c(0.2,0.8, 0.9),
#        control_specificity = NULL)
#  
#  ## get covariate-adjusted ROC curve with curve-based monotonizing method
#  curveROC <- caROC(diseaseData,controlData,formula,
#              mono_resp_method = "curve",
#              verbose = FALSE)

## ----getdata, eval = TRUE, warning=FALSE, message=FALSE-----------------------
library(caROC)
### n1: number of cases
### n0: number of controls
n1 = n0 = 1000

## Z_D and Z_C are the covariates in the disease and control groups
Z_D1 <- rbinom(n1, size = 1, prob = 0.3)
Z_D2 <- rnorm(n1, 0.8, 1)

Z_C1 <- rbinom(n0, size = 1, prob = 0.7)
Z_C2 <- rnorm(n0, 0.8, 1)

Y_C_Z0 <- rnorm(n0, 0.1, 1)
Y_D_Z0 <- rnorm(n1, 1.1, 1)
Y_C_Z1 <- rnorm(n0, 0.2, 1)
Y_D_Z1 <- rnorm(n1, 0.9, 1)

## M0 and M1 are the outcome of interest (biomarker to be evaluated) in the control and disease groups
M0 <- Y_C_Z0 * (Z_C1 == 0) + Y_C_Z1 * (Z_C1 == 1) + Z_C2
M1 <- Y_D_Z0 * (Z_D1 == 0) + Y_D_Z1 * (Z_D1 == 1) + 1.5 * Z_D2

diseaseData <- data.frame(M = M1, Z1 = Z_D1, Z2 = Z_D2)
controlData <- data.frame(M = M0, Z1 = Z_C1, Z2 = Z_C2)

## we are interested in evaluating biomarker M while adjusting for covariate Z
userFormula = "M~Z1+Z2"

## ----controlspec, warning=FALSE, message=FALSE--------------------------------
caROC(diseaseData,controlData,userFormula,
      control_sensitivity = c(0.2,0.8, 0.9),
      control_specificity = NULL,
      mono_resp_method = "ROC",
      whichSE = "bootstrap",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)

## ----controlsens, warning=FALSE, message=FALSE--------------------------------
caROC(diseaseData,controlData,userFormula,
      control_sensitivity = NULL,
      control_specificity = c(0.7,0.8, 0.9),
      mono_resp_method = "none",
      whichSE = "sample",nbootstrap = 100,
      CI_alpha = 0.95, logit_CI = TRUE)

## ----controlspecss, warning=FALSE, message=FALSE------------------------------
target_covariates = c(1, 0.7, 0.9)

sscaROC(diseaseData,controlData,
               userFormula = userFormula,
               control_sensitivity = c(0.2,0.8, 0.9),
               target_covariates = target_covariates,
               control_specificity = NULL,
               mono_resp_method = "none",
               whichSE = "sample",nbootstrap = 100,
               CI_alpha = 0.95, logit_CI = TRUE)

## ----multicontrolspecss, message=FALSE, warning=FALSE-------------------------
target_covariates = matrix(c(1, 0.7, 0.9,
                      1, 0.8, 0.8), 2, 3, byrow = TRUE)
sscaROC(diseaseData,controlData,
               userFormula = userFormula,
               control_sensitivity = c(0.2,0.8, 0.9),
               target_covariates = target_covariates,
               control_specificity = NULL,
               mono_resp_method = "none",
               whichSE = "sample",nbootstrap = 100,
               CI_alpha = 0.95, logit_CI = TRUE)

## ----ROC, warning=FALSE, message=FALSE----------------------------------------
### ROC with curve-based monotonicity restoration
curveROC <- caROC(diseaseData,controlData,userFormula,
                 mono_resp_method = "ROC", 
                 verbose = FALSE)

## ----plotROC, warning=FALSE, message=FALSE------------------------------------
oldpar <- par()
par(mar = c(3, 3, 2, 0.3), mgp = c(1.2, 0.3, 0))
plot_caROC(curveROC)
par(oldpar)

## ----ROC2, warning=FALSE, message=FALSE---------------------------------------
curveROC_CB <- caROC_CB(diseaseData,controlData,
						userFormula, 
						mono_resp_method = "ROC",
						CB_alpha = 0.95,
						nbin = 100,verbose = FALSE)

## ----plotROCband, warning=FALSE, message=FALSE--------------------------------
oldpar <- par()
par(mar = c(3, 3, 2, 0.3), mgp = c(1.2, 0.3, 0))
plot_caROC_CB(curveROC_CB, add = FALSE, lty = 2, col = "blue")  
par(oldpar)

## ----plotROCband2, warning=FALSE, message=FALSE-------------------------------
oldpar <- par()
par(mar = c(3, 3, 2, 0.3), mgp = c(1.2, 0.3, 0))
plot_caROC(curveROC)
plot_caROC_CB(curveROC_CB, add = TRUE, lty = 2, col = "blue")
par(oldpar)

## ----ssROC, warning=FALSE, message=FALSE--------------------------------------
target_covariates = c(1, 0.7, 0.9)
myROC <- sscaROC(diseaseData,
                 controlData,
                 userFormula,
                 target_covariates,
                 global_ROC_controlled_by = "sensitivity",
                 mono_resp_method = "none")
oldpar <- par()
par(mar = c(3, 3, 2, 0.3), mgp = c(1.2, 0.3, 0))
plot_sscaROC(myROC, lwd = 1.6)
par(oldpar)

## ----ssROCband, eval=FALSE----------------------------------------------------
#  myROCband <- sscaROC_CB(diseaseData,
#                          controlData,
#                          userFormula,
#                          mono_resp_method = "none",
#                          target_covariates,
#                          global_ROC_controlled_by = "sensitivity",
#                          CB_alpha = 0.95,
#                          logit_CB = FALSE,
#                          nbootstrap = 100,
#                          nbin = 100,
#                          verbose = FALSE)
#  oldpar <- par()
#  par(mar = c(3, 3, 2, 0.3), mgp = c(1.2, 0.3, 0))
#  plot_sscaROC_CB(myROCband, col = "purple", lty = 2)
#  par(oldpar)

## ----treshold, warning=FALSE, message=FALSE-----------------------------------
### this is the given covariates of interest
new_covariates <- data.frame(M = 1,
                      Z1 = 0.7,
                      Z2 = 0.9)
### controlling sensitivity levels
caThreshold(userFormula, new_covariates,
            diseaseData = diseaseData,
            controlData = NULL,
            control_sensitivity = c(0.7,0.8,0.9),
            control_specificity = NULL)
            
### controlling specificity levels
caThreshold(userFormula,new_covariates,
            diseaseData = NULL,
            controlData = controlData,
            control_sensitivity = NULL,
            control_specificity = c(0.7,0.8,0.9))

