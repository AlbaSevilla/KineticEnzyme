#' @title InhibitionComparisonKineticEnzyme: Compare Different Enzyme Kinetic Inhibition Models
#' @description
#' This function compares different enzyme kinetic inhibition models by fitting them to the given data.
#' The function requires the 'KineticEnzyme' package to be installed.
#' The function calculates the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)
#' for each model and provides additional analysis results such as the estimated parameters and their
#' corresponding standard errors.
#' The function also plots the data and the fitted curves for each model.
#'
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @import car
#' @import lmtest
#' @import carData
#'
#' @param sb A vector or column of substrate concentrations.
#' @param rate A vector or column of corresponding reaction rates.
#' @param inh A vector or column of inhibitor concentrations.
#'
#' @return A data frame containing the results of the three inhibition models.
#' The data frame includes the following columns:
#' \item{AIC}{Akaike Information Criterion for each model.}
#' \item{BIC}{Bayesian Information Criterion for each model.}
#' \item{Standard error for Vm}{Standard error of the maximum velocity parameter for each model.}
#' \item{Standard error for Km}{Standard error of the Michaelis constant parameter for each model.}
#' \item{Standard error for Kic}{Standard error of the inhibitor constant parameter for each model.}
#' \item{Standard error for Kiu}{Standard error of the uncompetitive inhibitor constant parameter for each model.}
#'
#' @export InhibitionComparisonKineticEnzyme
#' @encoding UTF-8

InhibitionComparisonKineticEnzyme <- function(sb, rate, inh) {
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest",repos="https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",type="source")
  library(car)
  library(minpack.lm)
  library(lmtest)

  par(mfrow = c(2, 2))

  noncompetitive <- KineticEnzyme::NonCompetitiveKineticEnzyme(sb,rate,inh)
  noncompetitiveAIC <- noncompetitive$AIC
  noncompetitiveBIC <- noncompetitive$BIC
  noncompetitivelogLike<-noncompetitive$logLike
  noncompetitiveKm <- noncompetitive$Km
  noncompetitiveVm <- noncompetitive$Vm
  noncompetitiveKic <- noncompetitive$Kic
  noncompetitiveKiu <- noncompetitive$Kiu
  noncompetitiveStandardErrorKm <- noncompetitive$StandardErrorkm
  noncompetitiveStandardErrorVm <- noncompetitive$StandardErrorvm
  noncompetitiveStandardErrorKic <- noncompetitive$StandardErrorkic
  noncompetitiveStandardErrorKiu <- noncompetitive$StandardErrorkiu

  competitive <- KineticEnzyme::CompetitiveKineticEnzyme(sb, rate, inh)
  competitiveAIC <- competitive$AIC
  competitiveBIC <- competitive$BIC
  competitivelogLike<-competitive$logLike
  competitiveKm <- competitive$Km
  competitiveVm <- competitive$Vm
  competitiveKic <- competitive$Kic
  competitiveKiu <- competitive$Kiu
  competitiveStandardErrorKm <- competitive$StandardErrorkm
  competitiveStandardErrorVm <- competitive$StandardErrorvm
  competitiveStandardErrorKic <- competitive$StandardErrorkic
  competitiveStandardErrorKiu <- competitive$StandardErrorkiu

  acompetitive <- KineticEnzyme::AcompetitiveKineticEnzyme(sb, rate, inh)
  acompetitiveAIC <- acompetitive$AIC
  acompetitiveBIC <- acompetitive$BIC
  acompetitivelogLike<-acompetitive$logLike
  acompetitiveKm <- acompetitive$Km
  acompetitiveVm <- acompetitive$Vm
  acompetitiveKic <- acompetitive$Kic
  acompetitiveKiu <- acompetitive$Kiu
  acompetitiveStandardErrorKm <- acompetitive$StandardErrorkm
  acompetitiveStandardErrorVm <- acompetitive$StandardErrorvm
  acompetitiveStandardErrorKic <- acompetitive$StandardErrorkic
  acompetitiveStandardErrorKiu <- acompetitive$StandardErrorkiu

  vAIC <- c(noncompetitiveAIC, competitiveAIC, acompetitiveAIC)
  vAICcolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorAIC <- min(vAIC)
  menorAICcolname <- vAICcolnames[which.min(vAIC)]

  vBIC <- c(noncompetitiveBIC, competitiveBIC, acompetitiveBIC)
  vBICcolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorBIC <- min(vBIC)
  menorBICcolname <- vBICcolnames[which.min(vBIC)]

  vlogLike <- c(noncompetitivelogLike, competitivelogLike, acompetitivelogLike)
  vlogLikecolnames <- c("Non competitive", "Competitive", "Acompetitive")
  mayorlogLike <- max(vlogLike)
  mayorlogLikecolname <- vlogLikecolnames[which.max(vlogLike)]

  StandardErrorVm <- c(noncompetitiveVm, competitiveVm, acompetitiveVm)
  StandardErrorVmcolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorStandardErrorVm <- min(StandardErrorVm)
  menorStandardErrorVmcolname <- StandardErrorVmcolnames[which.min(StandardErrorVm)]

  StandardErrorKm <- c(noncompetitiveKm, competitiveKm, acompetitiveKm)
  StandardErrorKmcolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorStandardErrorKm <- min(StandardErrorKm)
  menorStandardErrorKmcolname <- StandardErrorKmcolnames[which.min(StandardErrorKm)]

  StandardErrorKic <- c(noncompetitiveKic, competitiveKic, acompetitiveKic)
  StandardErrorKiccolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorStandardErrorKic <- min(StandardErrorKic)
  menorStandardErrorKiccolname <- StandardErrorKiccolnames[which.min(StandardErrorKic)]

  StandardErrorKiu <- c(noncompetitiveKiu, competitiveKiu, acompetitiveKiu)
  StandardErrorKiucolnames <- c("Non competitive", "Competitive", "Acompetitive")
  menorStandardErrorKiu <- min(StandardErrorKiu)
  menorStandardErrorKiucolname <- StandardErrorKiucolnames[which.min(StandardErrorKiu)]

  AIC <- c(noncompetitiveAIC, competitiveAIC, acompetitiveAIC)
  BIC <- c(noncompetitiveBIC, competitiveBIC, acompetitiveBIC)
  logLike <- c(noncompetitivelogLike, competitivelogLike, acompetitivelogLike)
  StandardErrorVm <- c(noncompetitiveStandardErrorVm, competitiveStandardErrorVm, acompetitiveStandardErrorVm)
  StandardErrorKm <- c(noncompetitiveStandardErrorKm, competitiveStandardErrorKm, acompetitiveStandardErrorKm)
  StandardErrorKic <- c(noncompetitiveStandardErrorKic, competitiveStandardErrorKic, acompetitiveStandardErrorKic)
  StandardErrorKiu <- c(noncompetitiveStandardErrorKiu, competitiveStandardErrorKiu, acompetitiveStandardErrorKiu)

  df <- matrix(c(AIC, BIC,logLike, StandardErrorVm, StandardErrorKm, StandardErrorKic, StandardErrorKiu), ncol = 7, byrow = FALSE)
  rownames(df) <- c("Non competitive", "Competitive", "Acompetitive")
  colnames(df) <- c("AIC", "BIC","logLike", "Standard error for Vm", "Standard error for Km",
                    "Standard error for Kic", "Standard error for Kiu")

  # Which model is better?
  cat("The smallest value of the AIC is:", menorAIC, "corresponding to the model of", menorAICcolname, "\n")
  cat("The smallest value of the BIC is:", menorBIC, "corresponding to the model of", menorBICcolname, "\n")
  cat("The biggest value of the Log Likelihood is:", mayorlogLike, "corresponding to the model of", mayorlogLikecolname, "\n")
  cat("The model with the smallest standard error of Vm is:", menorStandardErrorVm, "corresponding to the model of", menorStandardErrorVmcolname, "\n")
  cat("The model with the smallest standard error of Km is:", menorStandardErrorKm, "corresponding to the model of", menorStandardErrorKmcolname, "\n")
  cat("The model with the least standard error of Kic is:", menorStandardErrorKic, "corresponding to the model of", menorStandardErrorKiccolname, "\n")
  cat("The model with the smallest standard error of Kiu is:", menorStandardErrorKiu, "corresponding to the model of", menorStandardErrorKiucolname, "\n")
  return(df)
}
