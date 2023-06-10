#' @title NonCompetitiveKineticEnzyme: noncompetitive inhibition Fitting
#' @description
#' This function fits the noncompetitive inhibition model to the given data.
#' The noncompetitive inhibition model describes the relationship between substrate concentration,
#' reaction rate, and inhibitor concentration in enzyme kinetics.
#'
#' The function takes substrate, rate, and inhibitor data as input and performs the model fitting.
#' It estimates the parameters Km, Vm, Kic, and Kiu using the nonlinear least squares method.
#' The function also provides a plot of the fitted lines and returns a data frame containing the
#' estimated parameters and other analysis results.
#'
#' @param sb
#' Substrate data. A numeric vector containing the substrate concentrations.
#'
#' @param rate
#' Rate data. A numeric vector containing the corresponding reaction rates.
#'
#' @param inh
#' Inhibitor data. A numeric vector containing the inhibitor concentrations.
#'
#' @return A data frame containing the estimated parameters and other analysis results. The data frame includes the following columns:
#' \item{AIC}{The Akaike Information Criterion (AIC) value for model selection.}
#' \item{BIC}{The Bayesian Information Criterion (BIC) value for model selection.}
#' \item{Km}{The estimated Michaelis-Menten constant.}
#' \item{Vm}{The estimated maximum reaction rate.}
#' \item{Kic}{The estimated inhibition constant.}
#' \item{Kiu}{The estimated uncompetitive inhibition constant.}
#' \item{StandardErrorvm}{The standard error associated with the estimated maximum reaction rate (Vm).}
#' \item{StandardErrorkm}{The standard error associated with the estimated Michaelis-Menten constant (Km).}
#' \item{StandardErrorkic}{The standard error associated with the estimated inhibition constant (Kic).}
#' \item{StandardErrorkiu}{The standard error associated with the estimated uncompetitive inhibition constant (Kiu).}
#'
#' @export NonCompetitiveKineticEnzyme
#' @encoding UTF-8

NonCompetitiveKineticEnzyme <- function(sb,rate,inh){
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest",repos="https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",type="source")
  library(car)
  library(minpack.lm)
  library(lmtest)

  # DATA
  km<-unname(summary(sb)[2])
  vm<-max(rate)
  kic<-median(inh)
  kiu<-kic/3

  #Lineweaver Burk data
  dataset <- data.frame(sb,rate,inh)
  dataset$inv.sb <- 1/dataset$sb
  dataset$inv.rate <- 1/dataset$rate

  #How many unique of inhibitor are?
  inh.sb <- unique(dataset$inh)

  #non competitive fit
  noncompetitive_nls <- nls(rate ~ ((vm * sb) / (Km*(1+inh/Kic) +
                                                   sb*(1+inh/Kiu))), data=dataset, start=list(Km=km, vm=vm, Kic=kic, Kiu=kiu))

  #Coefficients
  noncompetitive_km <- unname(coef(noncompetitive_nls)["Km"])
  noncompetitive_vm <- unname(coef(noncompetitive_nls)["vm"])
  noncompetitive_kic <- unname(coef(noncompetitive_nls)["Kic"])
  noncompetitive_kiu <- unname(coef(noncompetitive_nls)["Kiu"])

  #Fitted Values
  fitted_values <- expand.grid(x=sb, inhib=inh.sb)
  fitted_values$inv.x <- 1/fitted_values$x
  fitted_values$noncompetitive.y <- (noncompetitive_vm*fitted_values$x)/(noncompetitive_km*(1+fitted_values$inhib/noncompetitive_kic)+fitted_values$x*
                                                                           (1+fitted_values$inhib/noncompetitive_kiu))
  fitted_values$inv.noncompetitive.y <- 1/fitted_values$noncompetitive.y

  # plot lines of best fit - noncompetitive
  # generate a blank plot and then plot the raw data
  plot(dataset$inv.sb,dataset$inv.rate, pch="", main="noncompetitive inhibition")
  for (i in 1:length(inh.sb)) {
    points(subset(dataset$inv.sb, dataset$inh==inh.sb[i]),
           subset(dataset$inv.rate, dataset$inh==inh.sb[i]),
           col="blue", pch=1)
  }
  for (i in 1:length(inh.sb)) {
    lines(subset(fitted_values$inv.x, fitted_values$inhib==inh.sb[i]),
          subset(fitted_values$inv.noncompetitive.y, fitted_values$inhib==inh.sb[i]),
          col="red")
  }
  Akaike<-AIC(noncompetitive_nls)
  Bayesian<-BIC(noncompetitive_nls)
  logLike<-logLik(noncompetitive_nls)
  summary <- (summary(noncompetitive_nls))
  km <- summary$coefficients[1,1]
  vm <- summary$coefficients[1,2]
  kic <- summary$coefficients[1,3]
  kiu <- summary$coefficients[1,4]
  stderrorkm <- summary$coefficients[2,1]
  stderrorvm <- summary$coefficients[2,2]
  stderrorkic <- summary$coefficients[2,3]
  stderrorkiu <- summary$coefficients[2,4]
  resumen <- summary(noncompetitive_nls)
  coeficientedeterminacion<-resumen$sigma

  # Return the estimated parameters
  Resultados <- list(AIC = Akaike, BIC = Bayesian,logLike=logLike,Km = km, Vm = vm, Kic = kic, Kiu = kiu, StandardErrorvm=stderrorvm, StandardErrorkm=stderrorkm,
                     StandardErrorkic=stderrorkic, StandardErrorkiu=stderrorkiu,ResidualStandardError=coeficientedeterminacion)
  return(Resultados)
}
