#' @title AcompetitiveKineticEnzyme: Acompetitive Inhibition Fitting
#'
#' @description This function fits the acompetitive inhibition model to the given data. Acompetitive inhibition occurs when an inhibitor binds to both the enzyme and the enzyme-substrate complex, reducing the reaction rate. This function estimates the parameters of the acompetitive inhibition model using non-linear least squares fitting.
#'
#' @param sb Substrate concentrations. A numeric vector representing the substrate concentrations used in the enzyme assay.
#' @param rate Reaction rates. A numeric vector representing the corresponding reaction rates observed at each substrate concentration.
#' @param inh Inhibitor concentrations. A numeric vector representing the inhibitor concentrations used in the enzyme assay.
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
#' @examples
#' f<-"https://www.ugr.es/~bioest/data/inhibicionnc.txt"
#' datos<-read.table(f,sep=",",header = T)
#' AcompetitiveKineticEnzyme(sb=datos$substrate,rate=datos$rate,inh=datos$inhibitor)
#'
#' @encoding UTF-8
#' @export AcompetitiveKineticEnzyme

AcompetitiveKineticEnzyme <- function(sb, rate, inh){
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
  km <- unname(summary(sb)[2])
  vm <- max(rate)
  kic <- median(inh)
  kiu <- kic/3

  # Lineweaver Burk data
  dataset <- data.frame(sb, rate, inh)
  dataset$inv.sb <- 1/dataset$sb
  dataset$inv.rate <- 1/dataset$rate

  # How many unique inhibitors are there?
  inh.sb <- unique(dataset$inh)

  # Competitive fit
  acompetitive_nls <- nls(rate ~ ((vm * sb) / ((Km + sb) * (1 + inh/Kiu))),
                          data = dataset, start = list(Km = km, vm = vm, Kiu = kiu))

  #Coefficientes
  acompetitive_km <- unname(coef(acompetitive_nls)["Km"])
  acompetitive_vm <- unname(coef(acompetitive_nls)["vm"])
  acompetitive_kiu <- unname(coef(acompetitive_nls)["Kiu"])

  #Fitted Values
  fitted_values <- expand.grid(x = sb, inhib = inh.sb)
  fitted_values$inv.x <- 1/fitted_values$x
  fitted_values$predict <- (acompetitive_vm * fitted_values$x) / (acompetitive_km + fitted_values$x * (1 + fitted_values$inhib/acompetitive_kiu))
  fitted_values$inv.predict <- 1/fitted_values$predict

  ######################################################################
  ################## PLOT ##############################################
  ######################################################################
  plot(dataset$inv.sb, dataset$inv.rate, pch = "", main = "acompetitive inhibition")
  for (i in 1:length(inh.sb)) {
    points(subset(dataset$inv.sb, dataset$inh == inh.sb[i]),
           subset(dataset$inv.rate, dataset$inh == inh.sb[i]),
           col = "blue", pch = 1)
  }
  for (i in 1:length(inh.sb)) {
    lines(subset(fitted_values$inv.x, fitted_values$inhib == inh.sb[i]),
          subset(fitted_values$inv.predict, fitted_values$inhib == inh.sb[i]),
          col = "red")
  }

  Akaike<-AIC(acompetitive_nls)
  Bayesian<-BIC(acompetitive_nls)
  logLike <- logLik(acompetitive_nls)
  summary <- (summary(acompetitive_nls))
  km <- summary$coefficients[1,1]
  vm <- summary$coefficients[1,2]
  kic <- summary$coefficients[1,3]
  kiu <- summary$coefficients[1,4]
  stderrorkm <- summary$coefficients[2,1]
  stderrorvm <- summary$coefficients[2,2]
  stderrorkic <- summary$coefficients[2,3]
  stderrorkiu <- summary$coefficients[2,4]
  resumen <- summary(acompetitive_nls)
  coeficientedeterminacion<-resumen$sigma

  # Return the estimated parameters
  Resultados <- list(AIC = Akaike, BIC = Bayesian,logLike=logLike,Km = km, Vm = vm, Kic = kic, Kiu = kiu, StandardErrorvm=stderrorvm, StandardErrorkm=stderrorkm,
                     StandardErrorkic=stderrorkic, StandardErrorkiu=stderrorkiu,ResidualStandardError=coeficientedeterminacion)
  return(Resultados)
}
