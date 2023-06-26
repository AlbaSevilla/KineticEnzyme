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
#' \item{kic}{The estimated inhibition constant.}
#' \item{StandardErrorvm}{The standard error associated with the estimated maximum reaction rate (Vm).}
#' \item{StandardErrorkm}{The standard error associated with the estimated Michaelis-Menten constant (Km).}
#' \item{StandardErrorkic}{The standard error associated with the estimated inhibition constant (kic).}
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
    install.packages("lmtest", repos = "https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz", type = "source")
  if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")

  library(car)
  library(minpack.lm)
  library(lmtest)
  library(dplyr)

  # DATA
  dataset <- data.frame(sb, inh,rate)

  # initial values for the parameters (considering the L-B's transformation)

  condition<-(dataset$sb>0 & dataset$rate>0)
  dataset2<-subset.data.frame(dataset,condition)
  dataset2$invsb <- 1/dataset2$sb
  dataset2$invrate <- 1/dataset2$rate
  lb<-lm(dataset2$invrate~dataset2$invsb)
  vm<-1/coefficients(lb)[1]
  km<-coefficients(lb)[2]*vm
  kic<-1

  # Lineweaver Burk data
  dataset$inv.sb <- 1/dataset$sb
  dataset$inv.rate <- 1/dataset$rate

  #Dataset without infinite data
  dataset <- dataset %>%
    filter_all(all_vars(!is.infinite(.)))

  # How many unique inhibitors are there?
  inh.sb <- unique(dataset$inh)

  # Competitive fit
  acompetitive_nls <- nls(rate ~ ((vm * sb) / (Km + sb * (1 + inh/kic))),
                          data = dataset, start = list(Km = km, vm = vm, kic = kic))

  #Coefficientes
  acompetitive_km <- unname(coef(acompetitive_nls)["Km"])
  acompetitive_vm <- unname(coef(acompetitive_nls)["vm"])
  acompetitive_kic <- unname(coef(acompetitive_nls)["kic"])

  #Fitted Values
  fitted_values <- expand.grid(x = sb, inhib = inh.sb)
  fitted_values$inv.x <- 1/fitted_values$x
  fitted_values$predict <- (acompetitive_vm * fitted_values$x) / (acompetitive_km + fitted_values$x * (1 + fitted_values$inhib/acompetitive_kic))
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
  logLike<-logLik(acompetitive_nls)
  summary <- (summary(acompetitive_nls))
  km <- summary$coefficients[1,1]
  vm <- summary$coefficients[2,1]
  kic <- summary$coefficients[3,1]
  stderrorkm <- summary$coefficients[1,2]
  stderrorvm <- summary$coefficients[2,2]
  stderrorkic <- summary$coefficients[3,2]
  resumen <- summary(acompetitive_nls)
  coeficientedeterminacion<-resumen$sigma

  # Return the estimated parameters
  Resultados <- list(AIC = Akaike, BIC = Bayesian,logLike=logLike,Km = km, Vm = vm, kic = kic,StandardErrorvm=stderrorvm, StandardErrorkm=stderrorkm,
                     StandardErrorkic=stderrorkic,ResidualStandardError=coeficientedeterminacion)
  return(Resultados)
}
