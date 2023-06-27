#' @title NonCompetitiveKineticEnzyme: noncompetitive inhibition Fitting
#' @description
#' This function fits the noncompetitive inhibition model to the given data.
#' The noncompetitive inhibition model describes the relationship between substrate concentration,
#' reaction rate, and inhibitor concentration in enzyme kinetics.
#'
#' The function takes substrate, rate, and inhibitor data as input and performs the model fitting.
#' It estimates the parameters Km, Vm, kic, using the nonlinear least squares method.
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
#' \item{kic}{The estimated inhibition constant.}
#' \item{StandardErrorvm}{The standard error associated with the estimated maximum reaction rate (Vm).}
#' \item{StandardErrorkm}{The standard error associated with the estimated Michaelis-Menten constant (Km).}
#' \item{StandardErrorkic}{The standard error associated with the estimated inhibition constant (kic).}
#'
#' @examples
#' f<-"https://www.ugr.es/~bioest/data/inhibicionnc.txt"
#' datos<-read.table(f,sep=",",header = TRUE)
#' head(datos)
#' NonCompetitiveKineticEnzyme(sb=datos$substrate,rate=datos$rate,inh=datos$inhibitor)
#'
#' @export NonCompetitiveKineticEnzyme
#' @encoding UTF-8

NonCompetitiveKineticEnzyme <- function(sb,rate,inh){
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest", repos = "https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz", type = "source")
  if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", repos = "https://cran.r-project.org/src/contrib/dplyr_1.1.2.tar.gz", type = "source")

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

  #How many unique of inhibitor are?
  inh.sb <- unique(dataset$inh)

  #non competitive fit
  noncompetitive_nls <- nls(rate ~ ((vm * sb) / (Km*(1+inh/kic) +
                                                   sb*(1+inh/kic))), data=dataset, start=list(Km=km, vm=vm, kic=kic))

  #Coefficients
  noncompetitive_km <- unname(coef(noncompetitive_nls)["Km"])
  noncompetitive_vm <- unname(coef(noncompetitive_nls)["vm"])
  noncompetitive_kic <- unname(coef(noncompetitive_nls)["kic"])

  #Fitted Values
  fitted_values <- expand.grid(x=sb, inhib=inh.sb)
  fitted_values$inv.x <- 1/fitted_values$x
  fitted_values$noncompetitive.y <- (noncompetitive_vm*fitted_values$x)/(noncompetitive_km*(1+fitted_values$inhib/noncompetitive_kic)+fitted_values$x*
                                                                           (1+fitted_values$inhib/noncompetitive_kic))
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
  vm <- summary$coefficients[2,1]
  kic <- summary$coefficients[3,1]
  stderrorkm <- summary$coefficients[1,2]
  stderrorvm <- summary$coefficients[2,2]
  stderrorkic <- summary$coefficients[3,2]
  resumen <- summary(noncompetitive_nls)
  coeficientedeterminacion<-resumen$sigma

  # Return the estimated parameters
  Resultados <- list(AIC = Akaike, BIC = Bayesian,logLike=logLike,Km = km, Vm = vm, kic = kic,StandardErrorvm=stderrorvm, StandardErrorkm=stderrorkm,
                     StandardErrorkic=stderrorkic,ResidualStandardError=coeficientedeterminacion)
  return(Resultados)
}
