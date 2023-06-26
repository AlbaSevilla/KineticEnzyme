#' @title CompetitiveKineticEnzyme: Competitive Inhibition Fitting
#'
#' @description This function fits the Competitive Inhibition model to the given data using the least squares method. The Competitive Inhibition model describes the enzyme kinetics in the presence of a competitive inhibitor.
#'
#' @param sb A numeric vector representing the substrate concentrations.
#' @param rate A numeric vector representing the corresponding reaction rates.
#' @param inh A numeric vector representing the concentrations of the competitive inhibitor.
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
#' @examples
#' f<-"https://www.ugr.es/~bioest/data/inhibicionnc.txt"
#' datos<-read.table(f,sep=",",header = T)
#' head(datos)
#' CompetitiveKineticEnzyme(sb=datos$substrate,rate=datos$rate,inh=datos$inhibitor)
#'
#' @export CompetitiveKineticEnzyme
#' @encoding UTF-8

CompetitiveKineticEnzyme <- function(sb,rate,inh){
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

  #How many unique of inhibitor are?
  inh.sb <- unique(dataset$inh)

  # Competitive fit
  competitive_nls <- nls(rate ~ ((vm * sb) / (Km*(1+inh/kic) + sb)),
                         data=dataset, start=list(Km=km, vm=vm, kic=kic))

  #Coefficients
  competitive_km <- unname(coef(competitive_nls)[1])
  competitive_vm <- unname(coef(competitive_nls)[2])
  competitive_kic <- unname(coef(competitive_nls)[3])

  #Fitted Values
  fittedvalues <- expand.grid(x=sb, inhib=inh.sb)
  fittedvalues$inv.x <- 1/fittedvalues$x
  fittedvalues$predict <- (competitive_vm*fittedvalues$x)/(competitive_km*(1+fittedvalues$inhib/competitive_kic)+fittedvalues$x)
  fittedvalues$inv.predict <- 1/fittedvalues$predict

  ######################################################################
  ################## PLOT ##############################################
  ######################################################################
  plot(dataset$inv.sb,dataset$inv.rate, pch="", main="Competitive inhibition",
       xlab="1/sb",ylab="1/vm")
  for (i in 1:length(inh.sb)) {
    points(subset(dataset$inv.sb, dataset$inh==inh.sb[i]),
           subset(dataset$inv.rate, dataset$inh==inh.sb[i]),
           col="blue", pch=1)
  }
  for (i in 1:length(inh.sb)) {
    lines(subset(fittedvalues$inv.x, fittedvalues$inhib==inh.sb[i]),
          subset(fittedvalues$inv.predict, fittedvalues$inhib==inh.sb[i]),col="red")
  }

  Akaike<-AIC(competitive_nls)
  Bayesian<-BIC(competitive_nls)
  logLike<-logLik(competitive_nls)
  summary <- (summary(competitive_nls))
  km <- summary$coefficients[1,1]
  vm <- summary$coefficients[2,1]
  kic <- summary$coefficients[3,1]
  stderrorkm <- summary$coefficients[1,2]
  stderrorvm <- summary$coefficients[2,2]
  stderrorkic <- summary$coefficients[3,2]
  resumen <- summary(competitive_nls)
  coeficientedeterminacion<-resumen$sigma

  # Return the estimated parameters
  Resultados <- list(AIC = Akaike, BIC = Bayesian,logLike=logLike,Km = km, Vm = vm, kic = kic,StandardErrorvm=stderrorvm, StandardErrorkm=stderrorkm,
                     StandardErrorkic=stderrorkic,ResidualStandardError=coeficientedeterminacion)
  return(Resultados)
}
