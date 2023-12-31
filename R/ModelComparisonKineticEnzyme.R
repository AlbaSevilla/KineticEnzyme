#' @title ModelComparisonKineticEnzyme: Compare Different Enzyme Kinetic Models
#' @description
#'     This function compares different enzyme kinetic models by fitting them to the given data.
#'     It calculates the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)
#'     for each model, as well as the estimated parameters and their standard errors.
#'     The function supports the following enzyme kinetic models: Michaelis-Menten, Lineweaver-Burk,
#'     Eadie-Scatchard, Woolf-Augustinson-Hofstee, and Hanes-Woolf.
#'
#' @param substrate A vector or column of substrate concentrations.
#' @param velocity A vector or column of corresponding velocity measurements.
#'
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @import car
#' @import lmtest
#' @import carData
#'
#' @return A data frame containing the following columns:
#' \item{AIC}{Akaike Information Criterion for each model.}
#' \item{BIC}{Bayesian Information Criterion for each model.}
#' \item{vmax}{Estimated maximum velocity parameter for each model.}
#' \item{km}{Estimated Michaelis-Menten constant parameter for each model.}
#' \item{StandardErrorvm_delta}{Standard error of the estimated maximum velocity parameter for each model using delta method.}
#' \item{StandardErrorkm_delta}{Standard error of the estimated Michaelis-Menten constant parameter for each model using delta method.}
#'
#' The row names of the data frame correspond to the names of the enzyme kinetic models.
#'
#' Additionally, the function prints the following information:
#' \item{Enzyme kinetic model with the lowest AIC value}{The enzyme kinetic model with the lowest AIC value.}
#' \item{Enzyme kinetic model with the lowest BIC value}{The enzyme kinetic model with the lowest BIC value.}
#' \item{Enzyme kinetic model with the lowest Log Likelihood value}{The enzyme kinetic model with the lowest Log Likelihood value.}
#' \item{Enzyme kinetic model with the lowest standard error of the maximum velocity (vm) parameter}{The enzyme kinetic model with the lowest standard error of the maximum velocity (vm) parameter.}
#' \item{Enzyme kinetic model with the lowest standard error of the Michaelis-Menten (km) constant parameter}{The enzyme kinetic model with the lowest standard error of the Michaelis-Menten (km) constant parameter.}
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' ModelComparisonKineticEnzyme(substrate=data$s,velocity=data$v)
#' @encoding UTF-8
#' @export ModelComparisonKineticEnzyme

ModelComparisonKineticEnzyme <- function(substrate,velocity){
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest",repos="https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",type="source")
  library(car)
  library(minpack.lm)
  library(lmtest)

  par(mfrow = c(3, 2))

  michaelismenten<-KineticEnzyme::MMKineticEnzyme(substrate,velocity)
  michaelismentenAIC<-michaelismenten$AIC
  michaelismentenBIC<-michaelismenten$BIC
  michaelismentenlogLike<-michaelismenten$logLike
  michaelismentenvmax<-michaelismenten$vmax
  michaelismentenkm<-michaelismenten$km
  michaelismentenStandardErrorvm_delta<-michaelismenten$StandardErrorvm
  michaelismentenStandardErrorkm_delta<-michaelismenten$StandardErrorkm

  lineweaverburk<-KineticEnzyme::LBKineticEnzyme(substrate,velocity)
  lineweaverburkAIC<-lineweaverburk$AIC
  lineweaverburkBIC<-lineweaverburk$BIC
  lineweaverburklogLike<-lineweaverburk$logLike
  lineweaverburkvmax<-lineweaverburk$vmax
  lineweaverburkkm<-lineweaverburk$km
  lineweaverburkStandardErrorvm_delta<-lineweaverburk$StandardErrorvm_delta
  lineweaverburkStandardErrorkm_delta<-lineweaverburk$StandardErrorkm_delta

  EadieScatchard<-KineticEnzyme::EAKineticEnzyme(substrate,velocity)
  EadieScatchardAIC<-EadieScatchard$AIC
  EadieScatchardBIC<-EadieScatchard$BIC
  EadieScatchardlogLike<-EadieScatchard$logLike
  EadieScatchardvmax<-EadieScatchard$vmax
  EadieScatchardkm<-EadieScatchard$km
  EadieScatchardStandardErrorvm_delta<-EadieScatchard$StandardErrorvm_delta
  EadieScatchardStandardErrorkm_delta<-EadieScatchard$StandardErrorkm_delta

  WoolfAugustinsonHofstee<-KineticEnzyme::WAHKineticEnzyme(substrate,velocity)
  WoolfAugustinsonHofsteeAIC<-WoolfAugustinsonHofstee$AIC
  WoolfAugustinsonHofsteeBIC<-WoolfAugustinsonHofstee$BIC
  WoolfAugustinsonHofsteelogLike<-WoolfAugustinsonHofstee$logLike
  WoolfAugustinsonHofsteevmax<-WoolfAugustinsonHofstee$vmax
  WoolfAugustinsonHofsteekm<-WoolfAugustinsonHofstee$km
  WoolfAugustinsonHofsteeStandardErrorvm_delta<-WoolfAugustinsonHofstee$StandardErrorvm_delta
  WoolfAugustinsonHofsteeStandardErrorkm_delta<-WoolfAugustinsonHofstee$StandardErrorkm_delta

  HanesWoolf<-KineticEnzyme::HWKineticEnzyme(substrate,velocity)
  HanesWoolfAIC<-HanesWoolf$AIC
  HanesWoolfBIC<-HanesWoolf$BIC
  HanesWoolflogLike<-HanesWoolf$logLike
  HanesWoolfvmax<-HanesWoolf$vmax
  HanesWoolfkm<-HanesWoolf$km
  HanesWoolfStandardErrorvm_delta<-HanesWoolf$StandardErrorvm_delta
  HanesWoolfStandardErrorkm_delta<-HanesWoolf$StandardErrorkm_delta


  vAIC<-c(michaelismentenAIC,lineweaverburkAIC,EadieScatchardAIC,WoolfAugustinsonHofsteeAIC,HanesWoolfAIC)
  vAICcolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorAIC <- vAIC[1]
  menorAICcolname <- vAICcolnames[1]
  for(i in seq_along(vAIC)){
    if(vAIC[i] < menorAIC){
      menorAIC = vAIC[i]
      menorAICcolname = vAICcolnames[i]
    }
  }

  vBIC<-c(michaelismentenBIC,lineweaverburkBIC,EadieScatchardBIC,WoolfAugustinsonHofsteeBIC,HanesWoolfBIC)
  vBICcolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorBIC <- vBIC[1]
  menorBICcolname <- vBICcolnames[1]
  for(i in seq_along(vBIC)){
    if(vBIC[i] < menorBIC){
      menorBIC = vBIC[i]
      menorBICcolname = vBICcolnames[i]
    }
  }

  vlogLike<-c(michaelismentenlogLike,lineweaverburklogLike,EadieScatchardlogLike,WoolfAugustinsonHofsteelogLike,HanesWoolflogLike)
  vlogLikecolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  mayorlogLike <- max(vlogLike)
  mayorlogLikecolname <- vlogLikecolnames[which.max(vlogLike)]

  StandardErrorvm_delta<-c(michaelismentenStandardErrorvm_delta,lineweaverburkStandardErrorvm_delta,EadieScatchardStandardErrorvm_delta,WoolfAugustinsonHofsteeStandardErrorvm_delta,HanesWoolfStandardErrorvm_delta)
  StandardErrorvm_deltacolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorStandardErrorvm_delta <- StandardErrorvm_delta[1]
  menorStandardErrorvm_deltacolname <- StandardErrorvm_deltacolnames[1]
  for(i in seq_along(StandardErrorvm_delta)){
    if(StandardErrorvm_delta[i] < menorStandardErrorvm_delta){
      menorStandardErrorvm_delta = StandardErrorvm_delta[i]
      menorStandardErrorvm_deltacolname = StandardErrorvm_deltacolnames[i]
    }
  }

  StandardErrorkm_delta<-c( michaelismentenStandardErrorkm_delta,lineweaverburkStandardErrorkm_delta,EadieScatchardStandardErrorkm_delta,WoolfAugustinsonHofsteeStandardErrorkm_delta,HanesWoolfStandardErrorkm_delta)
  StandardErrorkm_deltacolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorStandardErrorkm_delta <- StandardErrorkm_delta[1]
  menorStandardErrorkm_deltacolname <- StandardErrorkm_deltacolnames[1]
  for(i in seq_along(StandardErrorkm_delta)){
    if(StandardErrorkm_delta[i] < menorStandardErrorkm_delta){
      menorStandardErrorkm_delta = StandardErrorkm_delta[i]
      menorStandardErrorkm_deltacolname = StandardErrorkm_deltacolnames[i]
    }
  }



  AIC <- c(michaelismentenAIC,lineweaverburkAIC,EadieScatchardAIC,WoolfAugustinsonHofsteeAIC,HanesWoolfAIC)
  BIC <- c(michaelismentenBIC,lineweaverburkBIC,EadieScatchardBIC,WoolfAugustinsonHofsteeBIC,HanesWoolfBIC)
  logLike <- c(michaelismentenlogLike,lineweaverburklogLike,EadieScatchardlogLike,WoolfAugustinsonHofsteelogLike,HanesWoolflogLike)
  vmax <- c(michaelismentenvmax,lineweaverburkvmax,EadieScatchardvmax,WoolfAugustinsonHofsteevmax,HanesWoolfvmax)
  km <- c(michaelismentenkm,lineweaverburkkm,EadieScatchardkm,WoolfAugustinsonHofsteekm,HanesWoolfkm)
  StandardErrorvm_delta <- c(michaelismentenStandardErrorvm_delta,lineweaverburkStandardErrorvm_delta,EadieScatchardStandardErrorvm_delta,WoolfAugustinsonHofsteeStandardErrorvm_delta,HanesWoolfStandardErrorvm_delta)
  StandardErrorkm_delta <- c(michaelismentenStandardErrorkm_delta,lineweaverburkStandardErrorkm_delta,EadieScatchardStandardErrorkm_delta,WoolfAugustinsonHofsteeStandardErrorkm_delta,HanesWoolfStandardErrorkm_delta)


  df <- matrix(c(AIC, BIC,logLike, vmax, km, StandardErrorvm_delta, StandardErrorkm_delta), nrow = 5, ncol = 7, byrow = FALSE)
  rownames(df) <- c('Michaelis-Menten','Lineweaver-Burk', 'Eadie-Scatchard',
                    'Woolf-Augustinson-Hofstee', 'Hanes-Woolf')
  colnames(df) <- c('AIC','BIC','LogLikelihood','vmax','km','StandardErrorvm_delta','StandardErrorkm_delta')

  #¿Which model is better?
  cat("The smallest value of the AIC is:", menorAIC, "corresponding to the model of", menorAICcolname, "\n")
  cat("The smallest value of the BIC is:", menorBIC, "corresponding to the model of", menorBICcolname,"\n")
  cat("The biggest maximum likelihood value is:", mayorlogLike, "corresponding to the model of", mayorlogLikecolname,"\n")
  cat("The model with the least standard error in delta method of vm is:", menorStandardErrorvm_delta, "corresponding to the model of", menorStandardErrorvm_deltacolname,"\n")
  cat("The model with the smallest standard error in delta method of km is:", menorStandardErrorkm_delta, "corresponding to the model of", menorStandardErrorkm_deltacolname,"\n")

  return(df)

}

