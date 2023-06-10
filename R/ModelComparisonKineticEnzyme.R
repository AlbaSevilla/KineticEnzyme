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
#' \item{StdErrorvmax}{Standard error of the estimated maximum velocity parameter for each model.}
#' \item{StdErrorkm}{Standard error of the estimated Michaelis-Menten constant parameter for each model.}
#'
#' The row names of the data frame correspond to the names of the enzyme kinetic models.
#'
#' Additionally, the function prints the following information:
#' \item{Enzyme kinetic model with the lowest AIC value}{The enzyme kinetic model with the lowest AIC value.}
#' \item{Enzyme kinetic model with the lowest BIC value}{The enzyme kinetic model with the lowest BIC value.}
#' \item{Enzyme kinetic model with the lowest Log Likelihood value}{The enzyme kinetic model with the lowest Log Likelihood value.}
#' \item{Enzyme kinetic model with the lowest standard error of the maximum velocity (vm) parameter}{The enzyme kinetic model with the lowest standard error of the maximum velocity (vm) parameter.}
#' \item{Enzyme kinetic model with the lowest standard error of the Michaelis-Menten (km) constant parameter}{The enzyme kinetic model with the lowest standard error of the Michaelis-Menten (km) constant parameter.}
#'
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
  michaelismentenstderrorvmax<-michaelismenten$StandardErrorvm
  michaelismentenstderrorkm<-michaelismenten$StandardErrorkm


  lineweaverburk<-KineticEnzyme::LBKineticEnzyme(substrate,velocity)
  lineweaverburkAIC<-lineweaverburk$AIC
  lineweaverburkBIC<-lineweaverburk$BIC
  lineweaverburklogLike<-lineweaverburk$logLike
  lineweaverburkvmax<-lineweaverburk$vmax
  lineweaverburkkm<-lineweaverburk$km
  lineweaverburkstderrorvmax<-lineweaverburk$StandardErrorvm
  lineweaverburkstderrorkm<-lineweaverburk$StandardErrorkm

  EadieScatchard<-KineticEnzyme::EAKineticEnzyme(substrate,velocity)
  EadieScatchardAIC<-EadieScatchard$AIC
  EadieScatchardBIC<-EadieScatchard$BIC
  EadieScatchardlogLike<-EadieScatchard$logLike
  EadieScatchardvmax<-EadieScatchard$vmax
  EadieScatchardkm<-EadieScatchard$km
  EadieScatchardstderrorvmax<-EadieScatchard$StandardErrorvm
  EadieScatchardstderrorkm<-EadieScatchard$StandardErrorkm

  WoolfAugustinsonHofstee<-KineticEnzyme::WAHKineticEnzyme(substrate,velocity)
  WoolfAugustinsonHofsteeAIC<-WoolfAugustinsonHofstee$AIC
  WoolfAugustinsonHofsteeBIC<-WoolfAugustinsonHofstee$BIC
  WoolfAugustinsonHofsteelogLike<-WoolfAugustinsonHofstee$logLike
  WoolfAugustinsonHofsteevmax<-WoolfAugustinsonHofstee$vmax
  WoolfAugustinsonHofsteekm<-WoolfAugustinsonHofstee$km
  WoolfAugustinsonHofsteestderrorvmax<-WoolfAugustinsonHofstee$StandardErrorvm
  WoolfAugustinsonHofsteestderrorkm<-WoolfAugustinsonHofstee$StandardErrorkm

  HanesWoolf<-KineticEnzyme::HWKineticEnzyme(substrate,velocity)
  HanesWoolfAIC<-HanesWoolf$AIC
  HanesWoolfBIC<-HanesWoolf$BIC
  HanesWoolflogLike<-HanesWoolf$logLike
  HanesWoolfvmax<-HanesWoolf$vmax
  HanesWoolfkm<-HanesWoolf$km
  HanesWoolfstderrorvmax<-HanesWoolf$StandardErrorvm
  HanesWoolfstderrorkm<-HanesWoolf$StandardErrorkm


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

  stderrorvm<-c( michaelismentenstderrorvmax,lineweaverburkstderrorvmax,EadieScatchardstderrorvmax,WoolfAugustinsonHofsteestderrorvmax,HanesWoolfstderrorvmax)
  stderrorvmcolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorstderrorvm <- stderrorvm[1]
  menorstderrorvmcolname <- stderrorvmcolnames[1]
  for(i in seq_along(stderrorvm)){
    if(stderrorvm[i] < menorstderrorvm){
      menorstderrorvm = stderrorvm[i]
      menorstderrorvmcolname = stderrorvmcolnames[i]
    }
  }

  stderrorkm<-c( michaelismentenstderrorkm,lineweaverburkstderrorkm,EadieScatchardstderrorkm,WoolfAugustinsonHofsteestderrorkm,HanesWoolfstderrorkm)
  stderrorkmcolnames<-c("Michaelis Menten","Lineweaver Burk","Eadie Scatchard","Woolf Augustinson Hofstee","Hanes Woolf")
  menorstderrorkm <- stderrorkm[1]
  menorstderrorkmcolname <- stderrorkmcolnames[1]
  for(i in seq_along(stderrorkm)){
    if(stderrorkm[i] < menorstderrorkm){
      menorstderrorkm = stderrorkm[i]
      menorstderrorkmcolname = stderrorkmcolnames[i]
    }
  }

  AIC <- c(michaelismentenAIC,lineweaverburkAIC,EadieScatchardAIC,WoolfAugustinsonHofsteeAIC,HanesWoolfAIC)
  BIC <- c(michaelismentenBIC,lineweaverburkBIC,EadieScatchardBIC,WoolfAugustinsonHofsteeBIC,HanesWoolfBIC)
  logLike <- c(michaelismentenlogLike,lineweaverburklogLike,EadieScatchardlogLike,WoolfAugustinsonHofsteelogLike,HanesWoolflogLike)
  vmax <- c(michaelismentenvmax,lineweaverburkvmax,EadieScatchardvmax,WoolfAugustinsonHofsteevmax,HanesWoolfvmax)
  km <- c(michaelismentenkm,lineweaverburkkm,EadieScatchardkm,WoolfAugustinsonHofsteekm,HanesWoolfkm)
  stderrorvmax <- c(michaelismentenstderrorvmax,lineweaverburkstderrorvmax,EadieScatchardstderrorvmax,WoolfAugustinsonHofsteestderrorvmax,HanesWoolfstderrorvmax)
  stderrorkm <- c(michaelismentenstderrorkm,lineweaverburkstderrorkm,EadieScatchardstderrorkm,WoolfAugustinsonHofsteestderrorkm,HanesWoolfstderrorkm)

  df <- matrix(c(AIC, BIC,logLike, vmax, km, stderrorvmax, stderrorkm), nrow = 5, ncol = 7, byrow = FALSE)
  rownames(df) <- c('Michaelis-Menten','Lineweaver-Burk', 'Eadie-Scatchard',
                    'Woolf-Augustinson-Hofstee', 'Hanes-Woolf')
  colnames(df) <- c('AIC','BIC','LogLikelihood','vmax','km','StdErrorvmax','StdErrorkm')

  #¿Which model is better?
  cat("The smallest value of the AIC is:", menorAIC, "corresponding to the model of", menorAICcolname, "\n")
  cat("The smallest value of the BIC is:", menorBIC, "corresponding to the model of", menorBICcolname,"\n")
  cat("The biggest maximum likelihood value is:", mayorlogLike, "corresponding to the model of", mayorlogLikecolname,"\n")
  cat("The model with the least standard error of vm is:", menorstderrorvm, "corresponding to the model of", menorstderrorvmcolname,"\n")
  cat("The model with the smallest standard error of km is:", menorstderrorkm, "corresponding to the model of", menorstderrorkmcolname,"\n")

  return(df)

}

