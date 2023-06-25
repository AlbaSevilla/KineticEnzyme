#' @title MMKineticEnzyme: Michaelis Menten Model Fitting
#' @description
#' This function fits the Michaelis Menten model to the given data.
#' It estimates the parameters of the model using the nlsLM function from the minpack.lm library.
#' The function also calculates the standard errors of the estimated parameters and provides
#' additional analysis results such as the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC).
#'
#' @param substrate A vector or column of substrate concentrations.
#' @param velocity A vector or column of corresponding velocity measurements.
#' @param removeoutliers Logical value indicating whether to remove outliers from the data. Default is FALSE.
#'
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @import car
#' @import lmtest
#' @import minpack.lm
#' @import carData
#'
#' @return A data frame containing the following columns:
#' \item{AIC}{Akaike Information Criterion.}
#' \item{BIC}{Bayesian Information Criterion.}
#' \item{vmax}{Estimated maximum velocity parameter.}
#' \item{km}{Estimated Michaelis-Menten constant parameter.}
#' \item{StandardErrorvm}{Standard error of the estimated maximum velocity parameter.}
#' \item{StandardErrorkm}{Standard error of the estimated Michaelis-Menten constant parameter.}
#'
#' The function also plots the data and the fitted curve, and adds confidence bands to the plot.
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' MMKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE)
#' @export MMKineticEnzyme
#' @encoding UTF-8

MMKineticEnzyme <- function(substrate,velocity,removeoutliers=FALSE) {
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest",repos="https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",type="source")
  library(car)
  library(minpack.lm)
  library(lmtest)

  data <- data.frame(substrate, velocity)
  data <- subset(data, !rowSums(data == 0) > 0)
  km<-unname(summary(data$substrate)[2])
  vmax<-max(data$velocity)

  if(removeoutliers){
    ######################################################################
    #### ELIMINAMOS LOS OUTLIERS MEDIANTE EL CRITERIO DE TUKEY ###########
    ######################################################################
    boxplot(data$substrate, data$velocity, names = c("Substrate", "Velocity"),main="¿Outliers?", col = c("lightpink", "lightblue"))

    ############ SUBSTRATE ####################
    # identificar cuartiles y rango intercuartil
    Q <- quantile(substrate, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(substrate)

    up <-  Q[2]+1.5*iqr # Rango superior
    low <- Q[1]-1.5*iqr # Rango inferior

    # filtrar valores dentro del rango
    data<- subset(data, substrate > low & substrate < up)

    ############# VELOCITY ####################
    # identificar cuartiles y rango intercuartil
    Q <- quantile(velocity, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(velocity)

    up <-  Q[2]+1.5*iqr # Rango superior
    low <- Q[1]-1.5*iqr # Rango inferior

    # filtrar valores dentro del rango
    data<- subset(data,velocity > low & velocity < up)
    dim(data)
  }

  #######################################################################
  ##### DEFINIMOS EL MODELO DE MICHAELIS MENTEN #########################
  #######################################################################

  mm_model <- function(x, vmax, Km) {
    vmax * x / (Km + x)
  }

  # Estimate the parameters using the nlsLM function from the minpack.lm library
  fit <- nlsLM(velocity ~ mm_model(data$substrate, vmax, Km), data = data, start = list(vmax = vmax, Km = km))
  summary(fit)
  logLik(fit)

  # Get the estimated parameters
  vmax <- coef(fit)["vmax"]
  km <- coef(fit)["Km"]

  #Errores estándar de los parámetros
  std_error<-summary(fit)$coef[,2]
  errorvm<-std_error[1]
  errorkm<-std_error[2]

  # Calculate the fitted values
  vmod <- predict(fit)

  # Plot the data and the fitted curve
  plot(data$substrate, data$velocity, xlab = "S", ylab = "V",
       main = "Michaelis-Menten model", pch = 8, cex = 0.5, xlim=NULL, ylim=NULL)
  lines(data$substrate, vmod, col = "blue", lwd = 2)

  # Calculate the confidence bands
  ci <- confint(fit)
  vmax_ci <- ci["vmax", ]
  Km_ci <- ci["Km", ]

  s_seq <- seq(from = min(data$substrate), to = max(data$substrate), length.out = length(data$substrate))

  v_upper <- mm_model(s_seq, vmax_ci[2], Km_ci[1])
  v_lower <- mm_model(s_seq, vmax_ci[1], Km_ci[2])

  # Add the confidence bands to the plot
  lines(s_seq, v_upper, lty = 2, col = "red")
  lines(s_seq, v_lower, lty = 2, col = "red")

  AIC<- AIC(fit)
  BIC<-BIC(fit)
  logLike <- logLik(fit)

  errores_matriz_covarianza<-sqrt(diag(vcov(fit)))
  errorvm_linealizacion <- errores_matriz_covarianza[1]
  errorkm_linealizacion <- errores_matriz_covarianza[1]

  std_error<-summary(fit)$coef[,2]
  errorvm<-std_error[1]
  errorkm<-std_error[2]

  # Return the estimated parameters
  Resultados <- data.frame(AIC = AIC, BIC = BIC,
                           logLike=logLike,
                           vmax=vmax, km=km,
                           StandardErrorvm=errorvm, StandardErrorkm=errorkm,
                           StandardErrorvm_linealizacion=errorvm_linealizacion,
                           StandardErrorkm_linealizacion=errorkm_linealizacion)
  rownames(Resultados) <- c("Estimation")
  return(Resultados)
}
