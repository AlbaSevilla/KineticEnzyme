#' @title HWKineticEnzyme: Hanes Woolf Model Fitting
#' @description
#'     This function fits the Hanes Woolf model to the given data by estimating the kinetic parameters
#'     of an enzyme-catalyzed reaction. The Hanes Woolf model is a linearized form of the Michaelis-Menten
#'     equation, commonly used to analyze enzyme kinetics.
#'
#'     The function takes substrate concentration data and corresponding velocity data as input. It fits
#'     the Hanes Woolf model to the data using linear regression and provides estimates of the kinetic
#'     parameters, including the maximum velocity (Vmax) and the Michaelis constant (Km).
#'
#'     Optionally, you can choose to remove outliers from the data by setting the parameter
#'     `removeoutliers` to TRUE. This can be useful to improve the accuracy of the model fitting.
#'
#'     Additionally, you can enable the `deepening` parameter to perform in-depth analysis and
#'     diagnostics, which includes residual analysis, goodness-of-fit tests, and diagnostic plots.
#'
#' @param substrate Substrate concentration data.
#' @param velocity Velocity data corresponding to the substrate concentrations.
#' @param removeoutliers Logical value indicating whether to remove outliers from the data. Default is FALSE.
#' @param deepening Logical value indicating whether to perform in-depth analysis and diagnostics. Default is FALSE.
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @import lmtest
#' @import car
#' @import lmtest
#' @import carData
#' @return A data frame containing the estimated parameters (Vmax and Km) and other analysis results.
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' HWKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE,deepening=TRUE)

#' @export HWKineticEnzyme
#' @encoding UTF-8

HWKineticEnzyme <- function(substrate,velocity,removeoutliers=FALSE,deepening=FALSE) {
  if (!requireNamespace("car", quietly = TRUE))
    install.packages("car", repos = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-2.tar.gz", type = "source")
  if (!requireNamespace("minpack.lm", quietly = TRUE))
    install.packages("minpack.lm", repos = "https://cran.r-project.org/src/contrib/minpack.lm_1.2-3.tar.gz", type = "source")
  if (!requireNamespace("lmtest", quietly = TRUE))
    install.packages("lmtest",repos="https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",type="source")
  library(car)
  library(minpack.lm)
  library(lmtest)

  data <- data.frame(substrate,velocity)
  data <- subset(data, !rowSums(data == 0) > 0)

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

  ###########################################################################
  ############ Definimos el modelo de Hanes-Woolf #######################
  ###########################################################################
  hw_model <- function(x, Vmax_inv, one_divide_vmax_inv) {
    Vmax_inv * x + one_divide_vmax_inv
  }

  # Convertimos los datos
  hw_data <- data.frame(x = data$substrate, y = data$substrate/data$velocity)
  head(hw_data)

  ##########################################################################
  ############### PARÁMETROS ###############################################
  ##########################################################################

  # Estimamos los parámetros usando regresión lineal
  hw_lm <- lm(y ~ x, data = hw_data)
  hw_lm

  #Ajuste de polinomio de 2º grado
  y = data$substrate/data$velocity
  x = data$substrate
  x_2= (data$substrate)^2
  hw_lm2 <- lm(y ~ I(x)+I(x_2))
  summary(hw_lm2)
  coef2grado<-summary(hw_lm2)$coefficients[3,4]
  if(coef2grado<0.05){
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there is evidence to reject the null hypothesis and conclude that the coefficient of the second degree is significantly different from zero. \n")
  } else {
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there isn't evidence to reject the null hypothesis and conclude that the coefficient of the second degree isn't significantly different from zero. \n")
  }

  # Obtenemos los parametros estimados
  one_divide_vm <- coef(hw_lm)["x"]
  km_divide_vm <- coef(hw_lm)["(Intercept)"]
  one_divide_vm
  km_divide_vm

  ##########################################################################
  ################# VALORES PREDICHOS POR EL MODELO ########################
  ##########################################################################

  # Calculamos los valores predichos
  hw_data$ypredicho <- predict(hw_lm)
  head(hw_data)

  if(deepening){
    ##########################################################################
    ########## COEFICIENTE DE CORRELACIÓN ####################################
    ##########################################################################

    #Calcular el coeficiente de correlación
    hw_lm_summary <- summary(hw_lm)
    errores<-hw_lm_summary$coef[,2]
    errorvm<-errores[1]
    errorkm<-errores[2]

    R_squared <- hw_lm_summary$r.squared
    if(R_squared == 1){
      cat("Perfect Positive Correlation:", R_squared, "\n")
    }
    if(R_squared > 0.9 && R_squared < 1){
      cat("High Positive Correlation:", R_squared, "\n")
    }
    if(R_squared > 0.5 && R_squared < 0.9){
      cat("Low Positive Correlation:", R_squared, "\n")
    }
    if(R_squared > 0 && R_squared < 0.5){
      cat("Low Positive Correlation:", R_squared, "\n")
    }
    if(R_squared == 0 ){
      cat("No Correlation:", R_squared, "\n")
    }
    if(R_squared < 0 && R_squared > -0.5){
      cat("Low Negative Correlation:", R_squared, "\n")
    }
    if(R_squared < -0.5 && R_squared > -0.9){
      cat("Low Negative Correlation:", R_squared, "\n")
    }
    if(R_squared < -0.9 && R_squared > -1){
      cat("High Negative Correlation:", R_squared, "\n")
    }
    if(R_squared == -1){
      cat("Perfect Negative Correlation:", R_squared, "\n")
    }

    #############################################################################
    ################ CRITERIOS DE COMPARACIÓN DE MODELOS ########################
    #############################################################################
    ANOVA<-anova(hw_lm)
    AIC<-AIC(hw_lm)
    BIC<-BIC(hw_lm)

    #############################################################################
    ############ ANÁLISIS DE LOS RESIDUOS #######################################
    #############################################################################
    #Análisis de residuos
    ei<-residuals(hw_lm)
    plot(ei, main='Residuals')
    boxplot(ei, main="Residuals' Boxplot")
    par(mfrow=c(1, 2))
    plot(hw_lm, col='deepskyblue4', pch=19,2)
    plot(hw_lm, col='deepskyblue4', pch=19,1)
    plot(hw_lm, col='deepskyblue4', pch=19,3)
    plot(hw_lm, col='deepskyblue4', pch=19,4)
    plot(hw_lm, col='deepskyblue4', pch=19,5)

    # Obtener información sobre el diagrama de caja y bigotes
    bxp <- boxplot.stats(ei)

    #Los residuos que provocan la no normalidad son :
    out <- boxplot.stats(ei)$out
    indices_outliers <- which(ei %in% out)
    residuosnonormalidad <- data.frame(Valor = out, Indice = indices_outliers)

    print("The residual that alters normality is")
    print(residuosnonormalidad)

    # Seleccionar los valores que no son atípicos
    #ei <- subset(ei, ei >= bxp$stats[1] & ei <= bxp$stats[5])

    #Histograma de los residuos
    hist(ei,xlab='Residuals',freq=FALSE,main="'Residuals' Histogram")
    #Funcion de densidad estimada
    lines(density(ei),col=4)

    # Obtener el tamaño de la muestra
    n <- length(data$substrate)

    # Realizar la prueba de normalidad si TestNormality es verdadero
    # Si la muestra es menor que 50, utilizar el test de Shapiro-Wilk
    if (n<50){

    # Realizar la prueba de Shapiro-Wilk y obtener el p-valor
    ShapiroTestNormality <- shapiro.test(ei)
    p_valor_shapiro <- ShapiroTestNormality$p.value

    # Evaluar el p-valor y mostrar un mensaje en consecuencia
      if (p_valor_shapiro < 0.05) {
        cat("The p-value of the Shapiro Wilks Normality test is", p_valor_shapiro, ". The null hypothesis is rejected and we conclude that the residuals do not follow a normal distribution. \n")
      } else {
        cat("The p-value of the Shapiro Wilks Normality test is", p_valor_shapiro, ". The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals follow a normal distribution. \n")
      }

    }
    # Si la muestra es mayor o igual que 50, utilizar el test de Kolmogorov-Smirnov
    else{

      # Realizar la prueba de Kolmogorov-Smirnov y obtener el p-valor
      KolmogorovSmirnovNormality <- ks.test(ei, pnorm, mean=0, sd=sd(ei))
      p_valor_kolmogorov <- KolmogorovSmirnovNormality$p.value

      # Evaluate the p-value and display a message accordingly
      if (p_valor_kolmogorov < 0.05) {
        cat("The p-value of the Kolmogorov Smirnov Normality test is", p_valor_kolmogorov, ". The null hypothesis is rejected and we conclude that the residuals do not follow a normal distribution. \n")
      } else {
        cat("The p-value of the Kolmogorov Smirnov Normality test is", p_valor_kolmogorov, ". The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals follow a normal distribution. \n")
      }
    }

    #Heterocedasticidad?
    x = data$substrate
    y = data$substrate/data$velocity
    plot(x,ei,ylab="Residuals",main='Homocedasticity or Heterocedasticity?')
    abline(h=0,col=2)
    # Comprobación
    pvalueBreushPagan<-bptest(hw_lm)$p.value
    if (pvalueBreushPagan < 0.05) {
      cat("The p-value of the Breush-Pagan test is less than 0.05. The null hypothesis is rejected and we conclude that the residuals are heteroscedastic. \n")
    } else {
      cat("The p-value of the Breush-Pagan test is greater than or equal to 0.05. The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals are homoscedastic. \n")
    }

    #Incorrelación o autocorrelación
    DurbinWatson<-dwtest(hw_lm,alternative='two.sided')

    #Residuos estandarizados y estudentizados
    rstandard(hw_lm)

    #Test de valores anómalos
    outlierTest(hw_lm)
  }

  #############################################################################
  ############# INTERVALOS DE CONFIANZA PARA LOS PARÁMETROS ###################
  #############################################################################
  # Calculate the confidence intervals for the estimated parameters
  ci <- confint(hw_lm)
  Vmax_inv_ci <- ci["x", ]
  one_divide_vmax_inv_ci <- ci["(Intercept)", ]

  ###############################################################################
  ########################## GRÁFICO Y BANDAS DE CONFIANZA ######################
  ###############################################################################
  x = data$substrate
  y = data$substrate/data$velocity
  #Intervalos respuesta media y predicción
  predict(hw_lm,newdata=data.frame(x=length(x)+100),
          interval='confidence',se.fit=T)
  predict(hw_lm,newdata=data.frame(x=length(x)+100),
          interval='prediction',se.fit=T)
  errorestandarconfianza<-predict(hw_lm,newdata=data.frame(x=length(x)+100),
                                  interval='confidence',se.fit=T)$se.fit
  errorestandarprediccion<-predict(hw_lm,newdata=data.frame(x=length(x)+100),
                                   interval='prediction',se.fit=T)$se.fit

  #Bandas de confianza
  x0 <- data.frame(x=seq(min(hw_data$x), max(hw_data$x), length.out=length(hw_data$x)+100))
  pred.m<-predict(hw_lm,newdata=x0,interval='confidence',se.fit=T)
  pred.p<-predict(hw_lm,newdata=x0,interval='prediction',se.fit=T)
  matplot(x0$x,cbind(pred.m$fit,pred.p$fit[,-1]),lty=c(1,2,2,3,3), col=c("seagreen3",2,2,4,4),type='l', xlab = "s", ylab = "s/v",
          main = "Hanes-Woolf model",xlim=range(hw_data$x),ylim=range(hw_data$y))
  legend('topleft',c('95% Confidence Band for the mean response','95% Confidence Band for the prediction'), lty=c(2,3),col=c(2,4),bty='n')
  points(hw_data$x, hw_data$y)

  #plot(hw_data$x, hw_data$y)
  #abline(hw_lm, col="seagreen3")

  ##########################################################################################
  ############ Calcular los parámetros estimados en la escala original #####################
  ##########################################################################################
  vmax <- 1/one_divide_vm
  km <- km_divide_vm*vmax
  AIC<-AIC(hw_lm)
  BIC<-BIC(hw_lm)
  logLike<-logLik(hw_lm)
  hw_lm_summary <- summary(hw_lm)
  matrizcovarianza<-vcov(hw_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(hw_lm)))
  covarianza <- matrizcovarianza[1,2]
  errores<-hw_lm_summary$coef[,2]
  errorvm_delta<-(((errores[2])/(abs(one_divide_vm)))*(abs(vmax)))/(1)-2*covarianza
  errorkm_delta<- (((errores[1]/km_divide_vm^2)-(errorvm_delta^2/vmax^2))*km^2)-2*covarianza
  ANOVA<-anova(hw_lm)
  DurbinWatson<-dwtest(hw_lm,alternative='two.sided')
  #Errores estándares
  matrizcovarianza<-vcov(hw_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(hw_lm)))
  covarianza <- matrizcovarianza[1,2]
  errores_matriz_covarianza<-sqrt(diag(vcov(hw_lm)))
  errorvm_linealizacion <- errores_matriz_covarianza[1]
  errorkm_linealizacion <- errores_matriz_covarianza[2]

  # Return the estimated parameters
  Resultados <- list(AIC = AIC, BIC = BIC,logLike=logLike, vmax=vmax, km=km,
                     StandardErrorvm_delta=errorvm_delta, StandardErrorkm_delta=errorkm_delta,
                     StandardErrorvm_linealizacion=errorvm_linealizacion,
                     StandardErrorkm_linealizacion=errorkm_linealizacion,
                     ANOVA,DurbinWatson,Resumen=hw_lm_summary)
  return(Resultados)
}

