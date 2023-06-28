#' @title EAKineticEnzyme: Eadie Scatchard Model Fitting
#' @description This function fits the Eadie Scatchard model to the given data
#' and provides a comprehensive analysis of the model. The Eadie Scatchard model
#' is a kinetic model used to analyze enzyme-substrate interactions and determine
#' important parameters such as the maximum velocity (Vmax) and the Michaelis constant (Km).

#' The function takes two input vectors, substrate and velocity,
#' representing the substrate concentrations and the corresponding
#' velocity measurements, respectively. It also includes optional
#' parameters: removeoutliers and deepening to control outlier
#' removal and in-depth analysis, respectively.
#'
#' @param substrate
#' Sustrate data. A numeric vector containing the substrate concentrations.
#' @param velocity
#' Velocity data. A numeric vector containing the corresponding reaction velocities.
#' @param removeoutliers
#' Logical value indicating whether to remove outliers from the data.
#' If set to TRUE, the function will remove outliers using a suitable method
#' (e.g., the Hampel identifier or the Tukey's fences method) before fitting the model.
#' Default is FALSE.
#' @param deepening
#' Logical value indicating whether to perform in-depth analysis and diagnostics.
#' If set to TRUE, the function will provide additional diagnostic plots and statistics
#' to assess the goodness of fit and identify potential issues with the model.
#' Default is FALSE.
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @import lmtest
#' @import car
#' @import lmtest
#' @import carData
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' EAKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE,deepening=TRUE)
#' @return A data frame containing the estimated parameters and other analysis results.
#' @encoding UTF-8
#' @export EAKineticEnzyme


EAKineticEnzyme <- function(substrate,velocity,removeoutliers=FALSE,deepening=FALSE) {
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
  ############ Definimos el modelo de Eadie-Scatchard #######################
  ###########################################################################
  EadieS_model <- function(x, Vmax, Km) {
    -(1/Km)*x + (Vmax/Km)
  }

  # Convertimos los datos
  EadieS_data <- data.frame(x = data$velocity, y = data$velocity/data$substrate)
  head(EadieS_data)

  ##########################################################################
  ############### PARÁMETROS ###############################################
  ##########################################################################

  # Estimamos los parámetros usando regresión lineal
  EadieS_lm <- lm(y ~ x, data = EadieS_data)
  EadieS_lm

  #Ajuste de polinomio de 2º grado
  y = data$velocity/data$substrate
  x = data$velocity
  x_2 = data$velocity^2
  EadieS_lm2 <- lm(y ~ I(x)+I(x_2))
  summary(EadieS_lm2)
  coef2grado<-summary(EadieS_lm2)$coefficients[3,4]
  if(coef2grado<0.05){
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there is evidence to reject the null hypothesis and conclude that the coefficient of the second degree is significantly different from zero. \n")
  } else {
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there isn't evidence to reject the null hypothesis and conclude that the coefficient of the second degree isn't significantly different from zero. \n")
  }

  # Obtenemos los parámetros estimados
  Km_est <- coef(EadieS_lm)["x"]
  vmax_km <- coef(EadieS_lm)["(Intercept)"]
  Km_est #-1/km
  vmax_km #vm/km

  ##########################################################################
  ################# VALORES PREDICHOS POR EL MODELO ########################
  ##########################################################################

  # Calculamos los valores predichos
  EadieS_data$ypredicho <- predict(EadieS_lm)
  head(EadieS_data)

  if(deepening){
    ##########################################################################
    ########## COEFICIENTE DE CORRELACIÓN ####################################
    ##########################################################################

    #Calcular el coeficiente de correlación
    EadieS_lm_summary <- summary(EadieS_lm)
    errores<-EadieS_lm_summary$coef[,2]
    errorvm<-errores[1]
    errorkm<-errores[2]

    R_squared <- EadieS_lm_summary$r.squared
    if(R_squared == 1){
      cat("Perfect Positive Correlation \n")
    }
    if(R_squared > 0.9 && R_squared < 1){
      cat("High Positive Correlation \n")
    }
    if(R_squared > 0.5 && R_squared < 0.9){
      cat("Low Positive Correlation \n")
    }
    if(R_squared > 0 && R_squared < 0.5){
      cat("Low Positive Correlation \n")
    }
    if(R_squared == 0 ){
      cat("No Correlation \n")
    }
    if(R_squared < 0 && R_squared > -0.5){
      cat("Low Negative Correlation \n")
    }
    if(R_squared < -0.5 && R_squared > -0.9){
      cat("Low Negative Correlation \n")
    }
    if(R_squared < -0.9 && R_squared > -1){
      cat("High Negative Correlation \n")
    }
    if(R_squared == -1){
      cat("Perfect Negative Correlation \n")
    }

    #############################################################################
    ################ CRITERIOS DE COMPARACIÓN DE MODELOS ########################
    #############################################################################
    ANOVA<-anova(EadieS_lm)
    AIC<-AIC(EadieS_lm)
    BIC<-BIC(EadieS_lm)

    #############################################################################
    ############ ANÁLISIS DE LOS RESIDUOS #######################################
    #############################################################################
    #Análisis de residuos
    ei<-residuals(EadieS_lm)
    plot(ei, main='Residuals')
    boxplot(ei, main="Residuals' Boxplot")
    par(mfrow=c(1, 2))
    plot(EadieS_lm, col='deepskyblue4', pch=19,2)
    plot(EadieS_lm, col='deepskyblue4', pch=19,1)
    plot(EadieS_lm, col='deepskyblue4', pch=19,3)
    plot(EadieS_lm, col='deepskyblue4', pch=19,4)
    plot(EadieS_lm, col='deepskyblue4', pch=19,5)

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
    x = data$velocity
    y = data$velocity/data$substrate
    plot(x,ei,ylab="Residuals",main='¿Homocedasticity o Heterocedasticity?')
    abline(h=0,col=2)
    # Comprobación
    pvalueBreushPagan<-bptest(EadieS_lm)$p.value
    if (pvalueBreushPagan < 0.05) {
      cat("The p-value of the Breush-Pagan test is less than 0.05. The null hypothesis is rejected and we conclude that the residuals are heteroscedastic. \n")
    } else {
      cat("The p-value of the Breush-Pagan test is greater than or equal to 0.05. The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals are homoscedastic. \n")
    }

    #Incorrelación o autocorrelación
    DurbinWatson<-dwtest(EadieS_lm,alternative='two.sided')

    #Residuos estandarizados y estudentizados
    rstandard(EadieS_lm)

    #Test de valores anómalos
    outlierTest(EadieS_lm)
  }

  #############################################################################
  ############# INTERVALOS DE CONFIANZA PARA LOS PARÁMETROS ###################
  #############################################################################
  # Calculate the confidence intervals for the estimated parameters
  ci <- confint(EadieS_lm)
  Vmax_ci <- ci["x", ]
  Km_ci <- ci["(Intercept)", ]

  ###############################################################################
  ########################## GRÁFICO Y BANDAS DE CONFIANZA ######################
  ###############################################################################
  x = data$velocity
  y = data$velocity/data$substrate
  #Intervalos respuesta media y predicción
  predict(EadieS_lm,newdata=data.frame(x=length(x)+100),
          interval='confidence',se.fit=T)
  predict(EadieS_lm,newdata=data.frame(x=length(x)+100),
          interval='prediction',se.fit=T)
  errorestandarconfianza<-predict(EadieS_lm,newdata=data.frame(x=length(x)+100),
                                  interval='confidence',se.fit=T)$se.fit
  errorestandarprediccion<-predict(EadieS_lm,newdata=data.frame(x=length(x)+100),
                                   interval='prediction',se.fit=T)$se.fit

  #Bandas de confianza
  x0 <- data.frame(x=seq(min(EadieS_data$x), max(EadieS_data$x), length.out=length(EadieS_data$x)+100))
  pred.m<-predict(EadieS_lm,newdata=x0,interval='confidence',se.fit=T)
  pred.p<-predict(EadieS_lm,newdata=x0,interval='prediction',se.fit=T)
  matplot(x0$x,cbind(pred.m$fit,pred.p$fit[,-1]),lty=c(1,2,2,3,3), col=c("seagreen3",2,2,4,4),type='l', xlab = "v", ylab = "v/s",
          main = "Eadie-Scatchard model",xlim=range(EadieS_data$x),ylim=range(EadieS_data$y))
  legend('topleft',c('95% Confidence Band for the mean response','95% Confidence Band for the prediction'), lty=c(2,3),col=c(2,4),bty='n')
  points(EadieS_data$x, EadieS_data$y)

  # Calculate the estimated parameters in original scale
  km <- -1/Km_est #-1/km
  vmax <- vmax_km*km #vm/km
  AIC<-AIC(EadieS_lm)
  BIC<-BIC(EadieS_lm)
  logLike<-logLik(EadieS_lm)
  EadieS_lm_summary <- summary(EadieS_lm)
  errores<-EadieS_lm_summary$coef[,2]
  matrizcovarianza<-vcov(EadieS_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(EadieS_lm)))
  covarianza <- matrizcovarianza[1,2]
  errorkm_delta<-((errores[2]/abs(Km_est))*abs(km))/1
  errorvm_delta<-((errores[1]/vmax_km^2)-((errorkm_delta)^2/(km^2)))*(vmax^2)
  ANOVA<-anova(EadieS_lm)
  DurbinWatson<-dwtest(EadieS_lm,alternative='two.sided')
  #Errores estándares
  matrizcovarianza<-vcov(EadieS_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(EadieS_lm)))
  covarianza <- matrizcovarianza[1,2]
  errorvm_linealizacion <- errores_matriz_covarianza[1]
  errorkm_linealizacion <- errores_matriz_covarianza[2]

  # Return the estimated parameters
  Resultados <- list(AIC = AIC, BIC = BIC,logLike=logLike, vmax=vmax, km=km,
                     StandardErrorvm_delta=errorvm_delta, StandardErrorkm_delta=errorkm_delta,
                     StandardErrorvm_linealizacion=errorvm_linealizacion,
                     StandardErrorkm_linealizacion=errorkm_linealizacion,
                     ANOVA,DurbinWatson,Resumen=EadieS_lm_summary)
  return(Resultados)
}
