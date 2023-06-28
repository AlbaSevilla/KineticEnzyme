#' @title WAHKineticEnzyme: Woolf Augustinson Hofstee Model Fitting
#' @description
#' This function fits the Woolf Augustinson Hofstee model to the given data.
#' The model is a non-linear regression model commonly used in enzyme kinetics
#' to describe the relationship between substrate concentration and reaction velocity.
#' The model is given by V = (Vmax * [S]) / (Km + [S]), where V is the reaction velocity,
#' Vmax is the maximum velocity, [S] is the substrate concentration, and Km is the Michaelis-Menten constant.
#' The function takes substrate and velocity data as input and performs the model fitting.
#' It also provides options to remove outliers from the data and perform in-depth analysis and diagnostics.
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
#' @import lmtest
#' @import utils
#' @import car
#' @import lmtest
#' @import carData
#' @return A data frame containing the estimated parameters and other analysis results.
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' WAHKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE,deepening=TRUE)
#' @export WAHKineticEnzyme
#' @encoding UTF-8

WAHKineticEnzyme <- function(substrate,velocity,removeoutliers=FALSE,deepening=FALSE) {
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
  ############ Definimos el modelo de Woolf-Augustinson-Hofstee #############
  ###########################################################################
  wah_model <- function(x, Vmax, Km) {
    -Km*x + Vmax
  }

  #Convertimos los datos
  wah_data <- data.frame(x = data$velocity/data$substrate, y = data$velocity)
  head(wah_data)

  ##########################################################################
  ############### PARÁMETROS ###############################################
  ##########################################################################

  # Estimamos los parámetros usando regresión lineal
  wah_lm <- lm(y ~ x, data = wah_data)
  wah_lm

  #Ajuste de polinomio de 2º grado
  y = data$velocity
  x= data$velocity/data$substrate
  x_2= (data$velocity/data$substrate)^2
  wah_lm2 <- lm(y ~ I(x)+I(x_2))
  summary(wah_lm2)
  coef2grado<-summary(wah_lm2)$coefficients[3,4]
  if(coef2grado<0.05){
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there is evidence to reject the null hypothesis and conclude that the coefficient of the second degree is significantly different from zero. \n")
  } else {
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there isn't evidence to reject the null hypothesis and conclude that the coefficient of the second degree isn't significantly different from zero. \n")
  }


  # Obtenemos los parámetros estimados
  menos_km <- coef(wah_lm)["x"]
  vmax <- coef(wah_lm)["(Intercept)"]
  menos_km
  vmax

  ##########################################################################
  ################# VALORES PREDICHOS POR EL MODELO ########################
  ##########################################################################

  # Calculamos los valores predichos
  wah_data$ypredicho <- predict(wah_lm)
  head(wah_data)


  if(deepening){
    ##########################################################################
    ########## COEFICIENTE DE CORRELACIÓN ####################################
    ##########################################################################

    #Calcular el coeficiente de correlación
    wah_lm_summary <- summary(wah_lm)
    wah_lm_errors <- summary(wah_lm)$coef[,2]
    std_error_vm<-wah_lm_errors[1]
    std_error_km<-wah_lm_errors[2]
    R_squared <- wah_lm_summary$r.squared

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
    ANOVA<-anova(wah_lm)
    AIC<-AIC(wah_lm)
    BIC<-BIC(wah_lm)

    #############################################################################
    ############ ANÁLISIS DE LOS RESIDUOS #######################################
    #############################################################################
    #Análisis de residuos
    ei<-residuals(wah_lm)
    plot(ei, main='Residuals')
    boxplot(ei, main="Residuals' Boxplot")
    par(mfrow=c(1, 2))
    plot(wah_lm, col='deepskyblue4', pch=19,2)
    plot(wah_lm, col='deepskyblue4', pch=19,1)
    plot(wah_lm, col='deepskyblue4', pch=19,3)
    plot(wah_lm, col='deepskyblue4', pch=19,4)
    plot(wah_lm, col='deepskyblue4', pch=19,5)

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
    x = data$velocity/data$substrate
    y = data$velocity
    plot(x,ei,ylab="Residuals",main='Homocedasticity or Heterocedasticity?')
    abline(h=0,col=2)
    # Comprobación
    pvalueBreushPagan<-bptest(wah_lm)$p.value
    if (pvalueBreushPagan < 0.05) {
      cat("The p-value of the Breush-Pagan test is less than 0.05. The null hypothesis is rejected and we conclude that the residuals are heteroscedastic. \n")
    } else {
      cat("The p-value of the Breush-Pagan test is greater than or equal to 0.05. The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals are homoscedastic. \n")
    }

    #Incorrelación o autocorrelación
    DurbinWatson<-dwtest(wah_lm,alternative='two.sided')

    #Residuos estandarizados y estudentizados
    rstandard(wah_lm)

    #Test de valores anómalos
    outlierTest(wah_lm)
  }

  #############################################################################
  ############# INTERVALOS DE CONFIANZA PARA LOS PARÁMETROS ###################
  #############################################################################

  # Calculate the confidence intervals for the estimated parameters
  ci <- confint(wah_lm)
  Vmax_ci <- ci["x", ]
  Km_ci <- ci["(Intercept)", ]

  ###############################################################################
  ########################## GRÁFICO Y BANDAS DE CONFIANZA ######################
  ###############################################################################
  x = data$velocity/data$substrate
  y = data$velocity
  #Intervalos respuesta media y predicción
  predict(wah_lm,newdata=data.frame(x=length(x)+100),
          interval='confidence',se.fit=T)
  predict(wah_lm,newdata=data.frame(x=length(x)+100),
          interval='prediction',se.fit=T)
  errorestandarconfianza<-predict(wah_lm,newdata=data.frame(x=length(x)+100),
                                  interval='confidence',se.fit=T)$se.fit
  errorestandarprediccion<-predict(wah_lm,newdata=data.frame(x=length(x)+100),
                                   interval='prediction',se.fit=T)$se.fit

  #Bandas de confianza
  x0 <- data.frame(x=seq(min(wah_data$x), max(wah_data$x), length.out=length(wah_data$x)+100))
  pred.m<-predict(wah_lm,newdata=x0,interval='confidence',se.fit=T)
  pred.p<-predict(wah_lm,newdata=x0,interval='prediction',se.fit=T)
  matplot(x0$x,cbind(pred.m$fit,pred.p$fit[,-1]),lty=c(1,2,2,3,3), col=c("seagreen3",2,2,4,4),type='l', xlab = "v/s", ylab = "v",
          main = "Woolf-Augustinson-Hofstee model",xlim=range(wah_data$x),ylim=range(wah_data$y))
  legend('topleft',c('95% Confidence Band for the mean response','95% Confidence Band for the prediction'), lty=c(2,3),col=c(2,4),bty='n')
  points(wah_data$x, wah_data$y)

  #plot(wah_data$x, wah_data$y)
  #abline(wah_lm, col="seagreen3")

  ##########################################################################################
  ############ Calcular los parámetros estimados en la escala original #####################
  ##########################################################################################
  km <- -menos_km
  vmax <- vmax
  AIC<-AIC(wah_lm)
  BIC<-BIC(wah_lm)
  logLike<-logLik(wah_lm)
  DurbinWatson<-dwtest(wah_lm,alternative='two.sided')
  wah_lm_summary <- summary(wah_lm)
  errores<-wah_lm_summary$coef[,2]
  matrizcovarianza<-vcov(wah_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(wah_lm)))
  covarianza <- matrizcovarianza[1,2]
  errorvm_delta<-errores[1]
  errorkm_delta<-errores[2]
  ANOVA<-anova(wah_lm)
  #Errores estándares
  matrizcovarianza<-vcov(wah_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(wah_lm)))
  covarianza <- matrizcovarianza[1,2]
  errores_matriz_covarianza<-sqrt(diag(vcov(wah_lm)))
  errorvm_linealizacion <- errores_matriz_covarianza[1]
  errorkm_linealizacion <- errores_matriz_covarianza[2]

  # Return the estimated parameters
  Resultados <- list(AIC = AIC, BIC = BIC,logLike=logLike, vmax=vmax, km=km,
                     StandardErrorvm_delta=errorvm_delta, StandardErrorkm_delta=errorkm_delta,
                     ANOVA,DurbinWatson,Resumen=wah_lm_summary)
  return(Resultados)
}
