#' @title LBKineticEnzyme: Lineweaver-Burk Model Fitting
#' @description
#' This function fits the Lineweaver-Burk model to the given data.
#' The Lineweaver-Burk model is a linearized form of the Michaelis-Menten equation
#' used to analyze enzyme kinetics.
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
#'
#' @import grDevices
#' @import graphics
#' @import lmtest
#' @import stats
#' @import utils
#' @import car
#' @import lmtest
#' @import carData
#'
#' @examples
#' s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
#' v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
#' data <- cbind(s,v)
#' data<-data.frame(data)
#' LBKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE,deepening=TRUE)
#' @return A data frame containing the estimated parameters and other analysis results.
#'
#' @export LBKineticEnzyme
#' @encoding UTF-8

LBKineticEnzyme <- function(substrate,velocity,removeoutliers=FALSE,deepening=FALSE){
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
  }

  ###########################################################################
  ############ Definimos el modelo de lineweaver burk #######################
  ###########################################################################
  lb_model <- function(x, Vmax_inv, Km_Vmax_inv) {
    Vmax_inv * x + Km_Vmax_inv
  }

  # Convertimos los datos
  lb_data <- data.frame(x = 1/data$substrate, y = 1/data$velocity)
  lb_data



  ##########################################################################
  ############### PARÁMETROS ###############################################
  ##########################################################################

  # Estimamos los parámetros usando regresión lineal
  lb_lm <- lm(y ~ x, data = lb_data)

  #Ajuste de polinomio de 2º grado
  y = 1/data$velocity
  x= 1/data$substrate
  x_2= 1/(data$substrate^2)
  lb_lm2 <- lm(y ~ I(x)+I(x_2))
  summary(lb_lm2)
  coef2grado<-summary(lb_lm2)$coefficients[3,4]
  if(coef2grado<0.05){
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there is evidence to reject the null hypothesis and conclude that the coefficient of the second degree is significantly different from zero. \n")
  } else {
    cat("The p-value of the second grade coefficient is", coef2grado, ". So there isn't evidence to reject the null hypothesis and conclude that the coefficient of the second degree isn't significantly different from zero. \n")
  }

  #Errores estándares
  matrizcovarianza<-vcov(lb_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(lb_lm)))

  # Obtenemos los parámetros estimados
  km_vm <- coef(lb_lm)["x"]
  one_divide_vm <- coef(lb_lm)["(Intercept)"]

  ##########################################################################
  ################# VALORES PREDICHOS POR EL MODELO ########################
  ##########################################################################

  # Calculamos los valores predichos de 1/velocity
  lb_data$ypredicho <- predict(lb_lm)
  head(lb_data)

  if(deepening){
    ##########################################################################
    ########## COEFICIENTE DE CORRELACIÓN ####################################
    ##########################################################################

    #Calcular el coeficiente de correlación
    lb_lm_summary <- summary(lb_lm)
    errores<-lb_lm_summary$coef[,2]
    errorvm<-errores[1]
    errorkm<-errores[2]
    R_squared <- lb_lm_summary$r.squared
    cat("Correlation between variables is :", R_squared, "\n")

    #############################################################################
    ################ CRITERIOS DE COMPARACIÓN DE MODELOS ########################
    #############################################################################
    ANOVA<-anova(lb_lm)
    AIC<-AIC(lb_lm)
    BIC<-BIC(lb_lm)

    #############################################################################
    ############ ANÁLISIS DE LOS Residuals #######################################
    #############################################################################
    #Análisis de Residuals
    ei<-residuals(lb_lm)
    boxplot(ei, main="Residuals' Boxplot",col="orange")
    par(mfrow=c(1, 2))
    plot(lb_lm, col='deepskyblue4', pch=19,2)
    plot(lb_lm, col='deepskyblue4', pch=19,1)
    plot(lb_lm, col='deepskyblue4', pch=19,3)
    plot(lb_lm, col='deepskyblue4', pch=19,4)
    plot(lb_lm, col='deepskyblue4', pch=19,5)

    # Obtener información sobre el diagrama de caja y bigotes
    bxp <- boxplot.stats(ei)

    #Los Residuals que provocan la no normalidad son :
    out <- boxplot.stats(ei)$out
    indices_outliers <- which(ei %in% out)
    Residualsnonormalidad <- data.frame(Valor = out, Indice = indices_outliers)

    print("That residue that alters normality is")
    print(Residualsnonormalidad)

    # Seleccionar los valores que no son atípicos
    #ei <- subset(ei, ei >= bxp$stats[1] & ei <= bxp$stats[5])

    #Histograma de los Residuals
    hist(ei,xlab='Residuals',freq=FALSE,main="Residuals' Histogram")
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
      p_value_kolmogorov <- KolmogorovSmirnovNormality$p.value

      # Evaluar el p-valor y mostrar un mensaje en consecuencia
      if (p_value_kolmogorov < 0.05) {
        cat("The p-value of the Kolmogorov Smirnov Normality test is", p_value_kolmogorov, ". The null hypothesis is rejected and we conclude that the residuals do not follow a normal distribution. \n")
      } else {
        cat("The p-value of the Kolmogorov Smirnov Normality test is", p_value_kolmogorov, ". The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals follow a normal distribution. \n")
      }
    }

    #Heterocedasticidad?
    x = 1/data$substrate
    y = 1/data$velocity
    plot(x,ei,ylab="Residuals",main='Homocedasticity or Heterocedasticity?')
    abline(h=0,col=2)
    # Verification
    pvalueBreushPagan<-bptest(lb_lm)$p.value
    if (pvalueBreushPagan < 0.05) {
      cat("The p-value of the Breush-Pagan test is less than 0.05. The null hypothesis is rejected and we conclude that the residuals are heteroscedastic. \n")
    } else {
      cat("The p-value of the Breush-Pagan test is greater than or equal to 0.05. The null hypothesis is not rejected and we conclude that there is no evidence to rule out that the residuals are homoscedastic. \n")
    }

    #Incorrelación o autocorrelación
    DurbinWatson<-dwtest(lb_lm,alternative='two.sided')

    #Residuals estandarizados y estudentizados
    rstandard(lb_lm)

    #Test de valores anómalos
    outlierTest(lb_lm)
  }
  #############################################################################
  ############# INTERVALOS DE CONFIANZA PARA LOS PARÁMETROS ###################
  #############################################################################

  # Calculate the confidence intervals for the estimated parameters
  ci <- confint(lb_lm)
  Vmax_inv_ci <- ci["x", ]
  one_divide_vm_ci <- ci["(Intercept)", ]

  ###############################################################################
  ########################## GRÁFICO Y BANDAS DE CONFIANZA ######################
  ###############################################################################
  x = 1/data$substrate
  y = 1/data$velocity
  #Intervalos respuesta media y predicción
  predict(lb_lm,newdata=data.frame(x=length(x)+100),
          interval='confidence',se.fit=T)
  predict(lb_lm,newdata=data.frame(x=length(x)+100),
          interval='prediction',se.fit=T)
  errorestandarconfianza<-predict(lb_lm,newdata=data.frame(x=length(x)+100),
                                  interval='confidence',se.fit=T)$se.fit
  errorestandarprediccion<-predict(lb_lm,newdata=data.frame(x=length(x)+100),
                                   interval='prediction',se.fit=T)$se.fit

  #Bandas de confianza
  x0 <- data.frame(x=seq(min(lb_data$x), max(lb_data$x), length.out=length(lb_data$x)+100))
  pred.m <- predict(lb_lm, newdata=x0, interval='confidence', se.fit=T)
  pred.p <- predict(lb_lm, newdata=x0, interval='prediction', se.fit=T)
  matplot(x0$x, cbind(pred.m$fit, pred.p$fit[,-1]), lty=c(1,2,2,3,3), col=c("seagreen3",2,2,4,4), type='l', xlab = "1/substrate", ylab = "1/velocity", main = "Lineweaver-Burk model")
  legend('topleft', c('95% Confidence Band for the mean response','95% Confidence Band for the prediction'), lty=c(2,3), col=c(2,4), bty='n')
  points(lb_data$x, lb_data$y)


  #plot(lb_data$x, lb_data$y)
  #abline(lb_lm, col="seagreen3")

  ##########################################################################################
  ############ Calcular los parámetros estimados en la escala original #####################
  ##########################################################################################
  vmax <- 1/one_divide_vm
  km <- km_vm*vmax
  ANOVA<-anova(lb_lm)
  DurbinWatson<-dwtest(lb_lm,alternative='two.sided')
  AIC<-AIC(lb_lm)
  BIC<-BIC(lb_lm)
  logLike<-logLik(lb_lm)
  lb_lm_summary <- summary(lb_lm)
  matrizcovarianza<-vcov(lb_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(lb_lm)))
  covarianza <- matrizcovarianza[1,2]
  errores<-lb_lm_summary$coef[,2]
  errorvm_delta<-(((errores[1])/(abs(one_divide_vm)))*(abs(vmax)))/(1)
  errorkm_delta<-((errores[2]/km_vm^2)-(errorvm_delta^2/vmax^2))*km^2
  #Errores estándares
  matrizcovarianza<-vcov(lb_lm)
  errores_matriz_covarianza<-sqrt(diag(vcov(lb_lm)))
  covarianza <- matrizcovarianza[1,2]
  errores_matriz_covarianza<-sqrt(diag(vcov(lb_lm)))
  errorvm_linealizacion <- errores_matriz_covarianza[1]
  errorkm_linealizacion <- errores_matriz_covarianza[2]


  # MÉTODO DELTA
  #Inferencia sobre vm
  se_vm <- matrizcovarianza[1,1]
  errorvm_delta <- (1/(one_divide_vm)^2)*se_vm
  #Inferencia sobre km
  varianza_residual <- (summary(lb_lm)$sigma)**2
  x = 1/data$substrate
  n <- length(x)
  beta_0 <- one_divide_vm
  beta_1 <- km_vm
  media_x <- mean(x)
  varianza_x <- var(x)
  varianza_km <- varianza_residual*((1/n)+((media_x)^2/((n-1)*varianza_x))*(-beta_1/beta_0^2))+(varianza_residual/((n-1)*varianza_x))*(1/beta_0^2)-2*((-media_x*varianza_residual^2)/((n-1)*varianza_x))*(-beta_1/beta_0)*(1/beta_1)
  errorkm_delta <- sqrt(varianza_km)

  # Return the estimated parameters
  Resultados <- list(AIC = AIC, BIC = BIC,logLike=logLike, vmax=vmax, km=km,
                     StandardErrorvm_delta=errorvm_delta,
                     StandardErrorkm_delta=errorkm_delta,
                     ANOVA,DurbinWatson,Resumen=lb_lm_summary)
  return(Resultados)
}
