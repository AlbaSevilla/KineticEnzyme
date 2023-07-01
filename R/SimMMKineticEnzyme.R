#' @title SimMMKineticEnzyme: Simulate Samples for the M-M Equation
#'
#' @description This function generates samples for the Monod-Michel (M-M) equation, which describes the
#' relationship between the concentration of a substrate and the rate of an enzymatic reaction.
#'
#' @param n Number of points to generate. If a value greater than zero is specified,
#'   a linear sequence of length "n" will be generated within the specified range.
#'   If zero is specified, a sequence with an increment of "sinc" will be generated within the specified range.
#' @param vm Maximum velocity (Vmax) value in the M-M equation.
#' @param km Michaelis-Menten constant (Km) value in the M-M equation.
#' @param srange Range of substrate (S) values for which the samples will be generated.
#'   Specified as a vector with the minimum and maximum values of the range.
#' @param sinc Increment value for the linear sequence of substrate values. Applicable only when "n" is zero.
#' @param sigma Standard deviation of the random noise added to the generated samples.
#' @param het Logical value indicating whether to introduce heteroscedasticity in the samples. Default is TRUE.
#' @return A data frame with the generated samples, including the substrate (S) values, the observed velocity (v),
#'   and the model-predicted velocity (vmod) based on the M-M equation.
#' @examples
#' data<-SimMMKineticEnzyme(n=20,vm=1,km=0.3)
#' plot(data$s,data$v)
#'@encoding UTF-8
#'@export SimMMKineticEnzyme

# simula muestras para la ecuación de M-M
SimMMKineticEnzyme<-function(n=0,vm,km,srange=c(0,1), sinc=0.1,sigma=0.06,het=TRUE){

  if (n>0) {x<-seq(from=srange[[1]],to=srange[[2]],length.out=n)}
  else     {x<-seq(from=srange[[1]],to=srange[[2]],by=sinc)}
  lx<-length(x)

  vmod<-vm*x/(km+x)
  z<-rnorm(lx,0,1)

  if(het){
    v<-vmod*(1+z*sigma)
  }else {
    v<-vmod+(z*sigma)
  }

  df<-data.frame(s=x,v=v,vmod)

  result=df
  return(result)
}
