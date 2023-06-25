#' @title SimINHKineticEnzyme: Simulate Samples for the M-M Equation for inhibition data.
#'
#' @description This function generates samples for the Michaelis Menten equation for inhibition cases.
#' @param n Specifies the number of points to generate in the input range srange. If n is greater than 0, a sequence of n equidistant points in the range is generated. If n is 0, a sequence is generated using a spacing defined by the sinc parameter.
#' @param vm is the maximum rate of the enzymatic reaction.
#' @param ks Is the Michaelis-Menten constant.
#' @param ki Is the inhibition constant.
#' @param srange Specifies the range of input values for which the inhibit function will be calculated. It is a vector of two elements, the first represents the initial value and the second the final value of the range.
#' @param sinc Is the spacing between points in the generated sequence if n is 0.
#' @param inh Specifies the inhibition levels. It is a vector that contains the inhibition values to be considered.
#' @param inhtype Specifies the type of inhibition to use. It can be "com" (competitive), "aco" (antagonistic) or "noc" (non-competitive).
#' @param sigma This is the standard deviation used to introduce variability into the generated data.
#' @param het A boolean value indicating whether to generate heterogeneous data (with variability) or not. If TRUE, a variability is applied to the generated data. description
#' @examples
#' data <- SimINHKineticEnzyme(n = 10, vm = 1, ks = 0.3, ki = 1, srange = c(0, 1), sinc = 0.01, inh = c(0, 1, 2), inhtype = "com", sigma = 0.06, het = TRUE)
#' plot(data$s,data$v)
#'@encoding UTF-8
#'@export SimINHKineticEnzyme

SimINHKineticEnzyme<-function(n=0,vm=1,ks=0.3,ki=1,srange=c(0,1), sinc=0.1,inh=c(0,1,2),inhtype="com", sigma=0.06,het=TRUE){

  if (n>0) {x<-seq(from=srange[[1]],to=srange[[2]],length.out=n)}
  else     {x<-seq(from=srange[[1]],to=srange[[2]],by=sinc)}
  lx<-length(x)

  tipos<-as.factor(c("com","aco","noc"))
  if (which(tipos==inhtype,tipos)==0) stop("Tipo incorrecto")


  vmod<-c()
  inhdat<-c()
  inhlevels<-nlevels(as.factor(inh))
  validtype=FALSE

  if(inhtype=="com") {
     for(i in 1:inhlevels){vmod<-c(vmod, (vm * x / (ks*(1+inh[i]/ki) + x)))
     inhdat<-c(inhdat,inh[i])
     }
     validtype=TRUE
  }
  if(inhtype=="aco"){
    for(i in 1:inhlevels){vmod<-c(vmod, (vm * x / (ks+x(1+inh[i]/ki))))
    inhdat<-c(inhdat,inh[i])
    }
    validtype=TRUE
  }
  if(inhtype=="noc"){
    for(i in 1:inhlevels){vmod<-c(vmod, (vm * x / (ks*(1+inh[i]/ki) + x*(1+inh[i]/ki))))
    inhdat<-c(inhdat,inh[i])
    }
    validtype=TRUE
  }
  if (!validtype) stop("Tipo no especificado. Tipos válidos: 'com', 'aco','noc' \n")

  z<-rnorm(lx,0,1)

  if(het){
     v<-vmod*(1+z*sigma)
  }else {
     v<-vmod+(z*sigma)
  }

  df<-data.frame(s=x,v=v,inh=inhdat,vmod)

  result=df
  return(result)
}
