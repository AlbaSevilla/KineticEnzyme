% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimINHKineticEnzyme.R
\encoding{UTF-8}
\name{SimINHKineticEnzyme}
\alias{SimINHKineticEnzyme}
\title{SimINHKineticEnzyme: Simulate Samples for the M-M Equation for inhibition data.}
\usage{
SimINHKineticEnzyme(
  n = 0,
  vm = 1,
  ks = 0.3,
  ki = 1,
  srange = c(0, 1),
  sinc = 0.1,
  inh = c(1, 2, 3),
  inhtype = "com",
  sigma = 0.06,
  het = TRUE
)
}
\arguments{
\item{n}{Specifies the number of points to generate in the input range srange. If n is greater than 0, a sequence of n equidistant points in the range is generated. If n is 0, a sequence is generated using a spacing defined by the sinc parameter.}

\item{vm}{is the maximum rate of the enzymatic reaction.}

\item{ks}{Is the Michaelis-Menten constant.}

\item{ki}{Is the inhibition constant.}

\item{srange}{Specifies the range of input values for which the inhibit function will be calculated. It is a vector of two elements, the first represents the initial value and the second the final value of the range.}

\item{sinc}{Is the spacing between points in the generated sequence if n is 0.}

\item{inh}{Specifies the inhibition levels. It is a vector that contains the inhibition values to be considered.}

\item{inhtype}{Specifies the type of inhibition to use. It can be "com" (competitive), "aco" (antagonistic) or "noc" (non-competitive).}

\item{sigma}{This is the standard deviation used to introduce variability into the generated data.}

\item{het}{A boolean value indicating whether to generate heterogeneous data (with variability) or not. If TRUE, a variability is applied to the generated data. description}
}
\description{
This function generates samples for the Michaelis Menten equation for inhibition cases.
}
\examples{
data <- SimINHKineticEnzyme(n = 10, vm = 1, ks = 0.3, ki = 1, srange = c(0, 1), sinc = 0.01, inh = c(0, 1, 2), inhtype = "com", sigma = 0.06, het = TRUE)
plot(data$s,data$v)
}
