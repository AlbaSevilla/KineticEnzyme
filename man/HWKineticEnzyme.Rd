% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HWKineticEnzyme.R
\encoding{UTF-8}
\name{HWKineticEnzyme}
\alias{HWKineticEnzyme}
\title{HWKineticEnzyme: Hanes Woolf Model Fitting}
\usage{
HWKineticEnzyme(substrate, velocity, removeoutliers = FALSE, deepening = FALSE)
}
\arguments{
\item{substrate}{Substrate concentration data.}

\item{velocity}{Velocity data corresponding to the substrate concentrations.}

\item{removeoutliers}{Logical value indicating whether to remove outliers from the data. Default is FALSE.}

\item{deepening}{Logical value indicating whether to perform in-depth analysis and diagnostics. Default is FALSE.}
}
\value{
A data frame containing the estimated parameters (Vmax and Km) and other analysis results.
}
\description{
This function fits the Hanes Woolf model to the given data by estimating the kinetic parameters
    of an enzyme-catalyzed reaction. The Hanes Woolf model is a linearized form of the Michaelis-Menten
    equation, commonly used to analyze enzyme kinetics.

    The function takes substrate concentration data and corresponding velocity data as input. It fits
    the Hanes Woolf model to the data using linear regression and provides estimates of the kinetic
    parameters, including the maximum velocity (Vmax) and the Michaelis constant (Km).

    Optionally, you can choose to remove outliers from the data by setting the parameter
    `removeoutliers` to TRUE. This can be useful to improve the accuracy of the model fitting.

    Additionally, you can enable the `deepening` parameter to perform in-depth analysis and
    diagnostics, which includes residual analysis, goodness-of-fit tests, and diagnostic plots.
}
\examples{
s<-c(0.5,1.0,2.4,4.4,7.0,8.0,13.1,17.6,20.2,30.0)
v<-c(2.12,3.33,4.90,6.22,7.01,7.51,7.87,8.44,8.76,8.98)
data <- cbind(s,v)
data<-data.frame(data)
HWKineticEnzyme(substrate=data$s,velocity=data$v,removeoutliers=TRUE,deepening=TRUE)
}
