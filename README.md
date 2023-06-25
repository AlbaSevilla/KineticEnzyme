# KineticEnzyme

This package provides a comprehensive solution for data analysis
related to the Michaelis Menten model, including linearization of
the data and consideration of the three main cases of enzyme
inhibition. In addition, estimates of great relevance in the
biochemical field are offered, such as the estimated parameters,
the standard errors, and the residual analysis in those cases where
it is necessary.

The package offers various functionalities to handle the data and
perform analyses. First, it prepares the data by removing rows that
contain zeros. Additionally, if specified, outliers can be removed
using Tukey's criterion. The package then proceeds to estimate the
model parameters using linear regression and calculates the predicted
values by the model.

For advanced analysis and diagnostics, the package provides further
functionalities. This includes calculating correlations, performing
model comparisons using criteria such as ANOVA, AIC, and BIC, and
analyzing the residuals for normality, homoscedasticity, and
autocorrelation. Furthermore, the package enables the calculation of
confidence intervals for the estimated parameters and allows for
visually interpreting the results by plotting the data along with
corresponding confidence bands.

In summary, the KineticEnzyme package offers a complete and cohesive
tool for analyzing data related to the Michaelis Menten model. It
provides a wide range of functionalities, facilitates interpretation,
and aids in decision making in the biochemical field.
