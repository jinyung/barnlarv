# calculate shape score s as in Drake and Klingenberg (2008)
# pred argument specify whether to use raw Y value or predicted Y value...
# ... useful for drawing the regression line
shapesc <- function(model, pred = FALSE) {
  if (pred)
    Y <- predict(model)
  else
    Y <- model$y
  B <- model$coefficients[2, ] # 1st row intercept, 2nd row coeff
  Y %*% B %*% (t(B) %*% B)^(-0.5)
}
