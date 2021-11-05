# linmod_Regression-QRdecomposition
Functions to fitting linear regression using QR decomposition approach.
See code file P3s2258945.R

### The aim of the project
Write your own functions for fitting linear models (using the QR decomposition approach), printing a simple summary of the fit, plotting a default residual plot and predicting the expected response, given new values of the predictors. The educational aim of the practical is to deepen your understanding of linear models and of classes and matrix computation in R.

### Function Lists
1. linmod(formula,dat):
  Takes a model formula, *formula*, specifying a linear model, and a data frame containing the corresponding data, *dat*. 
  Estimates the specified linear model using the QR decomposition of the model matrix approach. 
  Returns an object of class "linmod", a list containing the following elements:
    1. *beta* the vector of least squares parameter estimates for the model. The elements of beta should be named to identify the model components that they relate to (as in R’s *lm* function). (Hint: *colnames*, *names*.)
    2. *V* the estimated covariance matrix of the least squares estimators.
    3. *mu* the vector of expected values of the response variable according to the estimated model (also known as ‘fitted values’).
    4. *y* the vector containing the response variable. y <- model.frame(formula,dat)[[1]] will obtain this (model frames are beyond the scope of this course, so you can treat this line of code as a ‘black box’).
    5. *yname* the name of the response variable, obtainable via yname <- all.vars(formula)[1].
    6. *formula* the model formula.
    7. *flev* a named list. For each factor variable in dat, flev contains an item which is the vector of levels of the factor - the item’s name should be the name of the factor variable. You need flev to predict properly with factors.
    8. *sigma* the estimated standard deviation, σˆ, of the response (and residuals).
    
    
2. print.linmod(x,...) 
    A print method function. 
    x is an object of class "linmod". 
    The function should give the model formula defining the model and report the parameter estimates and their standard deviations (hint: *cbind*, *colnames*).
    
3. plot.linmod — a plot method function. 

4. predict.linmod(x,newdata) — a predict method function. 
