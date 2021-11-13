linmod <- function(formula, dat){
    ## This is a function to fit linear models as lm() in R does, using QR decomposition approach.
    ## Input: formula - a specified linear model formula, e.g. y ~ x1 + x2 ; dat - a data frame containing corresponding data
    ## Output: re - an object of class 'linmod', a list containing following elements:
    ##  beta - the vector of least squares parameter estimates      
    ##  V - estimated covariance matrix 
    ##  mu - the vector of expected values      
    ##  y - the vector of response variables
    ##  yname - the name of response variables      
    ##  formula - model formula
    ##  flev - a named list for factor variables and their levels     
    ##  fitted_value - a vector of fitted values of response variable using parameters estimates
    ##  residuals - a vector of residuals between the true values of response variable and its fitted values
    ##  sigma - estimated standard deviation
    
    y <- model.frame(formula, dat)[[1]]   ## the vector of response variable
    yname <- all.vars(formula)[1]   ## the name of response variable
    
    xname <- all.vars(formula)[-1]  ## the name(s) of explanatory variable(s)
    x <- model.matrix(~. , dat[xname])    ## model matrix
    
    p <- ncol(x)    ## number of parameters
    n <- nrow(x)    ## number of observations
    qrx <- qr(x)    ## get QR decomposition of model matrix
    beta <- backsolve(qr.R(qrx), qr.qty(qrx, y)[1 : p])   ## R^{-1} Q^T y
    names(beta) <- colnames(x)    ## name beta to identify the model components
    
    y_fitted <- beta %*% t(x)   ## the vector of 'fitted values'
    mu <- mean(y_fitted)    ## expected value of 'fitted values'
    names(mu) <- 'fitted'   ## name to identify fitted values
    
    residuals <- y - y_fitted     ## residuals
    var_residuals <- as.numeric(residuals %*% t(residuals) / (n - p))     ## the variance of residuals
    inv_R <- backsolve(qr.R(qrx), diag(p))    ## the inverse of R
    V <- (inv_R %*% t(inv_R)) * var_residuals    ## estimated covariance matrix of least square estimators, R^{-1} Q^{-T} sigma^2
    colnames(V) <- rownames(V) <- colnames(x)   ## name the rows and columns of estimated covariance matrix
    
    flev <- NULL    ## initial list flev
    for (i in xname){
        if (is.factor(dat[[i]]) == TRUE){
            flev[[i]] <- levels(dat[[i]])
        }
    } ## a loop to find factor variables and store their levels in list flev
    
    sigma <- c(sqrt(var_residuals))   ## the estimated standard deviation of response and residuals
    names(sigma) <- 'residuals'     ## name sigma to identify response and residuals
    
    re <- list(beta = beta, V = V, mu = mu, y = y, yname = yname, formula = formula, flev = flev, sigma = sigma, fitted_values = drop(y_fitted), residuals = drop(residuals)) 
    class(re) <- 'linmod'   ## return to an object of class 'linmod'
    return(re)
    
} ## Linear Regression Function


print.linmod <- function(x){
    ##  a print method function for x, an object of class 'linmod'
    ##  Output: model formula and a matrix contains parameters estimates and their standard errors 
    
    est_bind <- cbind(x$beta, sqrt(diag(x$V)))    ## calculate standard deviations and bind with parameters estimates 
    colnames(est_bind) <- c('Estimate', 's.e.')   ## name the columns of bound matrix    
    print(x$formula)    ## print formula 
    cat("\n")     ## print blank line
    print(est_bind)     ## print bound matrix
    
} ## Print Method Function


plot.linmod <- function(x){
    ##  a plot method function for x, an object of class 'linmod'
    ##  Output: a scatter plot of residuals of fitted values with a dashed line where residual is 0 
    
    plot(x = x$fitted_values, y = x$residuals, xlab = 'fitted values', ylab = 'residuals')    ## plot model residuals against fitted values
    abline(h = 0, col = 'red', lty = 'dashed')  ## dashed horizontal line
    
} ## Plot Method Function


predict.linmod <- function(x, newdata){
    ## a predict method for x, an object of class 'linmod', with newdata, a data frame containing values of predictor variables
    ## Output: an error message if newdata contains wrong factor levels, or a vector of predictions of response variable based on newdata 
    
    error_message <- NULL   ## initial error_message
    
    for (i in names(x$flev)){
        
        if (is.character(newdata[[i]]) == TRUE){
            newdata[[i]] <- as.factor(newdata[[i]])   
        } ## convert factor levels provided as character strings into factor class
        
        ## following condition statements check if newdata contains levels that are not included in flev
        if (all(levels(newdata[[i]]) %in% x$flev[[i]]) == TRUE){
            levels(newdata[[i]]) <- x$flev[[i]] 
        } ## if it doesn't, ensure that corresponding factor variables in newdata have same levels as those in fitting model
        else{
            error_message <- 'Error: level out of range'
            break
        } ## if it does, set up an error message and break the loop
        
    }   
    
    if (is.null(error_message) == TRUE){
        xname <- all.vars(x$formula)[-1]    ## the names of explanatory variable(s) 
        x_new <- model.matrix(~. , newdata[xname])    ## model matrix from newdata
        y_predict <- drop(x$beta %*% t(x_new))   ## calculate the vector of predictions of response variable
        print(y_predict)
    }   ## if newdata doesn't contain levels that are not included in flev, use fitting model to predict
    else{ 
        print(error_message)
    } ## if newdata contains levels that are not included in flev, print error message
    
} ## Predict Method Function