# Zheyi SHEN, s2258945

# A function that fits linear models in QR decomposition method:
linmod <- function(formula, dat){
    # Given a specified linear model formula and a data frame, 
    # it returns a list of class "linmod", in which the common results of lm funtion are provided.
    # Those items would be calculated first, and then be composed in a list to return.     
    if (is.data.frame(dat) == FALSE) stop("Please use Data Frame!")
    
    y <-  model.frame(formula, dat)[[1]]  ## vector of response variable 
    yname <- all.vars(formula)[1]
    
    X <- model.matrix(formula, dat) ## model matrix (response excluded)
    n <- nrow(X) ## number of observations
    p <- ncol(X) ## number of parameters
    qrx <- qr(X) ## QR decomposition of model matrix
    beta <- backsolve(qr.R(qrx), qr.qty(qrx, y)[1:p]) ## R^{-1} Q^T y
    names(beta) <- colnames(X) ## Name coefficients by corresponding variables
    
    mu <- drop(X %*% beta) ## vector of fitted values
    SSR <- sum((y - mu)^2) ## Sum of Squared Residuals
    sigma2 <- SSR/(n-p) ## Estimated error/response variance
    sigma <- sqrt(sigma2) ## Estimated standard deviation
    
    # In QR decomposition method, beta_hat = R^{-1} Q^T y
    # which is a linear transformation of y.
    # Then Var(beta_hat) = R^{-1} Q^T y {R^{-1} Q^T} ^T
    # = R^{-1} Q^T Q R^{-T} sigma2 = R^{-1}R^{-T} sigma^2
    # as Q is orthogonal matrix having Q^{T}Q = I (identity matrix)
    # That is, R Var(beta_hat) = R^{-T}*sigma2
    # R is an upper triangular matrix having R,
    # we could use backsolve(R, R^{-T} sigma^2) to solve for Var(beta_hat)
    RT_inv <- forwardsolve(t(qr.R(qrx)),diag(p))  ## forwardsolve R^T R^{-T} = diag(p) for R^{-T}
    V <- backsolve(qr.R(qrx),RT_inv) * sigma2 ## estimated covariance matrix of beta
    colnames(V) <- rownames(V) <- colnames(X)  ## Name the var-cov matrix
    
    flev <- list() ## initializing the list of levels
    flev_label <- c() ## initializing the vector of labels
    
    Xname <- all.vars(formula)[-1]   # explanatory variable names
    for (item in Xname){
        if (is.factor(dat[[item]]) == TRUE){
            flev[[item]] <- levels(dat[[item]])
            flev_label <- c(flev_label, item)
        }
    }
    names(flev) <- flev_label    # named list flev
 
    # Compose everything informative as a list, set it as class "linmod" and return.
    linmod_summary <- list(beta = beta,
                   V = V,
                   mu = mu,
                   y = y,
                   yname = yname,
                   formula = formula,
                   flev = flev,
                   sigma = sigma)
    class(linmod_summary) <- "linmod"
    return(linmod_summary)
}



# A print method function:
print.linmod<- function(x,...){
    # Input: an object of class "linmod"
    # Output: A compact report including model formula, parametres estimates, standard deviations
    se <- sqrt(diag(x$V))  ## Standard error of estimated parametres
    report <- cbind(x$beta, se)    ## Put up parameters with standard errors
    colnames(report) <- c("Estimate", "s.e.")   ## name the columns  
    print(x$formula)    ## print formula 
    cat("\n")     ## skip a line
    print(report) 
}



# A plot method function:
plot.linmod <- function(x, ...){
    # Input: an object of class 'linmod'.
    # Output: a scatterplot showing residuals against fitted values.
    # A dashed horizontal line indicates zero residual.
    resid <- x$y - x$mu
    plot(x = x$mu, y = resid, 
         xlab = "Fitted values", ylab = "Residuals")    ## plot model residuals against fitted values
    abline(h = 0, col = 'red', lty = 'dashed') ## Zero residual indicator line
    lines(lowess(x=x$mu,y=resid))  ## smooth curve
    mtext(text="Residuals vs Fitted") ## serve title
}



# A (more resilient) predict method function:
predict.linmod <- function(x, newdata, ...){
    # Input: an object of class 'linmod' + a data frame containing predictor variables values
    # Output: a vector of predictions for newdata variables
    # Error message: when levels are not supplied along with factors.
    
    # If the newdata doesn't contain response variable,
    # just add dummy response data to it
    if (is.null(newdata[[x$yname]])){
        newdata[[x$yname]]<-sample(0:1,nrow(newdata),TRUE)
    } ## randomly generating 0-1 response variable
    
    
    raise_error <- NULL ## initial error message
    for (item in names(x$flev)){
        # For factor variables provided as character,
        # convert them into factor.
        if (is.character(newdata[[item]])){
            newdata[[item]] <- as.factor(newdata[[item]])   
        } 
        
        # check if newdata contains levels that are not included in flev
        if (all(levels(newdata[[item]]) %in% x$flev[[item]])){
            levels(newdata[[item]]) <- x$flev[[item]] 
        } 
        # Then check that for all factor variables,  
        # levels are matched between newdata and regression model.
        else{
            raise_error <- "Error: unmatched levels."
            break
        } 
        
    }   
    
    # Raise error if there exits any, or predict the model:
    if (is.null(raise_error)){
        X_new <- model.matrix(x$formula, data=newdata, xlev=x$flev) ##model matrix for newdata
        mu_new <- drop(X_new %*% (x$beta)) ## Apply estimated parameters to newdata
        return(mu_new) ## returns the predicted y
    }   
    else{ 
        print(raise_error)
    } 
} 
