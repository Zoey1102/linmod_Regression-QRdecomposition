##Overview:
##------------------------------------------------------------------------------------------------------------------------
##The aim of this code is to create a function that can be use for fitting linear model by using the 
##QR decomposition. Overall,there are four task that this code will performed:
##
##1. Calculate and display the model summary after fitting the linear model by using linmod function
##
##2. Calculate and display the regression model parameter and standard error of its parameter using print.linmod function
##
##3. Display the residual vs predicted value plot using plot.linmod function
##
##4. Predict the value of response variable or corresponding error if a new data is given by using predict.linmod function
##   after through some data processing (conversion of factor variable and adding response column containing dummy variable)
##------------------------------------------------------------------------------------------------------------------------


linmod <- function(formula,dat){
    
    ##function to estimates the specified linear model using the QR decomposition of the model matrix approach.
    
    ##input: formula = regression model formula, dat = the data used to build regression model
    
    ##output: 
    ##list containing this following element:
    
    ##beta = vector of estimated regression parameter; V = covariance matrix;
    ##mu = vector of fitted values or predicted values of the response variable; y = vector containing the response variable data;
    ##yname = response variable name; formula = model formula;
    ##flev = a named list containing factor and level factor; sigma = the estimated standard deviation
    
    ##Notation:
    ##Let A be a matrix size nxn. We use this following notation in our code:
    ##A^{-1} = inverse of A; A^T= transpose of A; 
    ##Q = orthogonal matrix from QR decomposition calculation;
    ##R = upper triangular matrix from QR decomposition calculation
    ##X = design matrix(matrix that contains data of expanatory variable) 
    
    
    X <- model.matrix(formula,dat) ##model matrix
    qrx <- qr(X) ##get QR decomposition of model matrix
    
    
    y <- model.frame(formula,dat)[[1]] ##response variable data
    yname <- all.vars(formula)[1] ##response variable name
    p <- ncol(X) ##number of parameters
    n <- nrow(X) ##number of observations
    
    
    beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p]) ##regression parameter (R^{-1} Q^T y)
    names(beta) <- colnames(X) ##label beta column based on the column name of X
    
    
    mu <- X %*% beta ##fitted values or predicted values
    mu <- drop(mu) ##return fitted value result as a vector
    sse <- sum((mu-y)^2) ##calculate sum squared error
    sigma <- sqrt(sse/(n-p)) ##estimated standard deviation
    
    
    ##Recall that in QR decomposition R is an upper triangular matrix. Therefore, R^{T} (transpose of R) is an lower triangular matrix. Thus, we can calculate
    ##R^{-T} by solving R^{T}x=I using forwardsolve and calculate R^{-1}R^{-T} by solving Rx=R^{-T} because R is an upper triangular matrix. 
    
    
    inv_Rt <- forwardsolve(t(qr.R(qrx)),diag(p)) ##calculate R^{-T}
    inv_RtR <- backsolve(qr.R(qrx),inv_Rt) ##calculate R^{-1}R^{-T}
    V <- sigma^2 * inv_RtR ##covariance matrix (sigma^2R^{-1}R^{-T})
    
    
    rownames(V) <- colnames(V) <- colnames(X) ##name the row and columns of covariance matrix by the parameter name
    
    
    ##In this step, we identify factor variable in the data and check weather they included in the formula or not. 
    ##The reason we do this is to make the investigation much faster from the variable that satisfied this criteria
    ##in order to investigate the level factor instead searching from all variable in the data.
    
    
    flev <- list() ##initialize factor level list
    flev_var <- c() ##initialize vector of factor variable name
    
    
    ##identify factor variable
    factors <- names(which(sapply(dat,is.factor))) ##finding all factor variable in the data
    factor_index <- which(all.vars(formula) %in% factors) 
    factor_name <- all.vars(formula)[factor_index] ##factor variable in the regression formula
    
    
    ##find level factor and the corresponding factor variable name that was used in the formula
    for (i in factor_name){
        flev <- c(flev,list(levels(dat[,i])))
        flev_var <- c(flev_var, i)
    } 
    
    names(flev) <- flev_var ##label level factor vector with factor variable name
    
    ##store all required calculation result in a list and return object of class linmod
    model_summary <- list(beta=beta, V=V, mu=mu, y=y,yname=yname,formula=formula, flev=flev, sigma=sigma)
    class(model_summary) <- 'linmod'
    return(model_summary)
}


print.linmod <- function(x){
    
    ##function to display the model formula along with the estimated parameter and its standard error
    
    ##input: x = object of class linmod
    
    ##output: model summary containing regression formula, estimated value of each parameter with its standard error
    
    
    ##To calculate the estimated standard error of regression coefficient we can do it by calculating the square root
    ##of each entries in the diagonal of covariance matrix.
    
    
    se <- sapply(diag(x$V),sqrt) ##standard error of regression parameter
    
    ##store estimated regression parameter and its standard error in a column matrix
    model_report <- cbind(x$beta,se)
    colnames(model_report) <- c('Estimate','s.e.')##give name to the matrix column
    
    ##displaying the regression model formula along with the parameter and standard error of the parameter
    print(x$formula)
    cat('\n')
    print(model_report)
}


plot.linmod <- function(x){
    
    ##function to plot residual vs fitted value
    
    ##input: x = object of class linmod
    
    ##output: residual vs fitted value plot
    
    
    resid <- x$y-x$mu ##residual
    fitted <- x$mu ##fitted value
    
    plot(x=fitted,y=resid,xlab='Fitted values',ylab='Residuals')##create plot
    lines(lowess(x=fitted,y=resid),col="red")##smooth line of residual pattern
    abline(h=0,lty=3)##dashed line (residual equal zero)
    mtext(text="Residuals vs Fitted") ## plot title
}


predict.linmod <- function(x,newdata){
    
    ##function to estimate the response variable value given newdata
    
    ##input: x = object of class linmod , newdata = data without response variable
    
    ##ouput:input error or vector containing the predicted value or response variable by using newdata
    
    
    ##Recall that newdata is not containing response variable. So, we have to add the dummy response
    ##variable(0 and 1) so when we specified the formula from x class model matrix can detect it in
    ##newdata and perform the correct calculation that we desired. Because in calculation of predicted
    ##values using newdata didn't require response variable, so it is safe to set any values in the 
    ##dummy response variable.
    
    
    ##give dummy response variable to newdata if this data doesn't contain response variable
    if (is.null(newdata[[x$yname]])==TRUE){
        newdata[[x$yname]]<-sample(0:1,nrow(newdata),replace=TRUE)
    }
    
    ##During newdata processing we should reconsider error such as there exist an explanatory variable in the formula that wasn't included
    ##in newdata and a level factor has a different level than in newdata (if we included a factor variable). So, it would be better if 
    ##the error output displayed from predict.linmod function to make the user now what the cause of error. 
    
    ##checking error associated with variable name that was input in newdata
    ii <- which(!all.vars(x$formula)[-1] %in% colnames(newdata)) ##find explanatory variable in formula that wasn't included in newdata
    
    if (length(ii)!=0){
        output <- 'Error : Explanatory Variable In Formula Not Found In Input Data!'
    }
    
    ##ignored error if no factor level
    else if(length(x$flev)==0){
        output <- 'No Error' 
    } 
    
    else{
        for (i in names(x$flev)){
            
            ##convert a variable in newdata that suppose to be a factor variable in dat
            if (is.factor(newdata[[i]])==FALSE){
                newdata[[i]] <- factor(newdata[[i]]) 
            }
            
            ##checking error associated with level factor problem
            jj <- which(!levels(newdata[[i]]) %in% x$flev[[i]]) ##find level factor in newdata that wasn't included in flev
            
            if (length(jj)!=0){
                output <- 'Error : Level Factor Is Difference!'
                break
            }
            
            else{
                output <-'No Error'
            }
            
        }
    }
    
    ##When constructing a model.matrix for newdata, we specified the level factor list (either empty or not) of factor variable in dat in the 
    ##formula to make the model.matrix for newdata do the same factor encoding (if there exist factor variable in the formula) like in dat, 
    ##even though the factor is same but incomplete. Empty list of level factor indicates no factor variable included in the formula, which we 
    ##use and specify if the formula doesn't contain any factor variable for this code.
    
    
    ##display prediction result (error or predicted values)
    if (output != 'No Error'){
        result <- output ##cause of error
    }
    
    else{
        X_new <- model.matrix(x$formula, data=newdata, xlev=x$flev) ##model matrix for newdata
        mu_new <- X_new %*% x$beta ##predicted value for newdata
        mu_new <- drop(mu_new) ##return predicted value for newdata as a vector
        result <- mu_new 
    }
    return(result)
}