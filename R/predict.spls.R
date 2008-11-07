
# return fit or coefficients

"predict.spls" <-
function( object, newx, type = c("fit","coefficient"), ... )
{
    # newx: matrix of predictors
    # type = "fit" or "coefficient"
    
    type <- match.arg(type)    
    betahat <- object$betahat
    x <- object$x
    
    if ( type=="fit" )
    {        
        if ( missing(newx) )
        {
            pred <- x %*% betahat + matrix(1,nrow(x),1) %*% object$mu
        } else
        {
            newx <- scale( newx, object$meanx, object$normx )
            pred <- newx %*% betahat + matrix(1,nrow(newx),1) %*% object$mu
        }
    }
    if ( type=="coefficient" ) { pred <- betahat }
    
    invisible(pred)
}

"coef.spls" <-
function( object, ... )
{
    predict.spls( object, type="coefficient", ... )
}
