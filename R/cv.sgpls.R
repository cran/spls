
"cv.sgpls" <-
function( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
        br=TRUE, ftype='iden' )
{
    y.u <- unique(y)
    ngroups <- length(y.u)
    
    if ( ngroups >= 2 )
    {
        if ( ngroups == 2 )
        {
            cv <- cv.sgpls.binary( x=x, y=y, fold=fold, K=K, eta=eta,
                scale.x=scale.x, plot.it=plot.it, br=br, ftype=ftype )
        }
        if ( ngroups > 2 )
        {
            cv <- cv.sgpls.multi( x=x, y=y, fold=fold, K=K, eta=eta,
                scale.x=scale.x, plot.it=plot.it, br=br, ftype=ftype )
        }
    } else
    {
        stop( "Check the response vector!" )
    }
    
    invisible( cv )
}
