


#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the correlation
#' between the expression of the genes should be calculated
#'
#' @param exprMatrix A matrix (genes by samples) of expression values.
#' @param corMethod A string contating the correlation method to use.
#' @param exprName A string containg the name of the expression matrix.
#'


calcCorMatrix <- function(exprMatrix,corMethod,exprName,useMethod){
	#
	# Args:
	#
	# Returns:
	#
	
	#=======================================
	#=======================================
	
	# cor calculates correlations between *columns* => take transpose of exprMatrix
	corMatrix <- cor(t(exprMatrix),method=corMethod,use=useMethod)
	
	# Meta data for the for the corMatrix
	attr(corMatrix,"expr") <- exprName
	attr(corMatrix,"method") <- corMethod
	
	return(corMatrix)
	
}
