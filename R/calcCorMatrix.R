

calcCorMatrix <- function(exprMatrix,corMethod,exprName){
	#
	# Args:
	#
	# Returns:
	#
	
	#=======================================
	#=======================================
	
	# cor calculates correlations between *columns* => take transpose of exprMatrix
	corMatrix <- cor(t(exprMatrix),method=corMethod)
	
	# Meta data for the for the corMatrix
	attr(corMatrix,"expr") <- exprName
	attr(corMatrix,"method") <- corMethod
	
	return(corMatrix)
	
}
