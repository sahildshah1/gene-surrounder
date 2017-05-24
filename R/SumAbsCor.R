

#' Sum of vector elements.
#'
#' \code{SumAbsCor} returns the total observed correlation between gene i and its neighbors
#'
#' @param cor.vector This is a description. 
#' @param diameter This is a description.
#' @param distances.to.j This is a description.
#'
#' 

SumAbsCor <- function(cor.vector,diameter,distances.to.j){


	# Calculate sum of cor.vector ----------------------------
	sum.abs.cor <- vapply(1:diameter, function(distance){

		#
		sum( abs( cor.vector[distances.to.j <= distance] ) )


	},
	numeric(1))

	return(sum.abs.cor)


}