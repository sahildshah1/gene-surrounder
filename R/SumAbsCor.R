

#' Helper function for Sphere of Influence procedure 
#'
#' \code{SumAbsCor} returns the total observed correlation between gene i and its neighbors
#'
#' @param cor.vector A vector of correlations
#' @param diameter The diameter of the global network
#' @param distances.to.j A vector of distances to the gene to which GeneSurrounder is applied
#'
#' 

SumAbsCor <- function(cor.vector,diameter,distances.to.j){


	# Calculate sum of cor.vector ----------------------------
	sum.abs.cor <- vapply(1:diameter, function(distance){

		#
		sum( abs( cor.vector[distances.to.j <= distance] ), na.rm= TRUE)


	},
	numeric(1))

	return(sum.abs.cor)


}