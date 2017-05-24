# ===========================================================================
# file: Observed.DecayDE.R
# description:
# requires: 
# author: Sahil Shah <sahil.shah@u.northwestern.edu>
# ==========================================================================

#' Decay of Differential Expression step
#'
#' Decay of Differential Expression, tests whether the magnitude of 
#' differential expression of other genes j in the neighborhood is inversely 
#' related to the distance d(i,j) of gene j from gene i
#'
#' @param distance.matrix This is a description. 
#' @param gene.id This is a description. 
#' @param genes.assayedETnetwork This is a description. 
#' @param diameter This is a description. 
#' @param geneStats This is a description. 
#'
#'

Observed.DecayDE <- function(distance.matrix,
							 gene.id,
							 genes.assayedETnetwork,
							 diameter,
							 geneStats.observed){



	distances <- distance.matrix[gene.id,
							 genes.assayedETnetwork]


	observed.tau_b <- vapply(1:diameter,function(RADIUS){


		# Excludes gene j with distances > 0 
		igenes.distances <- distances[distances <= RADIUS
									 & distances > 0]
		igenes.names <- names(igenes.distances)


		return(

			cor.fk(abs(geneStats.observed[igenes.names]), igenes.distances)

			)

	},
	numeric(1))

}
