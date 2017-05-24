


#' Sphere of Influence Step 
#'
#' Sphere of Influence assesses if a 
#' candidate gene i meets the first criterion by testing if gene i 
#' is more strongly correlated with its network neighbors than with a 
#' random set of genes
#'
#' @param gene.id This is a description. 
#' @param distance.matrix This is a description. 
#' @param cor.matrix This is a description. 
#' @param diameter This is a description. 
#' @param genes.assayedETnetwork This is a description. 
#'
#'


Observed.SI <- function(gene.id,
						distance.matrix,
						cor.matrix,
						diameter,
						genes.assayedETnetwork){




	# Vector of cor between j and all other genes on network (excluding j)
	cor.with.j <- cor.matrix[gene.id,
							 setdiff(genes.assayedETnetwork,gene.id)]

	distances.to.j <- distance.matrix[gene.id,
									  setdiff(genes.assayedETnetwork,gene.id)]


	observed.cor <- SumAbsCor(cor.with.j,
			  				  diameter,
			  				  distances.to.j)


}