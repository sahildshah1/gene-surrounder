


#' Sphere of Influence Step 
#'
#' Sphere of Influence assesses if a 
#' candidate gene i meets the first criterion by testing if gene i 
#' is more strongly correlated with its network neighbors than with a 
#' random set of genes
#'
#' @param gene.id The name of the gene to which GeneSurrounder is applied  
#' @param distance.matrix A matrix of the distances on the global network
#' @param cor.matrix A matrix of correlations between the expression of the genes
#' @param diameter The diameter of the global network 
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
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