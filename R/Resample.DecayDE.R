# ===========================================================================
# file: Resample.DecayDE.R
# description:
# requires: 
# author: Sahil Shah <sahil.shah@u.northwestern.edu>
# ==========================================================================



#' Decay of Differential Expression step
#'
#' Decay of Differential Expression tests whether the magnitude of 
#' differential expression of other genes j in the neighborhood is inversely 
#' related to the distance d(i,j) of gene j from gene i
#'
#' @param distance.matrix A matrix of the distances on the global network
#' @param gene.id The name of the gene to which GeneSurrounder is applied 
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network 
#' @param perm.geneStats.matrix A matrix of the resampled differential expression
#' @param sizes The number of genes assayed and on the network in each neighborhood
#'


Resample.DecayDE <- function(distance.matrix,
							 gene.id,
							 genes.assayedETnetwork,
							 diameter,
							 perm.geneStats.matrix,
							 sizes){


	distances <- distance.matrix[gene.id,
							 genes.assayedETnetwork]





	num.genes <- length(genes.assayedETnetwork) - 1
	geneid.d <- which(sizes == num.genes)[1]


	# null.tau b matrix of 1000 row 34 columns b/c stacks columns
	null.tau_b <- vapply(1:geneid.d,function(RADIUS){


		# Excludes gene j with distances > 0 
		igenes.distances <- distances[distances <= RADIUS
									  & distances > 0]
		igenes.names <- names(igenes.distances)





		null <- vapply(1:nrow(perm.geneStats.matrix),function(resample.index){

			return ( cor.fk( abs(perm.geneStats.matrix[resample.index,igenes.names]) , 
							igenes.distances) 

				   )

		},
		numeric(1))

		# X <- cbind(igenes.distances,
		# 	  t( abs(geneStats$resampled[,igenes.names]) ) 
		# 	  )

		# Y <- cor.fk(X)


	},
	numeric(1000)) # numeric is length of one of the columns vapply stacks

	# 2:genid.d means number of rows in null.tau_b is not geneid.d
	
	null.tau_b <- cbind(null.tau_b,
					matrix(
						rep(null.tau_b[, geneid.d],diameter - geneid.d),
						nrow = 1000)

						)



	# cf. ?vapply  
	# If FUN.VALUE is not an array, the result is a matrix with 
	# length(FUN.VALUE) rows and length(X) columns,

	# return matrix distance by resample number a la Resample.SI
	return(t(null.tau_b))


}


# 2: gene id. : cbind first row.
# preallocate matrix size instead of cbind?



