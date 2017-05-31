
# ===========================================================================
# file: geneNIDG.R
# description:
# requires: 
# author: Sahil Shah <sahil.shah@u.northwestern.edu>
# ==========================================================================

#' GeneSurrounder wrapper
#'
#' The GeneSurrounder method consists of two tests that are run 
#' independently of each other and then combined to determine if 
#' the putative disease gene is a “disruptive” candidate disease gene 
#' meeting both criteria.
#'
#' @param distance.matrix A matrix of the distances on the global network
#' @param cor.matrix A matrix of correlations between the expression of the genes
#' @param geneStats.observed A vector of the observed differential expression
#' @param perm.geneStats.matrix A matrix of the resampled differential expression
#' @param genes.assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network 
#' @param num.Sphere.resamples The number of resamples when running the Sphere of Influence Procedure
#' @param gene.id The name of the gene to which GeneSurrounder is applied 


geneNIDG <- function(distance.matrix,
					 cor.matrix,
					 geneStats.observed,
					 perm.geneStats.matrix,
					 genes.assayedETnetwork,
					 diameter,
					 num.Sphere.resamples,
					 gene.id){


	# size -------------------------------------------------------------------


	distances <- distance.matrix[gene.id,
							 genes.assayedETnetwork]


	sizes <- vapply(1:diameter,function(RADIUS){

		igenes.distances <- distances[distances <= RADIUS
									  & distances > 0]

		length(igenes.distances)


	},
	numeric(1))



	# gene.DecayDE --------------------------------------------------------------

	# vector of observed tau b and null tau b a la observed.cor and resampled.cor

	observed.tau_b <- Observed.DecayDE(distance.matrix,
								 		gene.id,
								 		genes.assayedETnetwork,
								 		diameter,
								 		geneStats.observed)


	null.tau_b <- Resample.DecayDE(distance.matrix,
								 gene.id,
								 genes.assayedETnetwork,
								 diameter,
								 perm.geneStats.matrix,
								 sizes)


	num.Decay.resamples <- nrow(perm.geneStats.matrix)

	# proportion of null taub \LEQ observed i.e more discordant 
	p.Decay <- vapply(1:diameter,function(distance){

		 	observed.p <- observed.tau_b[distance]
			
			null.p <- null.tau_b[distance,1:num.Decay.resamples]

			return(length(null.p[null.p <= observed.p])/length(null.p))

		 },
		 numeric(1))

	p.Decay[p.Decay == 0] <- 1/(num.Decay.resamples+1)	


	# gene.SI ----------------------------------------------------------------


	observed.cor <- Observed.SI(gene.id,
									distance.matrix,
									cor.matrix,
									diameter,
									genes.assayedETnetwork)


	resampled.cor <- Resample.SI( gene.id,
									  distance.matrix,
									  cor.matrix,
									  diameter,
									  num.Sphere.resamples,
									  genes.assayedETnetwork)

	# proportion of null sumabscor \GEQ 
	p.Sphere <- vapply(1:diameter,function(distance){

	 	#observed.p <- observed.cor[[distance]]$cor
	 	observed.p <- observed.cor[distance]
		
		#null.p <- resampled.cor[[distance]]$cor
		null.p <- resampled.cor[distance,1:num.Sphere.resamples ]

		return(length(null.p[null.p >= observed.p])/length(null.p))


	 },
	 numeric(1))

	p.Sphere[p.Sphere == 0] <- 1/(num.Sphere.resamples+1)	


	# no base R combine p values ? combine after running --------------------

	# ------------------------------------------------------------------------

	# ?adply
	# The most unambiguous behaviour is achieved when .fun returns a data frame

	return(data.frame(gene.id = rep(gene.id,diameter),
					  radius = 1:diameter,
					  size = sizes,
					  observed.tau_b = observed.tau_b,
					  p.Decay = p.Decay,
					  observed.cor = observed.cor,
					  p.Sphere = p.Sphere)

		)


}


