

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
#' @param num.resamples default 1000. 
#' @param genes.assayedETnetwork This is a description. 
#'
#'

 Resample.SI <- function(gene.id,
						distance.matrix,
						cor.matrix,
						diameter,
						num.Sphere.resamples=1000,
						genes.assayedETnetwork){

	# Vector of cor between j and all other genes on network (excluding j)
	cor.with.j <- cor.matrix[gene.id,
							 setdiff(genes.assayedETnetwork,gene.id)]

	distances.to.j <- distance.matrix[gene.id,
									  setdiff(genes.assayedETnetwork,gene.id)]


	x <- rep(0,length(cor.with.j))		

	A <- vapply(1:diameter, function(distance){

		x[distances.to.j <= distance] <- 1

		return(x)


	},
	numeric(length(x)) )

	A <- t(A)
								  


	# For each reasampling calculating sum of abs cor at each distance
	# Under this new implementaiotn, cor added up in d are also included in d + 1
	# and we only resample 1000 times instead of 34000 times.
	resampled.sum.abs.cor <- replicate(num.Sphere.resamples,


										# Simplify to [1:34] from [1:34,1]
										c( A %*% abs( sample(cor.with.j,
											   			 size = length(cor.with.j),
								                         replace = TRUE)
											    )
										  )

										# SumAbsCor(sample(cor.with.j,
										# 	   			 size = length(cor.with.j),
								  #                        replace = TRUE),
										# 		  diameter,
										# 		  distances.to.j
										# 		  )

									 )


	# apply(1:num.resamples, 1, function(resample.num){

	# 	# Sample with replacement form cor.with.j ----------------------------
	# 	resampled.cor <- sample(cor.with.j,
	# 							size = length(cor.with.j),
	# 							replace = TRUE)

	# 	# Calculate sum of resampled cor ----------------------------
	# 	sum.abs.cor <- apply(1:diameter, 1, function(distance){

	# 		#
	# 		sum( abs( resampled.cor[distances.to.j <= distance] ) )


	# 		})


	# 	return(sum.abs.cor)

		
	# })



}
