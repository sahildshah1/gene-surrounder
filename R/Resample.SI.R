#' Sphere of Influence Step 
#'
#' Sphere of Influence assesses if a 
#' candidate gene i meets the first criterion by testing if gene i 
#' is more strongly correlated with its network neighbors than with a 
#' random set of genes
#'
#' @param gene_id The name of the gene to which GeneSurrounder is applied 
#' @param distance_matrix A matrix of the distances on the global network
#' @param cor_matrix A matrix of correlations between the expression of the genes
#' @param diameter The diameter of the global network 
#' @param num_Sphere_resamples defaults to 1000. The number of resamples when running the Sphere of Influence Procedure
#' @param genes_assayedETnetwork The names of the genes that are assayed and on the network
#'
Resample.SI <- function(gene_id, distance_matrix, cor_matrix, diameter,
                        num_Sphere_resamples = 1000, genes_assayedETnetwork) {
  # Vector of cor between j and all other genes on network (excluding j)
  cor_with_j <- cor_matrix[gene_id, setdiff(genes_assayedETnetwork,gene_id)]
  distances_to_j <- distance_matrix[gene_id, setdiff(genes_assayedETnetwork,gene_id)]

  x <- rep(0, length(cor_with_j))
  A <- vapply(1:diameter, function(distance){
    x[distances_to_j <= distance] <- 1
    x
  }, numeric(length(x)))

  A <- t(A)

  sampled_cors <- sample(cor_with_j, size = length(cor_with_j), replace = TRUE)								  
  # For each reasampling calculating sum of abs cor at each distance
  # Under this new implementaiotn, cor added up in d are also included in d + 1
  # and we only resample 1000 times instead of 34000 times.
  resampled_sum_abs_cor <- replicate(num_Sphere_resamples,
    # Simplify to [1:34] from [1:34,1]
    c(A %*% abs(sampled_cors))
  )

  resampled_sum_abs_cor										
  # SumAbsCor(sample(cor_with_j,
  # size = length(cor_with_j),
  # replace = TRUE),
  # diameter,
  # distances_to_j
  # )
  #)

  # apply(1:num_resamples, 1, function(resample_num){
  # 	# Sample with replacement form cor_with_j ----------------------------
  # 	resampled_cor <- sample(cor_with_j,
  # 							size = length(cor_with_j),
  # 							replace = TRUE)
  # 	# Calculate sum of resampled cor ----------------------------
  # 	sum_abs_cor <- apply(1:diameter, 1, function(distance){
  # 		#
  # 		sum( abs( resampled_cor[distances_to_j <= distance] ) )
  # 		})
  # 	return(sum_abs_cor)
  # })
}
