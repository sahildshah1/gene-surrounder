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
#' @param genes_assayedETnetwork The names of the genes that are assayed and on the network
#'
Observed.SI <- function(gene_id, distance_matrix, cor_matrix,
                        diameter, genes_assayedETnetwork){
  # Vector of cor between j and all other genes on network (excluding j)
  cor_with_j <- cor_matrix[gene_id, setdiff(genes_assayedETnetwork, gene_id)]
  distances_to_j <- distance_matrix[gene_id, setdiff(genes_assayedETnetwork, gene_id)]
  observed_cor <- SumAbsCor(cor_with_j, diameter, distances_to_j)
}
