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
#' @param distance_matrix A matrix of the distances on the global network
#' @param gene_id The name of the gene to which GeneSurrounder is applied 
#' @param genes_assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network 
#' @param perm_geneStats_matrix A matrix of the resampled differential expression
#' @param sizes The number of genes assayed and on the network in each neighborhood
#' @importFrom pcaPP cor.fk
#'
Resample.DecayDE <- function(distance_matrix, gene_id, genes_assayedETnetwork,
                             diameter, perm_geneStats_matrix, sizes) {

  numResamples <- nrow(perm_geneStats_matrix)
  distances <- distance_matrix[gene_id, genes_assayedETnetwork]

  num_genes <- length(genes_assayedETnetwork) - 1
  geneid_d <- which(sizes == num_genes)[1]

  # null_tau b matrix of numResamples (default 1000) rows and
  # length(distance_matrix) columns b/c stacks columns
  null_tau_b <- vapply(1:geneid_d, function(RADIUS) {
    # Excludes gene j with distances > 0 
    igenes_distances <- distances[distances <= RADIUS & distances > 0]
    igenes_names <- names(igenes_distances)

    null <- vapply(1:numResamples, function(resample_index) {
      tau_b <- abs(perm_geneStats_matrix[resample_index, igenes_names])
      pcaPP::cor.fk(tau_b, igenes_distances)
    }, numeric(1))

    # X <- cbind(igenes_distances,
    # t( abs(geneStats$resampled[,igenes_names]) ) 
    # )

    # Y <- cor.fk(X)
  }, numeric(numResamples)) # numeric is length of one of the columns vapply stacks

  # 2:genid_d means number of rows in null_tau_b is not geneid_d
  null_tau_b <- cbind(null_tau_b, matrix(
    rep(null_tau_b[, geneid_d], diameter - geneid_d),
    nrow = numResamples)
  )

  # cf. ?vapply  
  # If FUN.VALUE is not an array, the result is a matrix with 
  # length(FUN.VALUE) rows and length(X) columns,
  # return matrix distance by resample number a la Resample.SI
  return(t(null_tau_b))
}


# 2: gene id. : cbind first row.
# preallocate matrix size instead of cbind?



