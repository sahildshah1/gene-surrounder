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
#' @param distance_matrix A matrix of the distances on the global network
#' @param gene_id The name of the gene to which GeneSurrounder is applied 
#' @param genes_assayedETnetwork The names of the genes that are assayed and on the network
#' @param diameter The diameter of the global network 
#' @param geneStats_observed A vector of the observed differential expression
#'
Observed.DecayDE <- function(distance_matrix, gene_id, genes_assayedETnetwork,
                             diameter, geneStats_observed) {

  distances <- distance_matrix[gene_id, genes_assayedETnetwork]

  observed_tau_b <- vapply(1:diameter, function(RADIUS) {
    # Excludes gene j with distances > 0 
    igenes_distances <- distances[distances <= RADIUS & distances > 0]
    igenes_names <- names(igenes_distances)
    return(cor.fk(abs(geneStats_observed[igenes_names]), igenes_distances))
  }, numeric(1))
}
