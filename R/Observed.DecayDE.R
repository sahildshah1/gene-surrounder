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
#' @importFrom pcaPP cor.fk
#'
Observed.DecayDE <- function(distance_matrix, gene_id, genes_assayedETnetwork,
                             diameter, geneStats_observed) {
  if (class(genes_assayedETnetwork) %in% c("character", "factor")) {
    if (!all(genes_assayedETnetwork %in% colnames(distance_matrix))) {
      bad_genes <- genes_assayedETnetwork[!(genes_assayedETnetwork %in% colnames(distance_matrix))]
      bad_genes_msg <- paste(bad_genes, collapse = ", ")
      bad_genes_num <- length(bad_genes)
      stop(paste("There are ", bad_genes_num, " genes that are not in 'distance_matrix'.",
                  "Here is a list of those genes not included:", bad_genes_msg, sep = "\n"))
    }
  } else if (class(genes_assayedETnetwork) %in% c("numeric", "integer")) {
    genes_assayedETnetwork <- as.integer(genes_assayedETnetwork)
    if (any(genes_assayedETnetwork <= 0 | genes_assayedETnwork > ncol(distance_matrix))) {
      stop("At least one of the supplied column numbers was less than zero or greater than ",
           "the total number of columns in 'distance_matrix'.")
    }
  }

  if (gene_id %in% rownames(distance_matrix)) {
    stop("The supplied 'gene_id', '", gene_id, "', is not found in 'distance_matrix'.")
  }

  distances <- distance_matrix[gene_id, genes_assayedETnetwork]

  observed_tau_b <- vapply(1:diameter, function(RADIUS) {
    # Excludes gene j with distances > 0 
    igenes_distances <- distances[distances <= RADIUS & distances > 0]
    igenes_names <- names(igenes_distances)
    tau_b <- abs(geneStats_observed[igenes_names])
    pcaPP::cor.fk(tau_b, igenes_distances)
  }, numeric(1))
}
