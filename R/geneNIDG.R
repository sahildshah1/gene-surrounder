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
#' @param distance_matrix A matrix of the distances on the global network
#' @param cor_matrix A matrix of correlations between the expression of the genes
#' @param geneStats_observed A vector of the observed differential expression
#' @param perm_geneStats_matrix A matrix of the resampled differential expression
#' @param genes_assayedETnetwork The names of the genes that are assayed and on the network.
#'  if 'all' is used, then all genes in the network are tested.
#' @param diameter The diameter of the global network 
#' @param num_Sphere_resamples The number of resamples when running the Sphere of Influence Procedure
#' @param gene_id The name of the gene to which GeneSurrounder is applied 
#' @export
geneNIDG <- function(distance_matrix, cor_matrix, geneStats_observed,
                     perm_geneStats_matrix, genes_assayedETnetwork,
                     diameter, num_Sphere_resamples, gene_id) {
  ### Sanity Checks ---------------------------------
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

  if (genes_assayedETnetwork == "all") {
    genes_assayedETnetwork <- colnames(distance_matrix)
  }

  # size -------------------------------------------------------------------
  distances <- distance_matrix[gene_id, genes_assayedETnetwork]


  sizes <- vapply(1:diameter, function(RADIUS) {
    igenes_distances <- distances[distances <= RADIUS & distances > 0]
    length(igenes_distances)
  }, numeric(1))

  # gene_DecayDE --------------------------------------------------------------
  # vector of observed tau b and null tau b a la observed_cor and resampled_cor
  observed_tau_b <- Observed.DecayDE(distance_matrix,
    gene_id,
    genes_assayedETnetwork,
    diameter,
    geneStats_observed
  )


  null_tau_b <- Resample.DecayDE(distance_matrix,
    gene_id,
    genes_assayedETnetwork,
    diameter,
    perm_geneStats_matrix,
    sizes
  )

  num_Decay_resamples <- nrow(perm_geneStats_matrix)

  # proportion of null taub \LEQ observed i.e more discordant 
  p_Decay <- vapply(1:diameter, function(distance) {
    observed_p <- observed_tau_b[distance]
    null_p <- null_tau_b[distance, 1:num_Decay_resamples]
    return(length(null_p[null_p <= observed_p])/length(null_p))
  }, numeric(1))

  p_Decay[which(p_Decay == 0)] <- 1/(num_Decay_resamples + 1)	


  # gene_SI ----------------------------------------------------------------
  observed_cor <- Observed.SI(gene_id,
    distance_matrix,
    cor_matrix,
    diameter,
    genes_assayedETnetwork
  )


  resampled_cor <- Resample.SI(gene_id,
    distance_matrix,
    cor_matrix,
    diameter,
    num_Sphere_resamples,
    genes_assayedETnetwork
  )

  # proportion of null sumabscor \GEQ 
  p_Sphere <- vapply(1:diameter, function(distance) {
    #observed_p <- observed_cor[[distance]]$cor
    observed_p <- observed_cor[distance]
    #null_p <- resampled_cor[[distance]]$cor
    null_p <- resampled_cor[distance, 1:num_Sphere_resamples]
    return(length(null_p[null_p >= observed_p])/length(null_p))
  }, numeric(1))

  p_Sphere[which(p_Sphere == 0)] <- 1 / (num_Sphere_resamples + 1)	

  # no base R combine p values ? combine after running --------------------
  # ------------------------------------------------------------------------
  # ?adply
  # The most unambiguous behaviour is achieved when .fun returns a data frame

  fisher_stat <- -2*(log(p_Decay) + log(p_Sphere))
  p_Fisher <- pchisq(fisher_stat, 2, lower.tail = FALSE)

  res <- data.frame(
    gene_id = rep(gene_id, diameter),
    radius = 1:diameter,
    size = sizes,
    observed_tau_b = observed_tau_b,
    p_Decay = p_Decay,
    observed_cor = observed_cor,
    p_Sphere = p_Sphere,
    p_Fisher = p_Fisher
  )
  res
}


