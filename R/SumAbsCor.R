#' Helper function for Sphere of Influence procedure 
#'
#' \code{SumAbsCor} returns the total observed correlation between gene i and its neighbors
#'
#' @param cor_vector A vector of correlations
#' @param diameter The diameter of the global network
#' @param distances_to_j A vector of distances to the gene to which GeneSurrounder is applied
#'
#' 
SumAbsCor <- function(cor_vector,diameter,distances_to_j){
  # Calculate sum of cor_vector ----------------------------
  sum_abs_cor <- vapply(1:diameter, function(distance){
    cors <- cor_vector[distances_to_j <= distance]
    sum(abs(cors), na.rm= TRUE)
  }, numeric(1))

  sum_abs_cor
}
