#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the distances on the 
#' global network should be calculated
#'
#' @param network An igraph network on which to calculate the distances.
#' @param networkName A string containing the name of the network. Required.
#' @param directionPaths defaults to "all". A string indicating how the
#'   distances should be calculated.
#' @param weightVector defaults to NULL. A numeric vector giving edge weights.
#' @importFrom igraph shortest.paths V
#' @export
calcAllPairsDistances <- function(network, networkName,
                                  directionPaths = "all",
                                  weightVector = NULL){
  #
  # Args:
  #
  # Returns:
  #
  # weightVector = null => use weight attribute if it exists otherwise don't
  shortestPathsMatrix <- igraph::shortest.paths(network,
    v = igraph::V(network),
    to = igraph::V(network),
    mode = directionPaths,
    weights = weightVector
  )

  attr(shortestPathsMatrix, "networkName") <- networkName
  return(shortestPathsMatrix)
}
