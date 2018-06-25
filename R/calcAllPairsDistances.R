#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the distances on the 
#' global network should be calculated
#'
#' @param network An igraph network on which to calculate the distances.
#' @param directionPaths defaults to "all". A string indicating how the
#'   distances should be calculated.
#' @param weightVector defaults to NULL. A numeric vector giving edge weights.
#' @param networkName A string containg the name of the network.
#'
calcAllPairsDistances <- function(network, directionPaths = "all",
                                  weightVector = NULL, networkName){
  #
  # Args:
  #
  # Returns:
  #
  require(igraph)

  # weightVector = null => use weight attribute if it exists otherwise don't
  shortestPathsMatrix <- shortest.paths(network,
    v = V(network),
    to = V(network),
    mode = directionPaths,
    weights = weightVector)
	
  attr(shortestPathsMatrix,"networkName") <- networkName
  return(shortestPathsMatrix)
}
