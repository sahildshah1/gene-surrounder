
calcAllPairsDistances <- function(network,directionPaths="all",weightVector = NULL,networkName){
	#
	# Args:
	#
	# Returns:
	#
	
	require(igraph)


	# weightVector = null => use weight attribute if it exists otherwise don't
	shortestPathsMatrix <- shortest.paths(network, 
										 v=V(network), 
										 to=V(network),
										 mode = directionPaths,
										 weights = weightVector)
	
	 attr(shortestPathsMatrix,"networkName") <- networkName
	 return(shortestPathsMatrix)

}