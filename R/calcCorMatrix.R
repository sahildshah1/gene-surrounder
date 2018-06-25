#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the correlation
#' between the expression of the genes should be calculated
#'
#' @param exprMatrix A matrix (genes by samples) of expression values.
#' @param corMethod A string with the correlation method to use.
#' @param exprName A string with the name of the expression matrix.
#' @param useMethod A string to specify which values to use.
#'   See "use" option under ?cor
#'
calcCorMatrix <- function(exprMatrix, corMethod, exprName, useMethod){
  #
  # Args:
  #
  # Returns:
  #

  #=======================================
  #=======================================

  # cor calculates correlations between *columns* => take transpose of exprMatrix
  corMatrix <- cor(t(exprMatrix), method = corMethod, use = useMethod)

  # Meta data for the for the corMatrix
  attr(corMatrix, "expr") <- exprName
  attr(corMatrix, "method") <- corMethod

  corMatrix	
}
