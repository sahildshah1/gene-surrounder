#===========================================#
#===========================================#
#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the observed and resampled 
#' differential expression of the genes should be calculated 
#'
#' @param expr A matrix (genes by samples) of expression values.
#' @param classLabels A factor with levels corresponding to class labels.
#' @param numResamples defaults to 1000. The number of resamples when calculating resampled differential expression. 
#' @importFrom limma treat lmFit
#' @export
# Calc gene level statistics & a null set of gene level stats (shuffle phenotype labels)
calcGeneTStats <- function(expr, classLabels, numResamples = 1000){
  # Calc gene level statistics
  # This code is being written using CurOvGradeKEGGnets[[2]]
  # 
  # Args:
  # expr: is a matrix of genes by samples
  # classLabels: is a vector of class labels (e.g. high vs low)
  # numrResamples: number of times the phenotype labels are shuffled
  # Returns:
  # observedStats: a vector of observed moderated t-statistics
  # permStats: a matrix of resampled moderated t statistics (resamplings are rows)

  # cf. d715_timecourse_contrasts, network_review_GSEAhat ? 
  # Should I save the fit so I have the gene p values etc...? 

  desMat <- model.matrix(~factor(classLabels))
  # treat is a limma function
  # Given a microarray linear model fit, compute moderated t-statistic, etc
  fit <- limma::treat(limma::lmFit(expr, desMat))
  observedStats <- fit$t[, 2]


  permStats <- sapply(1:numResamples, function(resampleLoopIndex) {
    # Shuffle the phenotype labels
    permLabels <- sample(classLabels, replace = FALSE)

    #Refit and recalculcate gene level statistics using permLabels
    permDesMat <- model.matrix(~factor(permLabels))
    permFit <- limma::treat(limma::lmFit(expr, permDesMat))
    #print(head(permFit$t[, 2]))
    return(permFit$t[, 2])
  })

  # Transpose permStats so the rows are resamplings
  permStats <- t(permStats)

  # List and return
  geneTStats <- list(observed = observedStats, resampled = permStats)

  return(geneTStats)	
}
