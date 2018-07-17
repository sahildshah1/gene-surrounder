#' Plot results of GeneSurrounder 
#'
#' \code{plotRadiusVS} plots the results of GeneSurrounder using the output of
#' the geneNIDG function.
#'
#' @param geneNIDG_gene A dataframe that is the output of the geneNIDG function
#' @param plotTitle the title at the top of the plot
#'
#' @details All of the panels of this plot have as the x-axis the
#' 'radius', defined as the distance away from the gene along the network graph.
#' This plot has four panels:
#'   \itemize{
#'     \item{(1)}{the correlation sphere p-value versus the radius}
#'     \item{(2)}{the differential expression decay p-value versus the radius}
#'     \item{(3)}{the Fisher's combined p-value versus the radius}
#'     \item{(4)}{the number of genes tested versus the radius}
#'   }
#'
#' The first three panels have a line to indicate where the p-values of
#' 0.01 and 0.05 are located, which are the typical thresholds to use with
#' this method. Note that the reported p-values are not corrected for multiple
#' hypotheses because the p-values are not independent.
#'
#' The last panel has a line to indicate the maximum number of genes tested
#' on the network.
#'
#' @return no return value; this function is used for its side effect
#' @importFrom graphics par plot abline mtext 
#' @export
plotRadiusVS <- function(geneNIDG_gene, plotTitle = "GeneSurrounder Results") {
  par(mfcol = c(4, 1), mar = c(1, 4, 2, 1),
      oma = c(2, 2, 2, 2), cex.axis = 1.5)
  plot(geneNIDG_gene$radius, -log10(geneNIDG_gene$p.Sphere),
    type = 'b',
    cex = 2,
    pch = 21,
    lwd = 2,
    bg = 'lightblue',
    xlab = '',
    ylab = '-log10(p.Sphere)',
    ylim = c(0,5),
    main = plotTitle,
    xaxt = 'n'
  )

  abline(h = -log10(0.01),lty=3,lwd=2)
  abline(h = -log10(0.05),lty=4,lwd=2)

  # http://seananderson.ca/2013/10/21/panel-letters.html

  mtext("(A)", 
    side = 3, 
    adj = 0.05, 
    line = -2.25,
    cex = 1)

  # axis(4)
  # mtext("(A)", side=4, line=3)


  # radius vs p.Decay -------------------------------------------------------

  plot(geneNIDG_gene$radius, -log10(geneNIDG_gene$p.Decay),
    type = 'b',
    cex = 2,
    pch = 21,
    lwd = 2,
    bg = 'lightblue',
    xlab = '',
    ylab = '-log10(p.Decay)',
    ylim = c(0,5),
    xaxt = 'n')

  abline(h = -log10(0.01),lty=3,lwd=2)
  abline(h = -log10(0.05),lty=4,lwd=2)

  # http://seananderson.ca/2013/10/21/panel-letters.html

  mtext("(B)", 
    side = 3, 
    adj = 0.05, 
    line = -2.25,
    cex = 1)

  # radius vs p.Fisher -------------------------------------------------------
  plot(geneNIDG_gene$radius, -log10(geneNIDG_gene$p.Fisher),
    type = 'b',
    cex = 2,
    pch = 21,
    lwd = 2,
    bg = 'lightblue',
    xlab = '',
    ylab = '-log10(p.NIDG)',
    # ylim = c(0,5),
    xaxt = 'n')

  abline(h = -log10(0.01),lty=3,lwd=2)
  abline(h = -log10(0.05),lty=4,lwd=2)

  # http://seananderson.ca/2013/10/21/panel-letters.html

  mtext("(C)", 
    side = 3, 
    adj = 0.05, 
    line = -2.25,
    cex = 1)

  # radius vs size      -------------------------------------------------------
  plot(geneNIDG_gene$radius, geneNIDG_gene$size,
    type = 'b',
    cex = 2,
    pch = 21,
    bg = 'lightblue',
    xlab = 'Neighborhood Radius',
    ylab = 'Number of Assayed Genes',
    lwd = 2)

  abline(h = max(geneNIDG_gene$size))

  # http://seananderson.ca/2013/10/21/panel-letters.html

  mtext("(D)", 
    side = 3, 
    adj = 0.05, 
    line = -2.25,
    cex = 1)
}
