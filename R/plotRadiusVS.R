#' Plot results of GeneSurrounder 
#'
#' \code{plotRadiusVS} plots the results of GeneSurrounder using the output of the geneNIDG function
#'
#' @param geneNIDG.gene A dataframe that is the output of the geneNIDG function
#' 
plotRadiusVS <- function(geneNIDG.gene) {
  geneNIDG.hsa4171 <- geneNIDG.gene
  par(mfcol=c(4,1), mar=c(1,4,2,1), oma=c(2,2,2,2), cex.axis = 1.5)
  plot(geneNIDG.hsa4171$radius, -log10(geneNIDG.hsa4171$p.Sphere),
    type = 'b',
    cex = 2,
    pch = 21,
    lwd = 2,
    bg = 'lightblue',
    xlab = '',
    ylab = '-log10(p.Sphere)',
    ylim = c(0,5),
    main = 'MCM2 (DNA replication factor) identified as "mechanistic driver" ',
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

  plot(geneNIDG.hsa4171$radius, -log10(geneNIDG.hsa4171$p.Decay),
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
  plot(geneNIDG.hsa4171$radius, -log10(geneNIDG.hsa4171$p.Fisher),
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
  plot(geneNIDG.hsa4171$radius, geneNIDG.hsa4171$size,
    type = 'b',
    cex = 2,
    pch = 21,
    bg = 'lightblue',
    xlab = 'Neighborhood Radius',
    ylab = 'Number of Assayed Genes',
    lwd = 2)

  abline(h=2708)

  # http://seananderson.ca/2013/10/21/panel-letters.html

  mtext("(D)", 
    side = 3, 
    adj = 0.05, 
    line = -2.25,
    cex = 1)
}



