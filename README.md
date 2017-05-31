
**Network-based identification of disease genes in expression data: the GeneSurrounder method**
<!-- 
Sahil Shah sahil.shah AT u.northwestern DOT edu

--- -->

The GeneSurrounder package implements the method we previously developed to
identify disease- associated genes from expression data and an independent
network model of cellular interactions. We developed GeneSurrounder to find the
genes with neighbors on the network that are differentially expressed (with the
magnitude of the differential expression decreasing with distance from the
putative disease gene) and have correlated expression with the putative disease
gene. Since the differential expression of the neighbors of a putative disease
gene does not depend on their association with that gene, our algorithm consists
of two tests that are run independently of each other. Their results are then
combined to determine if the putative disease gene is a central candidate
disease gene.

--- 

In order to illustrate our method, we apply our algorithm to one study of high-
vs-low grade ovarian cancer from the publicly available and curated collection
curatedOvarianData (GEO accession GSE14764). We have constructed the global
network model from KEGG pathways. To follow the example:


1. Clone this repo and change directories to gene-surrounder/vignettes/figs. 
2. Open vignette-main.pdf and a R interpreter
3. In the R interpreter, change directories to gene-surrounder/vignettes/figs
and follow the example in vignette-main.pdf





<!-- Presented on Wed, 6/18/16 at Braun Research Group Meeting

---


View the **[slides](https://github.com/sahildshah1/shiny-groupmtg/blob/master/figs/main.pdf)**

This presentation introduces Shiny, building a Shiny app, and customizing /managing 
larger Shiny apps.

---

View the **[skeleton app](https://github.com/sahildshah1/shiny-groupmtg/tree/master/skeleton-app)**

Example directory / file structure  -->