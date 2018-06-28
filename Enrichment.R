# Mehrnoosh Oghbaie
# 06/27/2018


library("GOstats")
library("org.Mm.eg.db")

Enrichment <- function(gene, ontology){
  universe = Lkeys(org.Mm.egGO)
  # Ontology (BP, CC, MF)
  param <- new("GOHyperGParams", geneIds = gene,
               universeGeneIds=universe, annotation="org.Mm.eg.db", 
               ontology=ontology,pvalueCutoff = 0.05)
  hyp <- hyperGTest(param)
  return(summary(hyp,categorySize=2000))
}


#Enrichment(as.character(UpInFemale2$gene.id),"BP")

