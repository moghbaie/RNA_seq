# Mehrnoosh Oghbaie
# 06/21/2018
# Conversion

entrezID2symbol <- function(y){
  library(org.Mm.eg.db)
  eToSym <- select(org.Mm.eg.db,
                   keys = rownames(y),
                   keytype = "ENTREZID",
                   columns="SYMBOL")
  return(eToSym$SYMBOL)
}


