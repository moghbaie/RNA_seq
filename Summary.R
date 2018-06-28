# Mehrnoosh Oghbaie
# 06/20/2018

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicAlignments)
library(GenomicFeatures)
library(rstudioapi)
library(readxl)
library(Rsamtools)
library(edgeR)
library(limma)
library(DESeq)
library(DESeq2)
library(pheatmap)

rm(list=ls())
# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("conversion.R")
source("Enrichment.R")


bam.files <- list.files(path = "./data", pattern = ".bam", full.names = TRUE)
bamfiles <- BamFileList(bam.files, yieldSize=2000000)
seqinfo(bamfiles[1])

#1.  Sorting and Indexing bam files
for (x in bam.files[-1]){
  y <- paste0("./data/","Sorted_",tools::file_path_sans_ext( basename(x)))
  sortBam(x,y)
  indexBam( paste0("./data/","Sorted_",basename(x)))
}

bam.files.sorted <- list.files(path = "./data", pattern = "Sorted_.*.bam$", full.names = TRUE)


#2. Counting reads using summarizedoverlaps 

geneExons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,by="gene")
se <- summarizeOverlaps(features=geneExons, 
                        reads=bam.files.sorted,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )


#3. Adding experiment design as metadata to the summarized object
csvfile <- file.path("./data","Experiment_Design.xlsx")
sampleTable <- read_excel(csvfile)
colnames(sampleTable) <- c("num","filename","Sex","WT","AS","HBV")
sampleTable$WT <- as.character(sampleTable$WT)
sampleTable$AS <- as.character(sampleTable$AS)
sampleTable$HBV <- as.character(sampleTable$HBV)
sampleTable <- sampleTable[order(sampleTable$filename),]
sampleTable <- rbind(sampleTable[-8,],sampleTable[8,])
colData(se) <- DataFrame(sampleTable[,-1])
# Saving se
#save(se, file="se.RData")
load("se.RData")
# Choose one color per variable
col.pal <- c("AS"="green","noAS"="orange")
col.pal[as.vector(sampleTable$AS)]

## Dimensions
nb.samples <- ncol(assay(se))
head(summary(assay(se)))


#4. Plotting histogram of count per gene 
par(mfrow=c(3,1), mar=c(4,4,4,4))

hist(as.matrix(assay(se)), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(assay(se) + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

# Boxplots of gene count distributions per sample
boxplot(log2(assay(se) + epsilon), col=col.pal[as.vector(sampleTable$Sex)], pch=".", 
        horizontal=TRUE, cex.axis=0.5,main= "Gene count distribution per sample",
        las=1, ylab="Samples", xlab="log2(Counts +1)")

#5. Drawing Density plots between sample with and without acute stresw
par(mfrow=c(1,1))
library(affy)
plotDensity(log2(assay(se) + epsilon), lty=1, col=color, lwd=2)
grid()
legend("topright", legend=names(col.pal), col=col.pal, lwd=2)

#6. Creating a DGEList for use with edgeR
genetable <- data.frame(gene.id=rownames(se))

y <- DGEList(counts=assay(se), 
             samples=colData(se), 
             genes=genetable)
names(y)

#7. filter genes below X cpms in all samples 
y <- calcNormFactors(y)
cutoff <- 1
drop <- which(apply(cpm(y), 1, max) < cutoff)
y <- y[-drop,] 
dim(y) # number of genes left
#8. Drawing MDS plot
par(mar=c(4,4,4,4))
plotMDS(y, col = c(6,6,10,10, 13, 13, 13,16,16,16))
# Add the legend : associating samples with filename
par(new=TRUE, mar=c(0,0,1,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("top", c("group1: WT","group2: WT+AS","group3: HBV","group4: HBV+AS"),
       text.col=c(16,13,10,6), col=1:4, horiz=TRUE)


#Compare AS *Sex  between noWT/HBV  
# 01. Forming DESeqDataSet object
dds0 <- DESeqDataSetFromMatrix(countData = assay(se)[,1:4], colData = sampleTable[1:4,], design = ~ Sex*AS)
print(dds0)

##02. Normalizing using the method for an object of class"CountDataSet" 
dds.norm <-  estimateSizeFactors(dds0)
sizeFactors(dds.norm)

##03. Compare the raw and normalized count
## Checking the normalization
par(mfrow=c(2,2),cex.lab=0.7)
col.pal[as.vector(sampleTable$AS)]

boxplot(log2(counts(dds0)+epsilon),  col=col.pal[as.vector(sampleTable$AS[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=col.pal[as.vector(sampleTable$AS[1:4])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds0)+epsilon),  col=col.pal[as.vector(sampleTable$AS[1:4])], 
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=col.pal[as.vector(sampleTable$AS[1:4])], 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 

## Computing mean and variance
norm.counts <- counts(dds.norm, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
## Mean and variance relationship
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
## 04. Plotting mean to variance relationship
par(mfrow=c(1,1))
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")


## 05. estimation of dispersion parameter and plotting it
dds.disp <- estimateDispersions(dds.norm)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
plotDispEsts(dds.disp)
## 06. Performing differential expression call using wald test
alpha <- 0.0001
wald.test <- nbinomWaldTest(dds.disp)
res.DESeq2 <- results(wald.test, alpha=alpha, pAdjustMethod="none")
res.DESeq2 <- res.DESeq2[!is.na(res.DESeq2$pvalue),]
## What is the object returned by nbinomTest()
dim(res.DESeq2)
head(res.DESeq2)
res.DESeq2 <- res.DESeq2[order(res.DESeq2$padj),]
head(res.DESeq2)
##07. Draw an histogram of the p-values
hist(res.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")

##08. Volcano plot
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res.DESeq2$log2FoldChange, -log10(res.DESeq2$padj))
plot(res.DESeq2$log2FoldChange, -log10(res.DESeq2$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10( Adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")


gn.selected <- abs(res.DESeq2$log2FoldChange) > 1.5 & res.DESeq2$padj < alpha 
text(res.DESeq2$log2FoldChange[gn.selected],
     -log10(res.DESeq2$padj)[gn.selected],
     lab=entrezID2symbol(res.DESeq2[gn.selected ,]), cex=0.7)

##09. Plotting expression levels of the most differentially expressed gene
par(mfrow=c(2,2), mar=c(4,4,4,4))
barplot(counts(dds.norm, normalized=T)[rownames( res.DESeq2)[1],], names.arg=paste0(sampleTable[1:4,"Sex"]$Sex,"_",sampleTable[1:4,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[1:4])], main=entrezID2symbol( res.DESeq2[1,]), las=2, cex.names=0.7)
barplot(counts(dds.norm, normalized=T)[rownames( res.DESeq2)[2],], names.arg=paste0(sampleTable[1:4,"Sex"]$Sex,"_",sampleTable[1:4,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[1:4])], main=entrezID2symbol( res.DESeq2[2,]), las=2, cex.names=0.7)
barplot(counts(dds.norm, normalized=T)[rownames( res.DESeq2)[3],], names.arg=paste0(sampleTable[1:4,"Sex"]$Sex,"_",sampleTable[1:4,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[1:4])], main=entrezID2symbol(res.DESeq2[3,]), las=2, cex.names=0.7)
barplot(counts(dds.norm, normalized=T)[rownames( res.DESeq2)[4],], names.arg=paste0(sampleTable[1:4,"Sex"]$Sex,"_",sampleTable[1:4,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[1:4])], main=entrezID2symbol(res.DESeq2[4,]), las=2, cex.names=0.7)

# 10. Pheatmap Differential expression analysis based on the Negative Binomial
# apply a variance stabilizing transformation
dds0_esq <- DESeq(dds0)
assay(dds0_esq)
vsd0 <- vst(dds0_esq)
res0 <- results(dds0_esq)
summary(res0)
mat0 <- assay(vsd0)[ head(order(res0$padj),60),]
mat0 <- mat0 - rowMeans(mat0)
rownames(mat0) <- entrezID2symbol(mat0)
df <- as.data.frame(colData(vsd0)[,c("Sex","HBV","AS")])

pheatmap(mat0, annotation_col=df)
##11. Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
par(mfrow=c(1,1))
plotMA(res.DESeq2, colNonSig = "blue")
abline(h=c(-1:1), col="red")


#12. Performing Functional enrichment

library(gProfileR)

res.DESeq2.df <- na.omit(data.frame(res.DESeq2))
induced.sign <- res.DESeq2[abs(res.DESeq2$log2FoldChange) > 1.5 &  res.DESeq2$padj < alpha,]
# head(induced.sign)
# names(term.induced)

entrezID2symbol(induced.sign)
term.induced <- gprofiler(query=entrezID2symbol(induced.sign), organism="mmusculus")
term.induced <- term.induced[order(term.induced$p.value),]
# term.induced$p.value
knitr::kable(term.induced[1:10,c("term.name",
                          "term.size",
                          "query.size",
                          "overlap.size",
                          "recall",
                          "precision",
                          "p.value", 
                          "intersection")], 
      format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")


##Second part
#Compare AS *Sex  between WT/noHBV  
# 01. Forming DESeqDataSet object
dds1 <- DESeqDataSetFromMatrix(countData = assay(se)[,5:10], colData = sampleTable[5:10,], design = ~ Sex*AS)
print(dds1)

##02. Normalizing using the method for an object of class"CountDataSet" 
dds1.norm <-  estimateSizeFactors(dds1)
sizeFactors(dds1.norm)

##03. Compare the raw and normalized count
## Checking the normalization
par(mfrow=c(2,2),cex.lab=0.7)
col.pal[as.vector(sampleTable$AS)]

boxplot(log2(counts(dds1)+epsilon),  col=col.pal[as.vector(sampleTable$AS[5:10])], cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds1.norm, normalized=TRUE)+epsilon),  col=col.pal[as.vector(sampleTable$AS[5:10])], cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds1)+epsilon),  col=col.pal[as.vector(sampleTable$AS[5:10])], 
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds1.norm, normalized=TRUE)+epsilon), col=col.pal[as.vector(sampleTable$AS[5:10])], 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 

## Computing mean and variance
norm1.counts <- counts(dds1.norm, normalized=TRUE)
mean1.counts <- rowMeans(norm1.counts)
variance1.counts <- apply(norm1.counts, 1, var)
## Mean and variance relationship
mean1.var.col <- densCols(x=log2(mean1.counts), y=log2(variance1.counts))
## 04. Plotting mean to variance relationship
par(mfrow=c(1,1))
plot(x=log2(mean1.counts), y=log2(variance1.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")


## 05. estimation of dispersion parameter and plotting it
dds1.disp <- estimateDispersions(dds1.norm)

## A diagnostic plot which
## shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
plotDispEsts(dds1.disp)
## 06. Performing differential expression call using wald test
alpha <- 0.0001
wald.test1 <- nbinomWaldTest(dds1.disp)
res.DESeq3 <- results(wald.test1, alpha=alpha, pAdjustMethod="none")
res.DESeq3 <- res.DESeq3[!is.na(res.DESeq3$pvalue),]
## What is the object returned by nbinomTest()
dim(res.DESeq3)
head(res.DESeq3)
res.DESeq3 <- res.DESeq3[order(res.DESeq3$padj),]
head(res.DESeq3)
##07. Draw an histogram of the p-values
hist(res.DESeq3$padj, breaks=20, col="grey", main="DESeq3 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")

##08. Volcano plot
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(res.DESeq3$log2FoldChange, -log10(res.DESeq3$padj))
plot(res.DESeq3$log2FoldChange, -log10(res.DESeq3$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10( Adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")


gn.selected <- abs(res.DESeq3$log2FoldChange) > 1.5 & res.DESeq3$padj < alpha 
text(res.DESeq3$log2FoldChange[gn.selected],
     -log10(res.DESeq3$padj)[gn.selected],
     lab=entrezID2symbol(res.DESeq3[gn.selected ,]), cex=0.7)

##09. Plotting expression levels of the most differentially expressed gene
par(mfrow=c(3,3), mar=c(4,4,4,4))

barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[1],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol( res.DESeq3[1,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[2],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol( res.DESeq3[2,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[3],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol(res.DESeq3[3,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[4],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol(res.DESeq3[4,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[5],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol( res.DESeq3[5,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[6],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol( res.DESeq3[6,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[7],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol(res.DESeq3[7,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[8],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol(res.DESeq3[8,]), las=2, cex.names=0.7)
barplot(counts(dds1.norm, normalized=T)[rownames( res.DESeq3)[9],], names.arg=paste0(sampleTable[5:10,"Sex"]$Sex,"_",sampleTable[5:10,"AS"]$AS),
        col=col.pal[as.vector(sampleTable$AS[5:10])], main=entrezID2symbol(res.DESeq3[9,]), las=2, cex.names=0.7)



# 10. Pheatmap Differential expression analysis based on the Negative Binomial
# apply a variance stabilizing transformation
dds1_esq <- DESeq(dds1)
assay(dds1_esq)
vsd <- vst(dds1_esq)
res1 <- results(dds1_esq)
summary(res1)
mat <- assay(vsd)[ head(order(res1$padj),60),]
mat <- mat - rowMeans(mat)
rownames(mat) <- entrezID2symbol(mat)
df <- as.data.frame(colData(vsd)[,c("Sex","HBV","AS")])

pheatmap(mat, annotation_col=df)
##11. Draw a MA plot.
## Genes with adjusted p-values below 1% are shown
par(mfrow=c(1,1))
plotMA(res.DESeq3, colNonSig = "blue")
abline(h=c(-1:1), col="red")


#12. Performing Functional enrichment

library(gProfileR)

res.DESeq3.df <- na.omit(data.frame(res.DESeq3))
induced.sign <- res.DESeq3[abs(res.DESeq3$log2FoldChange) > 1.5 &  res.DESeq3$padj < alpha,]
# head(induced.sign)
# names(term.induced)

entrezID2symbol(induced.sign)
term.induced <- gprofiler(query=entrezID2symbol(induced.sign), organism="mmusculus")
term.induced <- term.induced[order(term.induced$p.value),]
# term.induced$p.value
knitr::kable(term.induced[1:10,c("term.name",
                                 "term.size",
                                 "query.size",
                                 "overlap.size",
                                 "recall",
                                 "precision",
                                 "p.value", 
                                 "intersection")], 
             format.args=c(engeneer=TRUE, digits=3), caption="**Table: functional analysis wit gProfileR. ** ")



