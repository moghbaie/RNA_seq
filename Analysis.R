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

# Sorting and Indexing bam files
for (x in bam.files[-1]){
  y <- paste0("./data/","Sorted_",tools::file_path_sans_ext( basename(x)))
  sortBam(x,y)
  indexBam( paste0("./data/","Sorted_",basename(x)))
}

bam.files.sorted <- list.files(path = "./data", pattern = "Sorted_.*.bam$", full.names = TRUE)

# Summarized Experiment

geneExons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,by="gene")
se <- summarizeOverlaps(features=geneExons, 
                        reads=bam.files.sorted,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

dim(se)
head(assay(se))
rowRanges(se)
length(rowRanges(se))
str(metadata(rowRanges(se)))
colData(se)
colnames(se)

# Adding meta data 
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
save(se, file="se.RData")

# Creating a DGEList for use with edgeR
genetable <- data.frame(gene.id=rownames(se))

y <- DGEList(counts=assay(se), 
             samples=colData(se), 
             genes=genetable)
names(y)
# calculate normalization factors
y <- calcNormFactors(y)

# filter genes below X cpms in all samples 
cutoff <- 1
drop <- which(apply(cpm(y), 1, max) < cutoff)
y <- y[-drop,] 
dim(y) # number of genes left

attach(y$samples)
# Making the group
group <- interaction(Sex, WT, AS, HBV)
group

# Quick MDS plot
colnames(y$counts) <- group

par(mar=c(4,4,4,4))
plotMDS(y, col = c(6,6,10,10, 13, 13, 13,16,16,16))
# Add the legend : associating samples with filename
par(new=TRUE, mar=c(0,0,1,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("top", c("group1: WT","group2: WT+AS","group3: HBV","group4: HBV+AS"),
       text.col=c(16,13,10,6), col=1:4, horiz=TRUE)

###################### Limma-voom analysis 
mm <- model.matrix(~0 + group) # specify model with no intercept for easier contrasts
d <- voom(y, mm, plot = T)
mtext(side = 3, line = 0.5, text = "1-Factor Model Without Intercept")
fit <- lmFit(d, mm)
head(coef(fit)[,!is.na(coef(fit))[1,]])

# Comparison between female and male in group 1 (Wild Type)
contr <- makeContrasts(groupM.WT.noAS.noHBV - groupF.WT.noAS.noHBV, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp1 <- topTable(tmp, sort.by = "P", n = Inf)
tmp1$SYMBOL <- entrezID2symbol(tmp1)
tmp1 <- tmp1[,c("gene.id","SYMBOL","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp1$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp1, file = "Wild_Type.txt", row.names = F, sep = "\t", quote = F)
head(tmp1, 10)

# Comparison between female and male in group 2 (Wild Type and Acute Stress)
# No significance between Female with/ without HBV
contr <- makeContrasts(groupM.WT.AS.noHBV - groupF.WT.AS.noHBV, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$SYMBOL <- entrezID2symbol(tmp2)
tmp2 <- tmp2[,c("gene.id","SYMBOL","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp2, file = "Wild_Type_Stress_noHBV.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 10)

# Comparison between female and male in group 3 (HBV)
contr <- makeContrasts(groupM.noWT.noAS.HBV - groupF.noWT.noAS.HBV, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp3 <- topTable(tmp, sort.by = "P", n = Inf)
tmp3$SYMBOL <- entrezID2symbol(tmp3)
tmp3 <- tmp3[,c("gene.id","SYMBOL","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp3$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp3, file = "HBV.txt", row.names = F, sep = "\t", quote = F)
head(tmp3, 10)

# Comparison between female and male in group 34 (HBV+AS)
contr <- makeContrasts(groupM.noWT.AS.HBV - groupF.noWT.AS.noHBV, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp4 <- topTable(tmp, sort.by = "P", n = Inf)
tmp4$SYMBOL <- entrezID2symbol(tmp4)
tmp4 <- tmp4[,c("gene.id","SYMBOL","logFC","AveExpr","P.Value","adj.P.Val")]
UpInAct4 <- tmp4$adj.P.Val < 0.05 & tmp4$logFC <0
UpInAct4 <- as.integer(UpInAct4)
names(UpInAct4) <- tmp4$SYMBOL
# Table of upregulated
table(UpInAct4)
length(which(tmp4$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp4, file = "HBV_Acute_Stress.txt", row.names = F, sep = "\t", quote = F)
head(tmp4, 10)



# Building RangedSummarizedExperiment
sampleTable$group <- group
dds <- DESeqDataSetFromMatrix(countData = assay(se),
                                 colData = sampleTable,
                                 design = ~Sex*HBV*AS)
dds <- DESeq(dds)
assay(dds)
vsd <- vst(dds)
plotPCA(vsd, "Sex")
res <- results(dds)
summary(res)
mat <- assay(vsd)[ head(order(res$padj),60),]
mat <- mat - rowMeans(mat)
rownames(mat) <- entrezID2symbol(mat)
df <- as.data.frame(colData(vsd)[,c("Sex","HBV","AS")])
par(mar=c(4,4,4,4))
pheatmap(mat[,1:4], annotation_col=df)


dds2 <- DESeqDataSetFromMatrix(countData = assay(se),
                              colData = sampleTable,
                              design = ~Sex*WT*AS)
dds2 <- DESeq(dds2)
vsd2 <- vst(dds2)
plotPCA(vsd2, "Sex")
res2 <- results(dds2)
summary(res2)
mat2 <- assay(vsd2)[ head(order(res2$padj),60),]
mat2 <- mat2 - rowMeans(mat2)
rownames(mat2) <- entrezID2symbol(mat2)
df2 <- as.data.frame(colData(vsd2)[,c("Sex","WT","AS")])
par(mar=c(4,4,4,4))
pheatmap(mat2[,5:10], annotation_col=df2)


# Enrichment Analysis with top 20 most expressed genes (Upregulated & Downregulated in Female)in each group
# Ontology (BP, CC, MF)
UpInFemale1 <- head(tmp1[tmp1$logFC<0,],30)
DownInFemale1 <- head(tmp1[tmp1$logFC>0,],30)
Up1BP <- Enrichment(as.character(UpInFemale1$gene.id),"BP")
Down1BP <- Enrichment(as.character(DownInFemale1$gene.id),"BP")
Up1CC <- Enrichment(as.character(UpInFemale1$gene.id),"CC")
Down1CC <- Enrichment(as.character(DownInFemale1$gene.id),"CC")
Up1MF <- Enrichment(as.character(UpInFemale1$gene.id),"MF")
Down1MF <- Enrichment(as.character(DownInFemale1$gene.id),"MF")

UpInFemale2 <- head(tmp2[tmp2$logFC<0,],30)
DownInFemale2 <- head(tmp2[tmp2$logFC>0,],30)
Up2BP <- Enrichment(as.character(UpInFemale2$gene.id),"BP")
Down2BP <- Enrichment(as.character(DownInFemale2$gene.id),"BP")
Up2CC <- Enrichment(as.character(UpInFemale2$gene.id),"CC")
Down2CC <- Enrichment(as.character(DownInFemale2$gene.id),"CC")
Up2MF <- Enrichment(as.character(UpInFemale2$gene.id),"MF")
Down2MF <- Enrichment(as.character(DownInFemale2$gene.id),"MF")

UpInFemale3 <- head(tmp3[tmp3$logFC<0,],30)
DownInFemale3 <- head(tmp3[tmp3$logFC>0,],30)
Up3BP <- Enrichment(as.character(UpInFemale3$gene.id),"BP")
Down3BP <- Enrichment(as.character(DownInFemale3$gene.id),"BP")
Up3CC <- Enrichment(as.character(UpInFemale3$gene.id),"CC")
Down3CC <- Enrichment(as.character(DownInFemale3$gene.id),"CC")
Up3MF <- Enrichment(as.character(UpInFemale3$gene.id),"MF")
Down3MF <- Enrichment(as.character(DownInFemale3$gene.id),"MF")

UpInFemale4 <- head(tmp4[tmp4$logFC<0,],30)
DownInFemale4 <- head(tmp4[tmp4$logFC>0,],30)
Up4BP <- Enrichment(as.character(UpInFemale4$gene.id),"BP")
Down4BP <- Enrichment(as.character(DownInFemale4$gene.id),"BP")
Up4CC <- Enrichment(as.character(UpInFemale4$gene.id),"CC")
Down4CC <- Enrichment(as.character(DownInFemale4$gene.id),"CC")
Up4MF <- Enrichment(as.character(UpInFemale4$gene.id),"MF")
Down4MF <- Enrichment(as.character(DownInFemale4$gene.id),"MF")

rowBP <- unique(c(Up1BP$Term, Up2BP$Term,Up3BP$Term, Up4BP$Term, Down1BP$Term, Down2BP$Term, Down3BP$Term, Down4BP$Term))

SummaryBP <- matrix(0, nrow =length(rowBP), ncol=8)
colnames(SummaryBP) <- c("Up1BP", "Up2BP","Up3BP", "Up4BP", "Down1BP", "Down2BP", "Down3BP", "Down4BP")
rownames(SummaryBP) <- rowBP

SummaryBP[,"Up1BP"] <- Up1BP$Count[match(rownames(SummaryBP),Up1BP$Term)]
SummaryBP[,"Up2BP"] <- Up2BP$Count[match(rownames(SummaryBP),Up2BP$Term)]
SummaryBP[,"Up3BP"] <- Up3BP$Count[match(rownames(SummaryBP),Up3BP$Term)]
SummaryBP[,"Up4BP"] <- Up4BP$Count[match(rownames(SummaryBP),Up4BP$Term)]
SummaryBP[,"Down1BP"] <- Down1BP$Count[match(rownames(SummaryBP),Down1BP$Term)]
SummaryBP[,"Down2BP"] <- Down2BP$Count[match(rownames(SummaryBP),Down2BP$Term)]
SummaryBP[,"Down3BP"] <- Down3BP$Count[match(rownames(SummaryBP),Down3BP$Term)]
SummaryBP[,"Down4BP"] <- Down4BP$Count[match(rownames(SummaryBP),Down4BP$Term)]

rowBP2 <- paste(unique(c(Up1BP$GOBPID, Up2BP$GOBPID,Up3BP$GOBPID, Up4BP$GOBPID, Down1BP$GOBPID, Down2BP$GOBPID, Down3BP$GOBPID, Down4BP$GOBPID)),":",
               unique(c(Up1BP$Term, Up2BP$Term,Up3BP$Term, Up4BP$Term, Down1BP$Term, Down2BP$Term, Down3BP$Term, Down4BP$Term)))
rownames(SummaryBP) <- rowBP2
SummaryBP[is.na(SummaryBP)]<- 0
pheatmap(SummaryBP, fontsize=9)
         
rowCC <- unique(c(Up1CC$Term, Up2CC$Term,Up3CC$Term, Up4CC$Term, Down1CC$Term, Down2CC$Term, Down3CC$Term, Down4CC$Term))
SummaryCC <- matrix(0, nrow =length(rowCC), ncol=8)
colnames(SummaryCC) <- c("Up1CC", "Up2CC","Up3CC", "Up4CC", "Down1CC", "Down2CC", "Down3CC", "Down4CC")
rownames(SummaryCC) <- rowCC

SummaryCC[,"Up1CC"] <- Up1CC$Count[match(rownames(SummaryCC),Up1CC$Term)]
SummaryCC[,"Up2CC"] <- Up2CC$Count[match(rownames(SummaryCC),Up2CC$Term)]
SummaryCC[,"Up3CC"] <- Up3CC$Count[match(rownames(SummaryCC),Up3CC$Term)]
SummaryCC[,"Up4CC"] <- Up4CC$Count[match(rownames(SummaryCC),Up4CC$Term)]
SummaryCC[,"Down1CC"] <- Down1CC$Count[match(rownames(SummaryCC),Down1CC$Term)]
SummaryCC[,"Down2CC"] <- Down2CC$Count[match(rownames(SummaryCC),Down2CC$Term)]
SummaryCC[,"Down3CC"] <- Down3CC$Count[match(rownames(SummaryCC),Down3CC$Term)]
SummaryCC[,"Down4CC"] <- Down4CC$Count[match(rownames(SummaryCC),Down4CC$Term)]

rowCC2 <- paste(unique(c(Up1CC$GOCCID, Up2CC$GOCCID,Up3CC$GOCCID, Up4CC$GOCCID, Down1CC$GOCCID, Down2CC$GOCCID, Down3CC$GOCCID, Down4CC$GOCCID)),":",
                unique(c(Up1CC$Term, Up2CC$Term,Up3CC$Term, Up4CC$Term, Down1CC$Term, Down2CC$Term, Down3CC$Term, Down4CC$Term)))
rownames(SummaryCC) <- rowCC2
SummaryCC[is.na(SummaryCC)]<- 0
pheatmap(SummaryCC)

rowMF <- unique(c(Up1MF$Term, Up2MF$Term,Up3MF$Term, Up4MF$Term, Down1MF$Term, Down2MF$Term, Down3MF$Term, Down4MF$Term))
SummaryMF <- matrix(0, nrow =length(rowMF), ncol=8)
colnames(SummaryMF) <- c("Up1MF", "Up2MF","Up3MF", "Up4MF", "Down1MF", "Down2MF", "Down3MF", "Down4MF")
rownames(SummaryMF) <- rowMF

SummaryMF[,"Up1MF"] <- Up1MF$Count[match(rownames(SummaryMF),Up1MF$Term)]
SummaryMF[,"Up2MF"] <- Up2MF$Count[match(rownames(SummaryMF),Up2MF$Term)]
SummaryMF[,"Up3MF"] <- Up3MF$Count[match(rownames(SummaryMF),Up3MF$Term)]
SummaryMF[,"Up4MF"] <- Up4MF$Count[match(rownames(SummaryMF),Up4MF$Term)]
SummaryMF[,"Down1MF"] <- Down1MF$Count[match(rownames(SummaryMF),Down1MF$Term)]
SummaryMF[,"Down2MF"] <- Down2MF$Count[match(rownames(SummaryMF),Down2MF$Term)]
SummaryMF[,"Down3MF"] <- Down3MF$Count[match(rownames(SummaryMF),Down3MF$Term)]
SummaryMF[,"Down4MF"] <- Down4MF$Count[match(rownames(SummaryMF),Down4MF$Term)]

rowMF2 <- paste(unique(c(Up1MF$GOMFID, Up2MF$GOMFID,Up3MF$GOMFID, Up4MF$GOMFID, Down1MF$GOMFID, Down2MF$GOMFID, Down3MF$GOMFID, Down4MF$GOMFID)),":",
                unique(c(Up1MF$Term, Up2MF$Term,Up3MF$Term, Up4MF$Term, Down1MF$Term, Down2MF$Term, Down3MF$Term, Down4MF$Term)))
rownames(SummaryMF) <- rowMF2
SummaryMF[is.na(SummaryMF)]<- 0
pheatmap(SummaryMF)
