library(clusterProfiler)
library(tidyverse)
library(limma)
library(GEOquery)
library(tinyarray)
library(org.Hs.eg.db)
library(openxlsx)
library(BiocManager)
library(DESeq2)
library(edgeR)

data <- read.csv("batch.csv",header = T)
rownames(data) <- data[,1]
library(FactoMineR)
library(sva)
library(factoextra)
pre.pca <- PCA(t(expr),graph = FALSE)
fviz_pca_ind(pre.pca,geom= "point",col.ind = data$batch, addEllipses =
               TRUE,legend.title="Group")
model <- model.matrix(~as.factor(data$type))
combat_Expr <- ComBat(dat = expr,batch = data$batch,mod = model)
af1.pca <- PCA(t(combat_Expr),graph = FALSE)
combat.pca <- PCA(t(combat_Expr),graph = FALSE)
fviz_pca_ind(combat.pca,geom= "point",col.ind = data$batch, addEllipses = TRUE,
             legend.title="Group")
write.csv(combat_Expr,file = "消除样本间差异.csv",quote = F)


if(!file.exists("GSE88966_eSet.Rdata")){GEO_file <- getGEO('GSE88966', destdir =
                                                             'E:\\GSE88966', getGPL = T)
save(GEO_file, file = "GSE88966_eSet.Rdata")}
load("GSE88966_eSet.Rdata")
GEO_file[[1]]
exp <- exprs(GEO_file[[1]])
plate <- fData(GEO_file[[1]])
clinical <- pData(GEO_file[[1]])
write.csv(clinical,"GSE88966_clin.csv")
exp <- as.matrix(exp)
boxplot(exp)
ID <- data.frame(ID_REF = plate$ID,Gene_Symbol = plate$`Gene symbol` )

ID$Gene_Symbol <- data.frame(sapply(ID$Gene_Symbol, function(x)unlist
                                    (strsplit(x,"///"))[1]),stringsAsFactors = F)[,1]
exp <- as.data.frame(exp)
exp$ID_REF <- rownames(exp)
exp <- merge(exp, ID, by = 'ID_REF')
exp[,grep("Gene_Symbol",colnames(exp))]<-trimws(exp[,grep("Gene_Symb ol",
                                                          colnames(exp))])
exp[exp==""] <- NA
exp <- na.omit(exp)
exp <- as.data.frame(exp)
table(duplicated(exp$Gene_Symbol))
exp1 <- avereps(exp,ID = exp$Gene_Symbol)
exp1 <- as.data.frame(exp1)
rownames(exp1) <- exp1$Gene_Symbol
exp1 <- exp1[,-c(1,ncol(exp1))]
write.csv(exp1,"GSE88966_average expression.csv")

distribution = combat_Expr[,c(1)]

hist(distribution)
subc1nM <- read.csv("GSE88966-subc-2^-1nM DESeq.csv")
rownames(subc1nM) <- subc1nM$X
subc1nM <- subc1nM[,-1]
subc1nM = round(subc1nM)
conditions = data.frame(conditions=factor(c(rep("1nM",3),rep("ctrl",3))))
subject <- c("1","2","3","1","2","3")
conditions$subject <- subject
rownames(conditions) = colnames(subc1nM)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = subc1nM, colData =
                                              conditions, design = ~ subject + conditions)
dds = DESeq(ddsFullCountTable)
contrast = c("conditions","1nM","ctrl")
res = results(dds, contrast)
res$padj[is.na(res$padj)] <- 1
res = as.data.frame(res[order(res$pvalue),])
res$change <- ifelse(res$pvalue<0.05 & abs(res$log2FoldChange)>=0.6,
                     ifelse(res$log2FoldChange>0.6,"up","down"), "stable")
write.csv(res,file = "GSE88966-subc-1nM DESeq-DEG.csv",quote = F)

library(org.Hs.eg.db)
library(clusterProfiler)
markers <- read.csv("D:\\DEGs.csv")
gene = markers$SYMBOL
gene = bitr(gene,fromType="SYMBOL", toType= "ENTREZID", OrgDb =
              "org.Hs.eg.db")
gene = dplyr::distinct(gene,SYMBOL,.keep_all=T)
data_all <- GC_markers %>% inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% arrange(desc(avg_log2FC))
geneList = data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$ENTREZID
GSEA_GOBP <- gseGO(geneList,ont = "BP",org.Hs.eg.db,keyType =
                     "ENTREZID",pvalueCutoff = 1)
GSEA_GOBP<-setReadable(GSEA_GOBP,OrgDb="org.Hs.eg.db", keyType =
                         "ENTREZID")
GSEA_GOCC <- gseGO(geneList,ont = "CC",org.Hs.eg.db,keyType =
                     "ENTREZID",pvalueCutoff = 1)
GSEA_GOCC<-setReadable(GSEA_GOCC,OrgDb="org.Hs.eg.db",keyType =
                         "ENTREZID")
GSEA_GOMF <- gseGO(geneList,ont = "MF",org.Hs.eg.db,keyType =
                     "ENTREZID",pvalueCutoff = 1)
GSEA_GOMF<-setReadable(GSEA_GOMF,OrgDb="org.Hs.eg.db", keyType =
                         "ENTREZID")
GSEA_KEGG <- gseKEGG(geneList, organism = "hsa",keyType = "kegg",pvalueCutoff
                     = 1)
GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb="org.Hs.eg.db", keyType =
                           "ENTREZID")
GSEA_WP <- gseWP(geneList,organism = "Homo sapiens",pvalueCutoff = 1)
GSEA_WP <-setReadable(GSEA_WP, OrgDb="org.Hs.eg.db", keyType = "ENTREZID")
library(ReactomePA)
GSEA_RET <- gsePathway(geneList,organism = "human",pvalueCutoff = 1)
GSEA_RET <-setReadable(GSEA_RET, OrgDb="org.Hs.eg.db", keyType =
                         "ENTREZID")

