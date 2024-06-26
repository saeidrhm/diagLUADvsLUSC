##source: https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/enrichplot/inst/doc/enrichplot.html
library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de) ##edo@result

edo2 <- gseNCG(geneList) 

barplot(edo, showCategory=20)

library(enrichplot)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p1
#p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
#plot_grid(p1, p2, ncol=2)


## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)

upsetplot(edo)

heatplot(edox)

heatplot(edox, foldChange=geneList)

#
x2 <- pairwise_termsim(edo)
emapplot(x2)


#####################################
library(org.Hs.eg.db)
library(DOSE)
data(geneList)
GeneList <- read.delim("data/SelectedCorrelatedGenes_v1.txt", header=F)
GeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= unlist(GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")

edo <- enrichDGN(GeneListENTREZID$ENTREZID) ##edo@result

#barplot(edo)


library(enrichplot)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p1
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox)

upsetplot(edo)

heatplot(edox)


x2 <- pairwise_termsim(edo)
emapplot(x2)

##############################################
GeneList <- read.delim("data/SelectedCorrelatedGenes_v1.txt", header=F)
PathwayTable <- read.delim("data/PathwayTable_v2.txt", header=T)

library(org.Hs.eg.db)
GeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= unlist(GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
	
NumofPW <- 12
##easy way: embed in result data
newresult <- data.frame(matrix(nrow=NumofPW,ncol=0))
newresult$ID <- c("D001","D002","D003","D004","D005","D006","D007","D008","D009","D010","D011","D012")
newresult$Description <- c("IL17 Signaling in Psoriasis",
"IL-17 signaling pathway",
"IL23-mediated signaling events Homo sapiens",
"Interleukin-23-mediated signaling events",
"Jak-STAT signaling pathway",
"IL-6/JAK/STAT3 Signaling",
"TNF signaling pathway",
"TNF-alpha Signaling via NF-kB",
"MAPK and NFkB signaling pathways inhibited by Yersinia YopJ WP3849",
"Keratinization R-HSA-6805567",
"IL-2 Receptor Beta Chain in T cell Activation Homo sapiens h il2rbPathway",
"IL2-mediated signaling events Homo sapiens")
newresult$GeneRatio <- rep("1/10",NumofPW)
newresult$BgRatio <- rep("1/10",NumofPW)
newresult$pvalue <- c(
1.21e-06,
1.52e-04,
2.67e-04,
2.67e-04,
0.001595447,
6.49e-03,
0.002699649,
2.83e-02,
4.48e-03,
7.45e-05,
8.35e-03,
1.15e-02)
newresult$p.adjust <- c(
4.17e-04,
0.008884007,
1.22e-02,
0.028391784,
0.059848917,
1.39e-01,
0.05055706,
0.209177052,
0.144570615,
0.009831843,
0.148106686,
0.122652397)
newresult$qvalue <- newresult$p.adjust#rep(0.01,NumofPW)

#newresult <- newresult[order(newresult$p.adjust),]
GetGeneCount <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(nrow(CurrGeneListENTREZID))
}
GetGeneIds <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(paste(c(CurrGeneListENTREZID$ENTREZID),collapse="/"))
}
GetGeneSymb <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(paste(c(CurrGeneListENTREZID$ALIAS),collapse="/"))
}

newresult$geneID <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneIds(x)))
newresult$Count <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneCount(x)))


edo@result <- newresult 
edo@pvalueCutoff <- 0.5
edo@gene <- c(GeneListENTREZID$ENTREZID)


barplot(edo, showCategory=20)


library(enrichplot)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p1
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
#cnetplot(edox)

#upsetplot(edo)

#heatplot(edox)


#x2 <- pairwise_termsim(edo)
#emapplot(x2)





######################################################
##read delim pathwaygenes
PathwayTable <- read.delim("data/PathwayTable_v2.txt", header=T)
GeneList <- read.delim("data/SelectedCorrelatedGenes_v1.txt", header=F)

PathwayTableAddrList <- c("data/KEGG_2021_Human_v2.txt",
"data/GO_Biological_Process_2023_v2.txt",
"data/GO_Molecular_Function_2023_v2.txt",
"data/GO_Cellular_Component_2023_v2.txt")


PathwayTableAddrList <- c("data/KEGG_2021_Human_v2.txt")
PathwayTableAddrList <- c("data/GO_Biological_Process_2023_v2.txt")
PathwayTableAddrList <- c("data/GO_Molecular_Function_2023_v2.txt")
PathwayTableAddrList <- c("data/GO_Cellular_Component_2023_v2.txt")

CurrPathwayDB <- "KEGG"#"GOBP" "GOMF" "GOCC"

PathwayGeneTableList <- lapply(PathwayTableAddrList, function(x) return(read.delim(x, header=F)))
PathwayGeneTableDF <- do.call(rbind, PathwayGeneTableList)
rownames_PathwayGeneTableDF <- PathwayGeneTableDF[,1]
PathwayGeneTableDF <- as.data.frame(PathwayGeneTableDF[,-1])
rownames(PathwayGeneTableDF) <- rownames_PathwayGeneTableDF

library(org.Hs.eg.db)
GeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= unlist(GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")


CurrPathwayTable <- PathwayTable[PathwayTable$DB==CurrPathwayDB,]

CurrPathwayGeneTableDF <- PathwayGeneTableDF[CurrPathwayTable$Term,]
NumofPW <- nrow(CurrPathwayTable)
##easy way: embed in result data
newresult <- data.frame(matrix(nrow=NumofPW,ncol=0))
newresult$ID <- sapply(seq(1,NumofPW,1),function(x)paste0("D0", x,collapse=""))##c("D001","D002","D003","D004","D005","D006")
newresult$Description <- CurrPathwayTable$Term
newresult$GeneRatio <- rep("1/100",NumofPW)
newresult$BgRatio <- rep("1/100",NumofPW)
newresult$pvalue <- CurrPathwayTable$P.value
newresult$p.adjust <- CurrPathwayTable$Adjusted.P.value
newresult$qvalue <- newresult$p.adjust#rep(0.01,NumofPW)

#newresult <- newresult[order(newresult$p.adjust),]
GetGeneCount <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(nrow(CurrGeneListENTREZID))
}
GetGeneIds <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(paste(c(CurrGeneListENTREZID$ENTREZID),collapse="/"))
}
GetGeneSymb <- function(x){
CurrGeneList <- strsplit(PathwayTable$Genes[x],";")
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (CurrGeneList[[1]]),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(paste(c(CurrGeneListENTREZID$ALIAS),collapse="/"))
}

newresult$geneID <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneIds(x)))
newresult$Count <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneCount(x)))


GetGeneIdsPW <- function(GeneList){
library(org.Hs.eg.db)
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(paste(c(CurrGeneListENTREZID$ENTREZID),collapse="/"))
}
GetGeneIdsSepPW <- function(GeneList){
library(org.Hs.eg.db)
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(c(CurrGeneListENTREZID$ENTREZID))
}
GetGeneCountPW <- function(GeneList){
library(org.Hs.eg.db)
CurrGeneListENTREZID <- AnnotationDbi::select(org.Hs.eg.db, keys= (GeneList),   
    columns=c("ENTREZID","GENENAME"), keytype="ALIAS")
return(nrow(CurrGeneListENTREZID))
}
##convert genes to entrez
##find number of genes in each pathway
PWCount <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneCountPW(unlist(strsplit(CurrPathwayGeneTableDF[x]," "))
)))
PWgeneID <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneIdsPW(unlist(strsplit(CurrPathwayGeneTableDF[x]," "))
)))
PWgeneIDSep <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneIdsSepPW(unlist(strsplit(CurrPathwayGeneTableDF[x]," "))
)))


geneSets <- list()
lapply(seq(1, length(newresult$ID), 1), function(x) {geneSets[[newresult$ID[x]]] <<- PWgeneIDSep[[x]]})
#geneSets[["D001"]] <- PWgeneIDSep[[1]]



PWgeneIDUnique <- unique(unlist(PWgeneIDSep))


newresult$GeneRatio <- sapply(seq(1,NumofPW,1),function(x) return(paste(c(newresult$Count[x],'/',PWCount[x]),collapse="")))
newresult$BgRatio <- sapply(seq(1,NumofPW,1),function(x) return(paste(c(newresult$Count[x],'/',length(PWgeneIDUnique)),collapse="")))
rownames(newresult) <- newresult$ID 

#newresult$geneID <- PWgeneID#sapply(seq(1,NumofPW,1),function(x) return(GetGeneIds(x)))
#newresult$Count <- PWCount#sapply(seq(1,NumofPW,1),function(x) return(GetGeneCount(x)))

#newresult$GeneRatio <- sapply(strsplit(PathwayTable$Overlap,"/"),function(x)as.numeric(x[1])/as.numeric(x[2]))
#newresult$BGRatio <- sapply(strsplit(PathwayTable$Overlap,"/"),function(x)as.numeric(x[1])/as.numeric(x[2]))

edo@result <- newresult 
edo@pvalueCutoff <- 0.05
edo@gene <- c(GeneListENTREZID$ENTREZID)
edo@geneSets <- geneSets
edo@universe <- as.character(PWgeneIDUnique)




#Pso5SelectedPathways_barplot_v10
barplot(edo, showCategory=50)


library(enrichplot)
p1 <- dotplot(edo, showCategory=50) #+ ggtitle("dotplot for Lung Cancer")
p1


## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
#edox@geneID <- sapply(seq(1,NumofPW,1),function(x) return(GetGeneSymb(x)))


cnetplot(edox)

upsetplot(edo)

heatplot(edox)


x2 <- pairwise_termsim(edo)
emapplot(x2)