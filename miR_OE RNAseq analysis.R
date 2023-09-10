setwd("")

#Packages
library("zoo")
library("tidyverse")
library(readxl)
library("DESeq2")
library("gplots")
library("MASS")
library("mitch")
library("limma")
library("kableExtra")
library("vioplot")
library(edgeR)
library(org.Hs.eg.db)
library(RColorBrewer)
library(FactoMineR)
library(pheatmap)
library(factoextra)
library(clusterProfiler)
library(enrichplot)
library(TissueEnrich)
library(multiMiR)
library(Biobase)
library(ggpubr)

#PRE-PROCESSING
###########
#Import gene names
ENSG <- read.table("Homo_sapiens.GRCh38.tx2gene.tsv", header = F, fill = T)
head(ENSG)

#Import in read counts
tmp <- read.table("3col.tsv",header=F)
#Merge df so that GeneName and geneID are included in dataframe.
colnames(tmp)
tmpALL<-merge(tmp, ENSG, by.x="V2", by.y="V1")
View(head(tmpALL))
tmpALL= tmpALL %>%
  rename("V2"="Transcript",
         "V1" = "ID",
         V3.x ="GeneCount",
         V2.y="GeneID",
         V3.y="GeneSYM")
head(tmpALL)

#Create Matrix and aggregate according to ENSG
colnames(tmpALL)
x <- as.matrix(acast(tmpALL, GeneSYM ~ ID, value.var= "GeneCount", fun.aggregate = sum))
x <- as.data.frame(x)
xx <- round(x)
xx[1:6,1:6]
dim(xx)
head(xx)
#Separate into EMSEMBL and SYMBOL into 2 columns
xx2<- rownames_to_column(xx)
xx3 <- tidyr::separate(xx2,rowname, into=c('ENSEMBL', 'SYMBOL'), sep = "\\_")
View(head(xx3))
#write.csv(xx3, "GeneCounts_ RAW DATA.csv")

#########
xx=read.csv("GeneCounts_ RAW DATA.csv", header = T)

#Remove version number (decimal point) from ENSG ID
xx$ENSEMBL = sub("\\..*", "", xx$ENSEMBL)
#xx2= xx[,-1]
colnames(xx)
xx3=xx%>%
  tidyr::unite(rowname, "ENSEMBL","SYMBOL") %>%
  tibble::column_to_rownames()

## QC analysis
par(mar=c(5,8,3,1))
barplot(colSums(xx3),horiz=TRUE,las=1,xlab="num reads")
sums <- colSums(xx3)
median(sums)
sums <- sums[order(sums)]
barplot(sums,horiz=TRUE,las=1,xlab="num reads",cex.names=0.8)
abline(v=20000000,col="red")

#Make Counts matrix for each miRNA
col= grep('SCRAM|30a', colnames(xx3))
m30a=xx3[,col]
col= grep('SCRAM|30c', colnames(xx3))
m30c=xx3[,col]

####### 
# START HERE
xx <- read_excel("GeneCounts_ RAW DATA.xlsx", 
                 sheet = 3, # Sheet 2 for 30a, Sheet 3 for 30c, sheet 4 for SCRAM only
                 trim_ws = T)

#Remove version number (decimal point) from ENSG ID
xx$ENSEMBL = sub("\\..*", "", xx$ENSEMBL)
colnames(xx)

xx3=xx%>%
  tidyr::unite(rowname, "ENSEMBL","SYMBOL") %>%
  tibble::column_to_rownames(var = "rowname")
head(xx3)
## QC analysis
par(mar=c(5,8,3,1))
barplot(colSums(xx3),horiz=TRUE,las=1,xlab="num reads")
sums <- colSums(xx3)
median(sums)
sums <- sums[order(sums)]
barplot(sums,horiz=TRUE,las=1,xlab="num reads",cex.names=0.8)
abline(v=20000000,col="red")

#Import in sample sheet
ss <- read_excel("PhenoTable.xlsx", 
                 sheet = 4, # Sheet 2 for 30a, Sheet 3 for 30c, sheet 4 for SCRAM only
                 trim_ws = T)
ss$Sex <- factor(ss$Sex, levels= c("MALE", "FEMALE"))
ss$miR <- factor(ss$miR, levels= c("CON", "30c")) #Don't forget to change this for 30a or 30c
ss$Sex.miR <- factor(ss$Sex.miR, levels = c("MALECON", "MALE30c", "FEMALECON", "FEMALE30c"))

#Make sure matrix (xx) is in same order as sample sheet (ss)
ss30a = ss %>% 
  tibble::column_to_rownames(var="CODE")
#write.csv(ss30a, "Pheno30a.csv")
m30a.2 <- xx3[,match(rownames(ss30a), colnames(xx3))]
#write.csv(m30a.2, "M30a.csv")
#Check they match
all(rownames(ss30a) == colnames(m30a.2))

# Logging data allows you to increase the distance between small measurements and decreases the distance between large measurements 
#log() is natural log (base 10)
#Log transform
#Create ExpressionSet Object

eset=ExpressionSet(assayData = as.matrix(m30a.2),
                   phenoData = AnnotatedDataFrame(ss30a))
exprs(eset) <- log(exprs(eset))
plotDensities(eset, legend =F)

dds <- DESeqDataSetFromMatrix(countData = m30a.2, 
                              colData = ss30a, 
                              design = ~ miR * Sex)
dds$miR #CON is reference level
dds$Sex #Male is reference level
dds$Sex.miR #MALE CON is reference level

#Pre-filter dds to remove poorly detected genes
#Minimal filtering rule- remove rows with at least 2 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
nrow(dds)

dds=DESeq(dds)
resultsNames(dds)

#Main effect of OE of miR (Male)
#What are the differences between OE vs control in males
res=results(dds, alpha = 0.05, contrast = c("miR", "30c", "CON"))
ix = which.min(res$pvalue) # most significant
res <- res[order(res$padj),] # sort
RES_MALE=as.data.frame(res)
RES_MALE2=RES_MALE%>%
  drop_na(padj) %>%
  filter(padj < 0.05)
summary(res)
barplot(assay(dds)[ix,],
        las=2, #direction of labels
        main=rownames(dds)[ ix  ]  )
#as.data.frame(res[1:20,]) %>% kbl() %>% kable_paper("hover", full_width = F)
#write.table(res,file="MALE.tsv",quote=FALSE,sep="\t")
#rownames(res) <- sapply( strsplit(rownames(res)," ") , "[[", 1)
#write.csv(RES_MALE2,"MALE OE v SCRAM.csv")

#Main effect plus interaction term
#Measuring the effect of OE in females
res2=results(dds, alpha=0.05, list(c("miR_30c_vs_CON","miR30c.SexFEMALE")))
ix = which.min(res2$padj) # most significant
res2 <- res2[order(res2$padj),] # sort
RES_FEMALE=as.data.frame(res2)
RES_FEMALE2=RES_FEMALE%>%
  drop_na(padj) %>%
  filter(padj < 0.05)
summary(res2)
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
#write.table(res2,file="FEMALE.tsv",quote=FALSE,sep="\t")
# Taking into account the overexpression what is the difference between Males and females
#MALE 30a VS FEMALE 30a
resultsNames(dds)
res4 = results(dds, alpha = 0.05, list( c("Sex_FEMALE_vs_MALE","miR30c.SexFEMALE") ))
ix = which.min(res4$padj) # most significant
res4 <- res4[order(res4$padj),] # sort
summary(res4)
RESDF4=as.data.frame(res4)
RESDF4.2=RESDF4%>%
  drop_na(padj)
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
as.data.frame(res4[1:20,]) %>% kbl() %>% kable_paper("hover", full_width = F)
#write.table(res4,file="OE COMPARISON.tsv",quote=FALSE,sep="\t")

# The difference in response to OE between male and females
# Is the effect of OE different across sexes (Interaction Term)
res5 = results(dds, alpha=0.05, name="miR30c.SexFEMALE")
ix = which.min(res5$padj) # most significant
res5 <- res5[order(res5$padj),] # sort
summary(res5)
as.data.frame(res5[1:20,]) %>% kbl() %>% kable_paper("hover", full_width = F)
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
# This gene is up-regulated by OE in males but down-regulated in females
#write.table(res5,file="INTERACTION.tsv",quote=FALSE,sep="\t")

#MALE SCRAM vs FEMALE SCRAM
#What is the difference between male and female without OE (i.e differences at baseline)
#CON (SCRAM) is reference so we just need to specify sex
res3=results(dds, alpha= 0.05, contrast = c("Sex", "MALE", "FEMALE"))
ix = which.min(res3$padj) # most significant
res3 <- res3[order(res3$padj),] # sort
summary(res3)
SCRAM=as.data.frame(res3)
SCRAM2=SCRAM%>%
  drop_na(padj) %>%
  filter(padj < 0.05)
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
#as.data.frame(res3[1:20,]) %>% kbl() %>% kable_paper("hover", full_width = F)
#write.table(res3,file="SCRAM comparison.tsv",quote=FALSE,sep="\t")

### FIGURES 
#####
#NORMALISE DATA FOR VISUALISATION PURPOSES (DO NOT USE FOR ANALYSIS)
NormDat=vst(dds, blind = F)
head(NormDat)
NormDat2= assay(NormDat)
#Intersect genes that were DE between males and females
RES_FEMALE2=rownames_to_column(RES_FEMALE2)
RES_MALE2=rownames_to_column(RES_MALE2)

#PLOTS
tiff("Male 30c.tiff", units = "in", width=10, height=7, res=600)
maplot <- function(de,contrast_name) {
  sig <-subset(de, padj < 0.05 )
  up <-rownames(subset(de, padj < 0.05 & log2FoldChange > 0))
  dn <-rownames(subset(de, padj < 0.05 & log2FoldChange < 0))
  GENESUP <- length(up)
  GENESDN <- length(dn)
  DET=nrow(de)
  SUBHEADER = paste(GENESUP, "Up-regulated, ", GENESDN, "Down-regulated, ", DET, "Detected")
  ns <-subset(de, padj > 0.05 )
  plot(log2(de$baseMean),
       de$log2FoldChange, 
       xlab="Log2 Base Mean", ylab="Log2 Fold-Change",
       pch=19, cex=1, col="dark gray",
       main=contrast_name, cex.main=2,
       cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
  points(log2(sig$baseMean),sig$log2FoldChange,
         pch=19, cex=1.5, col="red")
  mtext(SUBHEADER,cex = 1.5)
}
maplot(res,"DGE with miR-30c OE in Males")
dev.off()


## MAKE SURE TO REMOVE NA from DF (i.e. those with zero counts or no p val) for VOLCANO AND HEATMAP TO WORK
make_volcano <- function(de,name) {
  sig <- subset(de,padj<0.1)
  N_SIG=nrow(sig)
  N_UP=nrow(subset(sig,log2FoldChange>0))
  N_DN=nrow(subset(sig,log2FoldChange<0))
  DET=nrow(de)
  HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn", DET, "detected")
  plot(de$log2FoldChange,-log10(de$pval),cex=0.5,pch=19,col="darkgray",
       main=name, xlab="log2 FC", ylab="-log10 pval")
  mtext(HEADER)
  grid()
  points(sig$log2FoldChange,-log10(sig$pval),cex=0.5,pch=19,col="red")
}
make_volcano(res2,"DGE with miR-30c OE in males")

dev.off()

#HeatMap of SCRAM sex comparison
#HEATMAP.2
df2 <- as.data.frame(row.names(ss30a))
df2$type <- c("F1", "F2","F3","M1", "M2", "M3")
df2 <- df2$type[order(df2$type),]
labRow <- df2[match(row.names(), df2$'row.names(ss30a)') ]
heatmap.2(res2)

tiff("SCRAM comparison heatmap.tiff")
SCRAM_M= as.matrix(m30a.3)

tiff("myheatmap.tiff", units="in", width=7, height=5, res=600)
par(mar=c(1,1,1,1))
SCRAM2=make_heatmap <- function(de,name,ss30a,mx,n=30){
  colfunc <- colorRampPalette(c("blue", "white", "red"))
  values <- as.numeric(ss30a$Sex)
  f <- colorRamp(c("yellow", "orange"))
  rr <- range(values)
  svals <- (values-rr[1])/diff(rr)
  colcols <- rgb(f(svals)/255)
  mxn <- mx/rowSums(mx)*1000000
  x <- mxn[which(rownames(mxn) %in% rownames(head(de,n))),]
  heatmap.2(as.matrix(m30a.3),
            trace="none",
            col=colfunc(25),
            scale="row", 
            key=F,
            margins = c(8,9),  #first number is height, second is width, the smaller the bigger
            # The first number is the bottom margin, and the second is the left margin of the plot
            cexRow=1, 
            cexCol = 1,
            srtCol = 30,
            labRow = F,
            labCol = c("Female", "Female","Female","Male", "Male", "Male"),
            lwid=c(0.5,5), #make column of dendrogram and key very small and other column very big 
         #   lhei=c(0.2,5), #make row of key and other dendrogram very small and other row big.
            main=paste("DGE between male and female scramble condition"),
            ColSideColors = colcols)

}
make_heatmap(SCRAM2,"Control vs miR-30a overexpression",ss30a, m30a.3)

dev.off()

# Pathway Enrichment
##########
Female= read_excel("ORA genes.xlsx",
               sheet = 4,#Sheet 1 SCRAM, sheet 2 30a-bg, sheet 3 30a MALE, sheet 4 30a Female
               trim_ws = T)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
f=bitr(Female$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
m=bitr(Male$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
list=list(f=f$ENTREZID,m=m$ENTREZID)
xx= compareCluster(geneClusters = list, fun = enrichKEGG)
k <- setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

xx <- compareCluster(Gene ~ f+m,
                     data=DE_GSE8057, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")],
                     TERM2NAME=wp[,c("wpid", "name")])


BG= read_excel("ORA genes.xlsx",
                sheet = 2, #Background
                trim_ws = T)
gene_BG <- BG$ENSEMBL

GO <- enrichGO(gene         = SCRAM,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont           = "BP",
               pAdjustMethod = "BH",
               universe = gene_BG, 
               qvalueCutoff  = 0.05)

GO_RES=GO@result

## use simplify to remove redundant terms
ego = clusterProfiler::filter(GO, qvalue <0.05)
ego2 <- mutate(ego, richFactor = Count / as.numeric
               (sub("/\\d+", "", BgRatio)))
ego3 <- simplify(ego2, cutoff=0.7, by="qvalue", select_fun=min)


go=ggplot(ego3, showCategory = 15,
       aes(richFactor,
           fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=qvalue, size = Count)) +
  scale_color_gradientn(colours=c("red", "blue"),trans = "log10",
  guide=guide_colorbar(reverse=TRUE,
                       order=1)) +
  scale_size_continuous(range=c(2, 10)) +
 DOSE:: theme_dose(15) +
  xlab("Rich Factor") +
  ylab(NULL) +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0, face="bold")) +
  theme(axis.text.y = element_text(size = 12))+
  ggtitle("Top 15 Enriched GO Pathways: Biological Processes")

  ggpubr::ggarrange(go)+
  tiff('SCRAM GO PATHWAYS.tiff', width = 15, height = 10, units = 'in', res=600)
dev.off()
