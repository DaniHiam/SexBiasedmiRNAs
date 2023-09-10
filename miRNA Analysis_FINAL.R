setwd("C:/Users/dhiam/OneDrive - Deakin University/PostDoc_VU/miRNA_project")

library(tidyverse)
library(readxl)
library(limma)
library(superheat)
library(reshape2)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(ggrepel)
library('variancePartition')
library('edgeR')
library('BiocParallel')
library(mice)
library(miceadds)
library(pheatmap)
library(cluster)
library(purrr)
library(multiMiR)
library(clusterProfiler)
library("org.Hs.eg.db")
library(enrichplot)
library(TissueEnrich)
library("corrplot")
library("tableone")
library("PCAtools")
library(broom)
library(VIM)
library(naniar)
library(simputation)
library(skimr)
library(miRBaseConverter)

####################
# Participant Characteristics 
Pheno= read_excel("Pheno.xlsx",
                  sheet = 1,  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                  trim_ws = TRUE)
#RELEVEL
Pheno$Sex <- factor(Pheno$Sex,
                    levels = c( "Male", "Female"))
Pheno$Timepoint <- factor(Pheno$Timepoint,
                          levels = c("PRE","3HP"))
Pheno$Sex.Time <- factor(Pheno$Sex.Time,
                         levels = c("MALE.PRE", "MALE.3HP", "FEMALE.PRE", "FEMALE.3HP"))

Pheno2 = Pheno %>% column_to_rownames(var="FULL_CODE")
Pheno3= Pheno2 %>%
  dplyr:: select(2:13)

############## TABLE ONE ##################
#Create a variable list which we want in Table 1
Base_Pheno= Pheno %>%
  filter(Timepoint != "3HP")
listVars <- c("Age", "BMI", "Wpeak", "VO2", "LT", "TT", "fT", "SHBG", "Estrogen")
#Define categorical variables
catVars <- c("Sex")
#Total Population
table1 <- CreateTableOne(vars = listVars, 
                         factorVars = catVars,
                         strata = "Sex",
                         data = Base_Pheno)
table1 <- print(table1)
tab1 <- print(table1,  exact = "stage", 
              quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
#t.test(VO2 ~ Sex, data = Base_Pheno)
write.csv(tab1, file = "Table 1 Participant characteristics.csv")

############ PRE-PROCESSING ############
miR= read_excel("miRNA_RAW DATA_GlobalNorm.xlsx", 
                sheet = 2, 
                na = c("","NA", "?", "#VALUE!", "Undetermined"),
                trim_ws = TRUE)
colnames(miR)
# Identify duplicates (U6, RNU44 and RNU48 and calculate mean)
m2=aggregate(miR$Cq, by=list(ID=miR$ID,miRNA=miR$miRNA),data=miR,FUN=mean)

#Are there any duplicates? After aggregating there shouldn't be
Dup= m2 %>%
  dplyr::group_by(miRNA, ID) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

#Convert to wide
m3= pivot_wider(data = m2,
                      id_cols = miRNA,
                      names_from = ID,
                      values_from = "x")

#Exclude miRNAs where average of ALL samples are >35
HighCQ= function(m3) {
  m3[rowMeans(m3[,-1], na.rm=TRUE) < 35,]
}
m4=HighCQ(m3)
#Remove the rows that has miRNAs mean <35 

m5= m4 %>%
  filter(miRNA != "NA")

#Where individual cq is above 35 replace with NA
m6 = m5%>% 
  mutate(across(where(is.numeric), 
                    function(x) ifelse(x >= 35, NA, x)))

#Where individual cq is below 15 replace with NA 
#(biologically impossible- based on myomiR 133a cq)
m7 = m6%>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 15, NA, x)))
                
#20% NAs in x sample then exclude that miRNA
delete.na <- function(x, n=18) {
  x[rowSums(is.na(x)) <= n,]
}
m8= delete.na(m7) 

#Calculate Delta CT (0.5^x)(10^10)
DeltaFunc =  function(x) {
  return(0.5^x*10^10)
  }
Delta=data.frame(sapply(m8 [,-1],DeltaFunc, USE.NAMES = T))
sort=Pheno$FULL_CODE
m9= Delta[,sort]
miR.names= m8$miRNA
#DeltaFunc removes miRNA names therefore need to cbind
df

df_final=cbind(miR.names,m9)
#write.csv(df_final, "Processed Data before imputation.csv")

#PCA of data with missing NA
nb <- estim_ncpPCA(M, ncp.max=5)
nb
comp <- imputePCA(M,
                  ncp=nb$ncp,
                  scale=TRUE)
M_norm_t=t(scale(comp$completeObs))
res.pca = PCA(M_norm_t, graph = FALSE)
#coloured by Sex 
fviz_pca_ind(res.pca, 
                 geom = "point",
                 pointsize = 2,
                 habillage=Pheno$Sex,
                 addEllipses=T,
                 ellipse.level=0.95)+
  labs(  scale_color_brewer(palette="Set1") )+
  theme_minimal()

M_noNA= na.omit(M)
#Need to impute miRNA go from 308 to 111 miRs

#Imputation
#############
# IMPUTATION
Overview=skimr::skim(M)
hist(Overview$numeric.mean)
vis_miss(M)
gg_miss_var(M)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
ID_miss=as.data.frame(apply(M,2,pMiss)) #features (columns)
miR_miss=as.data.frame(apply(M,1,pMiss)) # rows
aggr_plot <- aggr(M, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(M), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

# KnN imputation with VIM package
# The optimal K value usually found is the square root of N, where N is the total number of samples.
sqrt(88)
M_KNN= df_final %>%
  kNN(k=9)
Overview2=skimr::skim(M_KNN)
hist(Overview2$numeric.mean)
M2_KNN=cbind(miR.names, M_KNN)
#write.csv(M2_KNN, "Matrix_KNN imputation.csv")

# MICE
set.seed(1)
m=5
MI <- mice(data = df_final, #dataset containing missing values
           m = m,
           print=F)

imp <- complete(MI,action="long")
#write.csv(imp, "M_MICE imputation.csv")
imp3 <- complete(MI,3)
#write.csv(imp3, "M_MICE imputation 3.csv")

############### START HERE ########
#READ IN miRNA matrix 
dat1=read_excel("miRNA_RAW DATA_GlobalNorm.xlsx",
                sheet = 4,  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                trim_ws = TRUE)
M=column_to_rownames(dat1, "miR.names")
# SHEET 3 = M before imp
# SHEET 4 = M with KNN **************

#Change array names to miRBASE ID
nameEDIT=read_excel("RESULTS_KNN LIMMA.xlsx",
                sheet = 2,  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                trim_ws = TRUE)
name.edit2 = merge(x = nameEDIT, y = dat1,
                 by.x='rowname', by.y='miR.names', all.y=T)
colnames(name.edit2)
M=name.edit2[,c(-1:-3, -5:-8)]
M2=column_to_rownames(M, "miRBaseID")
#M2=log(M2)

f=M$miRBaseID
miRNANames = f
version=checkMiRNAVersion(miRNANames, verbose = TRUE)
f=miRNA_NameToAccession(miRNANames, version = "v20")
colnames(f)
f2=column_to_rownames(f, "miRNAName_v20")

#Create ExpressionSet Object
eset=ExpressionSet(assayData = as.matrix(M2),
                   phenoData = AnnotatedDataFrame(Pheno3),
                   featureData = AnnotatedDataFrame(f2))

# Logging data allows you to increase the distance between small measurements and decreases the distance between large measurements 
plotDensities(eset, legend = F)

#log() is natural log (base 10)
#Log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, legend =F)
#Extract Matrix
Matrix =exprs(eset)
#Extract Pheno
P= pData(eset)
#Extract Anno
f=fData(eset)

# Identify any samples that appear as "outliers"
meltData2 <- melt(scale(Matrix))
q <- ggplot(meltData2, aes(factor(Var2), value)) 
q + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

# MDS and PCA plots
# gene.selection='common' makes 'plotMDS' perform a more traditional principal component analysis.  
plotMDS(eset, labels = pData(eset)[, 'Sex.Time'],
        gene.selection = "common")


# GROUP MEANS PARAMTERISATION. 
# The coefficients now represent the mean expression level for each of the groups of samples, so you need to specify a custom contrast to test for differences between groups

design <- model.matrix(~ 0 + Sex.Time  + BMI+  Wpeak, data= pData(eset)) 
#Make sure you have zero in model so you can run the comparisons "manually"

#The number of samples modeled by the coefficient should equal number of samples in each group
table(pData(eset)[,"Sex.Time"])
corfit <- duplicateCorrelation(eset,
                               design,
                               block = Pheno$ID) 
fit = lmFit(eset,
            design,
            block = Pheno$ID,
            correlation=corfit$consensus)

cont.matrix = makeContrasts(BASE="Sex.TimeMALE.PRE -  Sex.TimeFEMALE.PRE", 
                            POST="Sex.TimeMALE.3HP -  Sex.TimeFEMALE.3HP",
                            MALE="Sex.TimeMALE.3HP - Sex.TimeMALE.PRE", 
                            FEMALE="Sex.TimeFEMALE.3HP - Sex.TimeFEMALE.PRE", 
                            INTERACTION = "(Sex.TimeFEMALE.3HP - Sex.TimeFEMALE.PRE) - (Sex.TimeMALE.3HP - Sex.TimeMALE.PRE)",
                            levels = design)
cont.matrix
plotContrasts(cont.matrix)
fit2=contrasts.fit(fit, cont.matrix)
fit3=eBayes(fit2)
summary(decideTests(fit3))

stats=topTable(fit3, 
               coef = 'BASE',
               number = nrow(fit3), 
               sort.by = "none")
hist(stats[,"P.Value"])

#Volcano Plot
volcanoplot(fit3, highlight = 5, names = fit3$genes[,"Accession"])
#x-axis is the LogFC between contrasts and y-axis is the log-odds the more likely the gene is differentially expressed (The one highest on the y-axis).


########## PCA TOOLS #######
colnames(Pheno)
Pheno2 = Pheno %>% column_to_rownames(var="FULL_CODE")
Pheno3= Pheno2 %>%
  dplyr:: select(2:13)

#Match Matrix and Pheno order
all(colnames(M2) == rownames(Pheno3))

p=pca(M2, metadata=Pheno3, center = T,
      scale = F, removeVar=0.1)
SCREE=screeplot(p, components = getComponents(p, 1:5),
                    hline = 80, vline = 4, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = FALSE) +
  geom_label(aes(4.5, 70, label = '80% of variance explained', vjust = -1, size = 8))
SCREE

## What variables are important to add to model
CORR=eigencorplot(p, components = getComponents(p, 1:4),
             metavars = c("Age",
                          "BMI",
                          "Wpeak",
                          "VO2",
                          "LT",
                          "TT",
                          "fT",
                          "SHBG",
                          "Estrogen"),
            #col = c('white', 'red'),
             cexCorval = 1.0,
             fontCorval = 2,
             posLab = 'all', 
             rotLabX = 45,
             scale = T,
             main = "Correlation of phenotype variables with Principal components 1-4",
             cexMain = 1.2,
             plotRsquared = F,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('*', ''),
             signifCutpoints = c(0, 0.05, 1),
             returnPlot = FALSE)

library(cowplot)
library(ggplotify)
top_row <- plot_grid(SCREE,
                     ncol = 1,
                     labels = c('A'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 22,
                     align = 'h',
                     rel_widths = c(1.10))

bottom_row <- plot_grid(as.grob(CORR),
                        ncol = 1,
                        labels = c('B'),
                        label_fontfamily = 'serif',
                        label_fontface = 'bold',
                        label_size = 22,
                        align = 'h',
                        rel_widths = c(1.1))

p=plot_grid(top_row,bottom_row)
p
ggpubr::ggarrange(p, ncol = 1, nrow = 1)+
  tiff('PCA Analysis.tiff', width = 15, height = 5, units = 'in', res=600)

dev.off()


## PLOT results
M_row=rownames_to_column(M2)
M_row[which(M_row$rowname == "hsa-miR-628"),]
Pheno_with_miR <- cbind(Pheno,
                        M = as.numeric(M[225,]))

ggplot(Pheno_with_miR, aes(x=Timepoint, y=M, fill=Sex)) +
  geom_boxplot(lwd=1)+
  geom_point()+
  facet_wrap(~Sex)+
  scale_fill_brewer(palette="Accent")+
  theme(text = element_text(size=15),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "none",
        panel.background = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black", size=1),
        axis.title.x = element_blank())
ggpubr::ggarrange(miR30d, ncol = 1, nrow = 1)+
  tiff('miR30d.tiff', width = 5, height = 5, units = 'in', res=600)
dev.off()

library(miRBaseConverter)
SexDE= read_excel("RESULTS_KNN LIMMA.xlsx", 
                  sheet = 8, 
                  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                  trim_ws = TRUE)
SexDE2= as.data.frame(c(SexDE$Mature.Female))
SexDE3=na.omit(SexDE2)
names(SexDE3)[1] = "Mature"

miRNANames = SexDE3$Mature
version=checkMiRNAVersion(miRNANames, verbose = TRUE)
BASE=miRNA_NameToAccession(SexDE3$Mature, version = "v20")
BASE_PRE=miRNA_MatureToPrecursor(BASE$miRNAName_v20)
BASE2= cbind(BASE_PRE, BASE)
write.csv(BASE2, "miRNA Name conversion Female.csv")

### Background
# http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
all_results=topTable(fit3,
                      coef = 'INTERACTION',
                             number = nrow(M_log),
                             adjust.method = "BH",
                     #lfc = 1,
                             p.value = 1)
hist(all_results$P.Val)# Must be conducted on all miRs and RAW p val

#VOLCANO PLOT
  all_results = mutate(all_results, sig=ifelse(all_results$adj.P.Val<0.05, "FDR<0.05", "Not Sig"))
  all_results = mutate(all_results, coef=ifelse(all_results$t<0, "neg", "pos"))
  all_results = mutate(all_results, color=ifelse(all_results$adj.P.Val>0.05,"black",ifels(all_results$adj.P.Val<0.05&all_results$coef=="neg","blue","red")))
 VC_INTERACTION= ggplot(all_results, aes(logFC, -log10(P.Value))) +
    geom_point(aes(col=color), size=1.5)+
    scale_color_manual(values=c("black","#1F78B4", "#E31A1C"))+
    labs(x="DE miRNA (logFC) in exercise response between sexes",y="-log10(p-value)")+
    coord_cartesian(clip = 'off')+#for the one infinite point
    labs(text=element_text(family="Times", size = 12))+
    theme_classic()+
    theme(legend.position = "none")
VC_INTERACTION

all_results=topTable(fit3,
                     coef = 'BASE',
                     number = nrow(M_log),
                     adjust.method = "BH",
                     p.value = 1)
#VOLCANO PLOT
all_results = mutate(all_results, sig=ifelse(all_results$adj.P.Val<0.05, "FDR<0.05", "Not Sig"))
all_results = mutate(all_results, coef=ifelse(all_results$t<0, "neg", "pos"))
all_results = mutate(all_results, color=ifelse(all_results$adj.P.Val>0.05,"black",ifelse(all_results$adj.P.Val<0.05&all_results$coef=="neg","blue","red")))
VC_BASE= ggplot(all_results, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=color), size=1.5)+
  scale_color_manual(values=c("black","#1F78B4", "#E31A1C"))+
  labs(x="DE miRNA (logFC) between sexes at baseline",y="-log10(p-value)")+
  coord_cartesian(clip = 'off')+#for the one infinite point
  labs(text=element_text(family="Times", size = 12))+
  theme_classic()+
  theme(legend.position = "none")
VC_BASE

ggpubr::ggarrange(VC_BASE, VC_INTERACTION, ncol = 1, nrow =2)+ #labels = "AUTO")+
  tiff('Volcano Plots.tiff', width = 5, height = 10, units = 'in', res=600)
dev.off()

ggplot2::ggsave(file = "miRNA_Sex.Time Volcano .tiff", plot = VC_INTERACTION, 
                  width = 7, height = 7)
  dev.off() 
  
############## GENE TARGETS
  ###############
Sex.Time_DE= read_excel("RESULTS_KNN LIMMA.xlsx", 
                    sheet = 8, 
                    na = c("","NA", "?", "#VALUE!", "Undetermined"),
                    trim_ws = TRUE)  
colnames(Sex.Time_DE)

#MALE
Vector_MALE= c(Sex.Time_DE$Mature.Male)
Vector_MALE=na.omit(Vector_MALE)
TargetGenes= multiMiR::get_multimir(org = "hsa",
                                    mirna = Vector_MALE,
                                    table = "validated",
                                    summary = T)
validated=TargetGenes@data
genes_ensembl= c(validated$target_ensembl)
genes_2=unique(genes_ensembl) #character not integer?
write.csv(genes_2, file = "Validated_MALE.csv")

#BASE
Vector_base= c(Sex.Time_DE$Mature.BASE)
Vector_base=na.omit(Vector_base)
TargetGenes= multiMiR::get_multimir(org = "hsa",
                                   mirna = Vector_base,
                                   table = "validated",
                                   summary = T)
validated=TargetGenes@data
genes_ensembl= c(validated$target_ensembl)
genes_2=unique(genes_ensembl) #character not integer?
write.csv(genes_2, file = "Validated_BASE.csv")

#INTERACTION
Vector_interact= c(Sex.Time_DE$Mature.Interact)
Vector_interact=na.omit(Vector_interact)
TargetGenes= multiMiR::get_multimir(org = "hsa",
                                    mirna = Vector_interact,
                                    table = "validated",
                                    summary = T)
validated=TargetGenes@data
genes_ensembl= c(validated$target_ensembl)
genes_2=unique(genes_ensembl) #character not integer?
write.csv(genes_2, file = "Validated_INTERACT.csv")

#BACKGROUND
Vector_BG= c(Sex.Time_DE$Mature.BG)
Vector_BG=na.omit(Vector_BG)
TargetGenes= multiMiR::get_multimir(org = "hsa",
                                    mirna = Vector_BG,
                                    table = "validated",
                                    summary = T)
validated=TargetGenes@data
genes_ensembl= c(validated$target_ensembl)
genes_2=unique(genes_ensembl) #character not integer?
write.csv(genes_2, file = "Validated_BACKGROUND.csv")


# Pathway Enrichment
# PAPER ON VISUALISING ENRICHMENT PATHWAY
#https://www.cell.com/the-innovation/pdf/S2666-6758(21)00066-7.pdf 
##########
Val.BG= read.csv2("Validated_BACKGROUND.csv", sep = ",")
Val.BASE= read.csv2("Validated_BASE.csv", sep = ",")
Val.INTERACT= read.csv2("Validated_INTERACT.csv", sep = ",")
Val.FEMALE= read.csv2("Validated_FEMALE.csv", sep = ",")
Val.MALE= read.csv2("Validated_MALE.csv", sep = ",")
colnames(Val.BASE)
genes_BG <- Val.BG$x
genes_BASE <- Val.BASE$x
genes_INTERACT <- Val.INTERACT$x
genes_FEMALE <- Val.FEMALE$x
genes_MALE <- Val.MALE$x

### TISSUE ENRICH
#A tool to calculate tissue-specific gene enrichment
#https://www.bioconductor.org/packages/release/bioc/vignettes/TissueEnrich/inst/doc/TissueEnrich.html
#Genes from the Tissue Enriched, Group Enriched, and Tissue Enhanced groups are classified as tissue-specific genes.
gs_BASE=unique(genes_BASE)
gs2_BASE<-GeneSet(geneIds=gs_BASE,organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())

gs_INTERACT=unique(genes_INTERACT)
gs2_INTERACT<-GeneSet(geneIds=gs_INTERACT,organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())

gs_FEMALE=unique(genes_FEMALE)
gs2_FEMALE<-GeneSet(geneIds=gs_FEMALE,organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())

gs_MALE=unique(genes_MALE)
gs2_MALE<-GeneSet(geneIds=gs_MALE,organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())

gs_BG=unique(genes_BG)
gs2_BG<-GeneSet(geneIds=gs_BG,organism="Homo Sapiens",geneIdType=ENSEMBLIdentifier())


output<-teEnrichment(inputGenes = gs2_MALE,
                     rnaSeqDataset = 2, #GTEx portal
                     backgroundGenes = gs2_BG)

seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
View(enrichmentOutput)

#Retrieval of genes enriched in SM
seGroupInf<-output[[3]][["Muscle"]]
SM_genes<-data.frame(assay(seGroupInf))
SM_genes2 <- SM_genes$Gene
write.csv(SM_genes, "SM ENRICHED GENES_MALE.csv")

#### Gene ENRICHMENT 
##############
Genes= read_excel("RESULTS_KNN LIMMA.xlsx", 
                        sheet = 9, 
                        na = c("","NA", "?", "#VALUE!", "Undetermined"),
                        trim_ws = TRUE)  
colnames(Genes)
Genes_BASE= Genes$GENES_BASE_SM
Genes_INTERACT=Genes$GENES_INTERACTION_SM
Genes_FEMALE= unique(Genes$GENES_FEMALE_SM)
Genes_MALE= unique(Genes$GENES_MALE_SM)
Genes_BG = Genes$GENES_BACKGROUND

GO <- enrichGO(gene         = Genes_INTERACT,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont           = "BP",
               pAdjustMethod = "BH",
               universe = Genes_BG, #Background Genes
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1)
GOPathways=GO@result
#write.csv(GOPathways, "GO Pathway_MALE.csv")

ego2 <- simplify(GO, cutoff=0.7, by="p.adjust",
                 select_fun=min)

p1=barplot(ego2, showCategory=10, order=T,
        orderBy="GeneRatio")+
  ggtitle("A) Enriched biological processes")+
  theme_bw(base_size = 25)
p1

p2=barplot(ego2, showCategory=10, order=T,
           orderBy="GeneRatio")+
  ggtitle("B) Enriched biological processes")+
  theme_bw(base_size = 25)
p2

ggpubr::ggarrange(p1, p3, ncol = 2, nrow = 1)+
  tiff('Top 10 SM enriched pathways.tiff', width = 30, height = 10, units = 'in', res=600)
dev.off()


dotplot(ego2, showCategory=15,
        orderBy = "Count") + ggtitle("dotplot for ORA")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

#KEGG ORA PATHWAY
Genes= read_excel("RESULTS_KNN LIMMA.xlsx", 
                  sheet = 9, 
                  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                  trim_ws = TRUE)  
colnames(Genes)
Genes_BASE= Genes$GENES_BASE_SM
Genes_INTERACT=Genes$GENES_INTERACTION_SM
Genes_FEMALE= unique(Genes$GENES_FEMALE_SM)
Genes_MALE= unique(Genes$GENES_MALE_SM)
Genes_BG = Genes$GENES_BACKGROUND
?bitr
ENTREZ_BASE = bitr(Genes_BASE, fromType= "ENSEMBL"  , toType="ENTREZID", OrgDb="org.Hs.eg.db")
ENTREZ_INTERACT = bitr(Genes_INTERACT, fromType= "ENSEMBL"  , toType="ENTREZID", OrgDb="org.Hs.eg.db")

search_kegg_organism('hsa', by='kegg_code')
eg2np <- bitr_kegg(ENTREZ_BASE$ENTREZID, fromType='kegg', toType='ncbi-geneid', organism='hsa')
kk <- enrichKEGG(gene         = eg2np,
                 organism     = 'hsa',
                 keyType = "ncbi-geneid",
                 pvalueCutoff = 0.05)

WikkiP= enrichWP(ENTREZ_INTERACT$ENTREZID, organism = "Homo sapiens") 
WP= WikkiP@result
barplot(WikkiP, showCategory=15, order=T,
          orderBy="GeneRatio")+
  ggtitle("Top 10 enriched WikiPathways")+
  theme_bw(base_size = 25)

library(ReactomePA)
REACTOME <- enrichPathway(gene=ENTREZ_BASE$ENTREZID, 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          readable=TRUE)
Reactome=REACTOME@result
p3=barplot(REACTOME, showCategory=10, order=T,
        orderBy="GeneRatio")+
  ggtitle("B) Enriched Reactome Pathways")+
  theme_bw(base_size = 25)
p3


#############
#miRNA database 
Genome_COORD= read_excel("miRBase_Genome Location.xlsx", 
                  sheet = 4, 
                  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                  skip = ,
                  trim_ws = TRUE)

glimpse(Genome_COORD)

#Keep only precursor transcript
Coord_primary= Genome_COORD %>%
  filter(Mature_OR_Precurser == "miRNA_primary_transcript")

SexDE= read_excel("RESULTS_KNN LIMMA.xlsx", 
                  sheet = 8, 
                  na = c("","NA", "?", "#VALUE!", "Undetermined"),
                  trim_ws = TRUE)
SexDE2= as.data.frame(c(SexDE$Precursor.Interact))
SexDE3=na.omit(SexDE2)
names(SexDE3)[1] = "Precursor"
#Location of Sex-biased miRNAs
Chrom_BASE = left_join(SexDE3, Genome_COORD, by = "Precursor" )
SexFreq = as.data.frame(table(Chrom_BASE$chrom))
write.csv(SexFreq, "INTERACT Chromosome Freq.csv")

#Location of non sex-biased miRNAs
Chrom_BG = as.data.frame(c(SexDE$Precursor.BG))
Chrom_BG=na.omit(Chrom_BG)
names(Chrom_BG)[1] = "Precursor"
Chrom_BG2 = left_join(Chrom_BG, Genome_COORD, by = "Precursor" )
BG = as.data.frame(table(Chrom_BG2$chrom))
write.csv(BG, "BG Chromosome Freq.csv")

#Read in Chromosome Freq
ChromF = read_excel("Chromosome Freq.xlsx", sheet = 3,
             #       na = c("","NA", "?", "#VALUE!", "Undetermined"),
                    trim_ws = TRUE)
ChromF=ChromF %>% 
  column_to_rownames(var = "Chromosome")
ChromF= ChromF[c(-1,-4)]
glimpse(ChromF)

library("gplots")
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
dt <- as.table(as.matrix(ChromF))
View(dt)
chisq=chisq.test(dt) #Make there is only numbers not characters (i.e. rownames)
round(chisq$residuals, 3)
corrplot(chisq$residuals, is.cor = FALSE,
         tl.col="black",cl.offset=1,cl.ratio=0.7, cl.pos = 'n')

contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
colour <- colorRampPalette(brewer.pal(8, "RdBu"))(8)
corrplot(contrib, is.cor = FALSE, col = colour,
         tl.col="black",cl.offset=1,cl.ratio=0.7, cl.pos = 'n')+
colorlegend(xlim=c(3,4), ylim=c(1,23), colour, c(seq(0,20,2)), align="l", vertical=T, addlabels=TRUE)
# printing the p-value
chisq$p.value

fisher.test(ChromF, simulate.p.value = T)

#AUTO V Sex C
ChromF = read_excel("Chromosome Freq.xlsx", sheet = 6,
                                       trim_ws = TRUE)
ChromF=ChromF %>% 
  column_to_rownames(var = "Type")
ChromF2= ChromF[,-3]
mosaicplot(ChromF2,
           main = "Mosaic plot",
           color = TRUE
)
chisq.test(ChromF2)$expected
fisher.test(ChromF2, simulate.p.value = T)

CHROM_B=ggplot(data=ChromF, aes(x=Type, y=Prop, fill=Type)) + 
  geom_bar(stat="identity", position = "dodge")+
  theme_classic()+
  scale_fill_manual(values=c("#E31A1C","#1F78B4"))+
  xlab("Chromosome Type")+
  ylab("% of total miRNAs")+
  ggtitle("A) Proportion of sex-biased miRNAs at baseline")+
  ylim(0, 100)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))+
   theme(strip.background = element_rect(
    color="gray", fill="gray", size=1.5, linetype="solid"  ))

CHROM_INTER=ggplot(data=ChromF, aes(x=Type, y=Prop, fill=Type)) + 
  geom_bar(stat="identity", position = "dodge")+
  theme_classic()+
  scale_fill_manual(values=c("#E31A1C","#1F78B4"))+
  xlab("Chromosome Type")+
  ylab("% of total miRNAs")+
  ggtitle("B) Proportion of sex-biased miRNAs after exercise bout")+
  ylim(0, 100)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))+
  theme(strip.background = element_rect(
    color="gray", fill="gray", size=1.5, linetype="solid"
  ))
CHROM_INTER

ggpubr::ggarrange(CHROM_B, CHROM_INTER, ncol = 2, nrow = 1)+
  tiff('CHROM Location_Side.tiff', width = 12, height = 5, units = 'in', res=600)
dev.off()

########
# HORMONE ANALYSIS
# FEMALES
M_AR2= M_log %>%
  dplyr::select(-c(ends_with("3HP")))
M_AR2= M_AR2 %>%
  dplyr::select(-c("SG141_PRE", "FI5_PRE")) #No hormone data
M_females= M_AR2 %>%
  dplyr::select(starts_with("FI"))

Pheno_female= Pheno %>%
  dplyr::filter(Pheno$Sex == "Female") 
Pheno_female2= Pheno_female %>%
  dplyr::filter(Pheno_female$Timepoint == "PRE")

hist(Pheno_female$fT)

design <- model.matrix(~fT + Age, Pheno_female2)
fit = lmFit(M_females,
            design)
fit2 <- eBayes(fit)
summary(decideTests(fit2))
results=topTable(fit2, 
                 coef='fT',
                 adjust.method = "BH",
                 number=nrow(M_females),
                 p.value = 1)
View(results)
res_export=as.data.frame(results)
write.csv(res_export, "Females BASELINE.csv")


## MALES
M_males= M_AR2 %>%
  dplyr::select(starts_with("SG"))

Pheno_male= Pheno %>%
  dplyr::filter(Pheno$Sex == "Male") 
Pheno_male2= Pheno_male %>%
  dplyr::filter(Pheno_male$Timepoint == "PRE")
Pheno_male2= Pheno_male2 %>%
  dplyr::filter(Pheno_male2$ID != "SG141") #No hormone data

design <- model.matrix(~fT + Age, Pheno_male2)
fit = lmFit(M_males,
            design)
fit2 <- eBayes(fit)
summary(decideTests(fit2))
results=topTable(fit2, 
                 coef='fT',
                 adjust.method = "BH",
                 number=nrow(M_males),
                 p.value = 1)
View(results)
#res_export=as.data.frame(results)
#write.csv(res_export, "Males BASELINE.csv")

## PLOT results
M_row= rownames_to_column(M_males)
M_row[which(M_row$rowname == "hsa-miR-613"),]
Pheno_with_miR <- cbind(Pheno_male2,
                        miR = as.numeric(M_males[248,]))

miR=ggplot(Pheno_with_miR, aes(x=fT, y=miR)) +
  geom_point()+
  geom_smooth(method = "lm")+
labs(y=paste("Expression level of miR-613 (AU)"),
     x= paste("Free Testosterone (pM)"))+
  #ggtitle("hsa-miR-613 Expression in males is correlated with fT")+
  theme(text = element_text(size=12), 
        plot.title = element_text(face="bold", hjust = 0.5),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text= element_blank(), 
        axis.line = element_line(colour = "black", size=1))
miR

ggpubr::ggarrange(miR, ncol = 1, nrow = 1)+
  tiff('miR613_fT.tiff', width = 5, height = 5, units = 'in', res=600)

dev.off()

# Fitness Modelling 
# Timepoint and baseline fitness measures, does the fitness of the individual influence the miR response to training. 
design <- model.matrix(~ Sex+ Timepoint*Wpeak , data= pData(eset)) 

corfit <- duplicateCorrelation(eset,
                               design,
                               block = Pheno$ID) 
fit = lmFit(eset,
            design,
            block = Pheno$ID,
            correlation=corfit$consensus)

fit3=eBayes(fit)
summary(decideTests(fit3))

stats=topTable(fit3, 
               coef = 'Timepoint3HP:Wpeak',
               number = nrow(fit3), 
               sort.by = "none")
par(oma=c(1,1,1,1),
    mar=c(1.2,1.2,1.2,1.2))
hist(stats[,"P.Value"],
     main = "Hist of unadj p-values for baseline Peak Power",
     xlab="Unadjusted p values",
     ylab="Count",
     cex.main=1)
#Histogram of unadjusted p-values looking at the effect of baseline Wpeak changes in miRNA response to acute exercise bout.

view(stats)
results2= tibble::rownames_to_column(stats)
write.csv(results2, "miR_PP.csv")