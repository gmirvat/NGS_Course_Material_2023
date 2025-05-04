
#NGS workshop
#Data: 17th May 2023...2nd hour 
#Mirvat Surakhy

#In this script we will be doing DGE analysis using Deseq2 package and view the output using
#some plots
#References:
#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#in this  [https://yulab-smu.github.io/clusterProfiler-book/index.html]

setwd( "/mnt/beegfs/workshop/<SSO>/results/" )


library(dplyr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(apeglm)

#Read the count data generated in the previous code
#Load the Sample information (metadata). This is extracted from the SRA table when we #downloaded the data 

counts <- read.table(file = "counts/Count_fromtximport_Salmon.txt", header= TRUE, check.names = F)
metadata <- read.csv(file= "/mnt/beegfs/workshop/DGE_results_codes/SampleInfo_ngs.csv", row.names = 1)
metadata$condition <- factor(metadata$condition,
                             levels = c("Control", "Treatment"))

#check if the sample name is matching the metadata
all(rownames(metadata)%in% colnames(counts))

all(rownames(metadata)== colnames(counts))

#If the order of rows and columns is  not the same, try do the following
counts<- counts[, row.names(metadata)]  

###########DGE ######################
######## design the matrix
design <- as.formula(~condition)
model<- model.matrix(design, data= metadata)
keep <- rowSums(counts)>5 
countdata<- counts[keep,]
countdata<- as.matrix(countdata)
#### Create DESeq2Dataset object
dds.raw<- DESeqDataSetFromMatrix(countData = countdata, 
                                    colData = metadata,    
                                   design = design)  
#Perform the differential gene expression
dds <- DESeq(dds.raw)

### a. Obtain the results from DESeq object
##alpha indicates the value of padj. By default the argument alpha is set to 0.1
res_05<- results(dds,alpha= 0.05)

#lists the coefficients  #Use the out put in the lfcShrink
resultsNames(dds) 

resLFC_05 <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm", res= res_05)

#######Plot and view the results 
#####MA plot
##Plot MA plot before log fold shrinkage.
pdf(file ="plotMA_plot.pdf", width=5)
plotMA(res_05, ylim=c(-3,3))
 dev.off()

##Plot MA plot after log fold shrinkage. Can you see the effect of lfcShrink
pdf(file = "plotMA_plot_res_05.pdf", width=5)
plotMA(resLFC_05, ylim=c(-3,3))
dev.off()

########### Volcano Plot ####################
#Prepare the data
DEG<- as.data.frame(resLFC_05) # assign the results to data frame degenes

# assign symbol(common names) and ENTREZ id according to ENSEMBL id
DEG$symbol<- mapIds(org.Hs.eg.db,
                          keys= rownames(DEG),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

DEG$entrez<- mapIds(org.Hs.eg.db,
                          keys= rownames(DEG),
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          mutiVals= "first")

# remove genes that don't have a common name and those with duplicated gene name
DEG_symbol<- DEG[is.na(DEG$symbol)== FALSE,]
dim(DEG_symbol)

DEG_symbol<- DEG_symbol[!duplicated(DEG_symbol$symbol),]

DEG05_symbol<- subset(DEG_symbol, padj< 0.05  &abs(log2FoldChange)>1)

write.csv(DEG_symbol, "counts/DEGs_5uMaza_treatment_All.csv")

write.csv(DEG05_symbol, "counts/DEGs_5uMaza_treatment_significant.csv")


##################VolcanoPlot############

DEG_symbol$ENSEMBL.ID=row.names(DEG_symbol)
row.names(DEG_symbol) <- DEG_symbol$symbol

pdf(file = "Volcano_plot.pdf", width=7)
EnhancedVolcano(DEG_symbol,
    lab = rownames(DEG_symbol),
    x = 'log2FoldChange',
    y = 'padj', 
labSize = 3.0,
  shape = c(1, 4, 17, 18),
 col=c('black','gray','blue', 'red3'))
dev.off()
#####Principle component analysis (PCA)

#first extract log normalised counts from the dds object then plot the PCA

rld <- rlog(dds, blind=TRUE)

pdf(file = "plotPCA_plot.pdf", width=5)

plotPCA(rld, intgroup="condition")

dev.off()

sessionInfo()

