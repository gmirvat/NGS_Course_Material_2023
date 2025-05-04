
#Pathway Analysis
#Data 17thMay2023 ...3rd hour
#Mirvat Surakhy

#In this script we will be doing over represenation analysis using clusterProfiler package 
#This is based on the following tutorial

setwd("/mnt/beegfs/workshop/<SSO>/results/" )

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

######Pathway Analysis##################

DEG05_symbol <- read.csv("counts/DEGs_5uMaza_treatment_significant.csv"  )

geneset <- as.character(DEG05_symbol$entrez)

#GO over-representation analysis
ora_go <- enrichGO(gene          = geneset,
                   universe      = NULL,            # all available genes in db 
                   OrgDb         = org.Hs.eg.db,    # Hs: homo sapiens
                   ont           = "BP",            # One of MF, BP, CC*
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)


#The enriched pathway are saved in the results of the ora_go object
head(ora_go@result)

#Visualize enriched GO terms

pdf(file = "Pathway_GO_BP_dotplot.pdf", width=7)
dotplot(ora_go, showCategory = 10)+ ggtitle("Dot plot for ORA BP")
dev.off()

pdf(file = "Pathway_GO_BP_barplot.pdf", width=7)
barplot(ora_go, showCategory = 10)+ ggtitle("Bar plot for ORA BP")
dev.off()



logFC_de <- DEG05_symbol$log2FoldChange
names(logFC_de) <- DEG05_symbol$entrez
pdf(file = "Pathway_GO_BP_cnetplot.pdf", width=9, height = 9)
cnetplot(ora_go, foldChange = logFC_de, circular = FALSE)+ 
  ggtitle("CNE plot for ORA BP")

dev.off()


sessionInfo()


