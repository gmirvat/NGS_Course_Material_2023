
setwd("/mnt/beegfs/workshop/<SSO>/results/")

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

## GSEA with Gene Ontology Terms

All_genes <- read.table("counts/DEGs_5uMaza_treatment_All.csv", sep=",", header =T)

logFC_all <- All_genes$log2FoldChange
names(logFC_all) <- All_genes$entrez

logFC_all_sorted <- sort(logFC_all, decreasing = TRUE)

# results might vary depending on the laptop you are using.
set.seed(123456)

gsea_go <- gseGO(geneList     = logFC_all_sorted,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 1, #try 1 instead of 0.05 # # padj cutoff value
                 verbose      = FALSE)


dim(gsea_go)


# use the function gseaplot() to visualize the results of the GSEA. 


# GSEA results for the 2nd set
pdf(file = "Pathway_GSEA_BP.pdf", width=9, height = 9)
gseaplot(gsea_go, geneSetID = 1, title = gsea_go$Description[2])

dev.off()
