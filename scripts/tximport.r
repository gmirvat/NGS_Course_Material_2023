#NGS workshop
#Data 17th May 2023
#Mirvat Surakhy and  James Carrington
####This code is for changing the TPM values to counts that can be used for  differential gene expression using deseq2
#The codes are based on the following tutorials
#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

library(tximport)
library(dplyr)


setwd ("/mnt/beegfs/workshop/<SSO>/results/")

## List all directories containing data  
samples <- list.files(path = "./salmon_quants", full.names = T, pattern="_sample_quant$")
samples

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")
files

##assign a shorter name for each element
names(files) <- list.files("salmon_quants")


## read the annotation file 

tx2gene <- read.csv("/mnt/beegfs/workshop/DGE_results_codes/tx2gene_ens109.csv")

# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("TXNAME", "GENEID")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)
head(txi[["counts"]])
counts <- txi$counts %>% round()
counts <-data.frame(counts)

#change the column names
colnames (counts) <- c("5aza_rep1","5aza_rep2", "DMSO_rep1","DMSO_rep2" )
##view the output, can you notice the change 
head(counts)


write.table(counts, "counts/Count_fromtximport_Salmon.txt", sep ="\t", quote = F)

sessionInfo()



