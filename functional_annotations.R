library(affyPLM)
library(Biobase)
library(org.Hs.eg.db)
library(hgu133plus2.db)

glioData = readRDS('rma_glioData.rds')
eset = readRDS('eset.rds')

GO_genes = select(hgu133plus2.db, keys=featureNames(glioData), columns=c("SYMBOL", "GO"))


