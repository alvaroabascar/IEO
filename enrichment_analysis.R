library(GSEABase)
library(affyPLM)
library(Biobase)
library(genefilter)
library(hgu133plus2.db)

glioData = readRDS('rma_glioData.rds')
eset = readRDS('eset.rds')
annotation(eset)<-'hgu133plus2.db'

g28_ids = as.character(read.csv('cluster_g28.csv', header=FALSE)$V3)

g28 = GeneSet(AnnotationIdentifier("hgu133plus2.db"), geneIds=g28_ids, setName="G28")
details(g28)

cluster_collection = GeneSetCollection(g28)
cluster_collection
length(cluster_collection)
head(names(cluster_collection))

cluster_collection = mapIdentifiers(cluster_collection, AnnoOrEntrezIdentifier(annotation(eset)))
Im = incidence(cluster_collection)


GSeset = gsva(eset, cluster_collection, min.sz=10, max.sz=500, verbose=FALSE)$es.obs
