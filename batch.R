library(affyPLM)
library(Biobase)

MLL.B = readRDS('MLLB.rds')
eset = readRDS('eset.rds')
sampleNames(MLL.B) = letters[1:21]

scanDate = protocolData(MLL.B)$ScanDate
scanDate = gsub(" .*", "", scanDate)
scanDate = as.Date(scanDate, "%m/%d/%y")

# Group together sampleswith the same date

minscan = min(scanDate)
days = scanDate- minscan

# Just group samples by day
batch = days

treatment = eset$characteristics_ch1.5

# HIERARCHICAL CLUSTERING
# DENDROGRAM
d = as.dist(1 - cor(exprs(eset), method="spearman"))
sampleClustering = hclust(d)
sampleDendrogram = as.dendrogram(sampleClustering, hang=0.5)
names(batch) = sampleNames(eset)
treatment = as.character(treatment)
names(treatment) = sampleNames(eset)
sampleDendrogram = dendrapply(sampleDendrogram, function(x, batch, labels) {
    ## for every node in the dendrogram if it is a leaf node
    if (is.leaf(x)) {
        attr(x, "nodePar") = list(lab.col = as.vector(batch[attr(x, "label")]))  ## color by batch
        attr(x, "label") = as.vector(labels[attr(x, "label")])  ## label by treatment
    }
    x
}, batch, treatment)

plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

# MULTIDIMENSIONAL SCALING (MDS)
cmd = cmdscale(as.dist(1 - cor(exprs(eset), method = "spearman")))
plot(cmd, type = "n")
text(cmd, treatment, col = batch, cex = 0.9)
legend("topleft", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)

