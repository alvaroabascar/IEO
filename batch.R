library(affyPLM)
library(Biobase)
library(corpcor)

# glioData = ReadAffy(celfile.path = "./data")
# eset = readRDS('full_eset.rds')
glioData = readRDS('glioData.rds')
eset = readRDS('eset.rds')
sampleNames(glioData) = letters[1:21]

scanDate = protocolData(glioData)$ScanDate
scanDate = gsub(" .*", "", scanDate)
scanDate = as.Date(scanDate, "%m/%d/%y")

# Group together sampleswith the same date

minscan = min(scanDate)
days = scanDate- minscan

batch_days = data.frame(row.names=sort(unique(days)), batch=c(1:length(unique(days))))

# Just group samples by day
ourcolors = rainbow(length(unique(days)))

treatment = eset$characteristics_ch1.5

survival_time = eset$characteristics_ch1.7
survival_time = gsub('.*:', '', survival_time)
survival_time = as.character(survival_time)

# HIERARCHICAL CLUSTERING
# DENDROGRAM
d = as.dist(1 - cor(exprs(eset), method="spearman"))
sampleClustering = hclust(d)
sampleDendrogram = as.dendrogram(sampleClustering, hang=0.5)

batch = days
names(batch) = sampleNames(eset)
names(survival_time) = sampleNames(eset)
sampleDendrogram = dendrapply(sampleDendrogram, function(x, batch, labels) {
    ## for every node in the dendrogram if it is a leaf node
    if (is.leaf(x)) {
      print(x)
        attr(x, "nodePar") = list(lab.col = as.vector(ourcolors[batch_days[as.character(batch[attr(x, "label")]),]]))
        attr(x, "label") = as.vector(labels[attr(x, "label")])  ## label by survival_time
    }
    x
}, batch, survival_time)

plot(sampleDendrogram, main = "Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch))), fill = ourcolors)

# MULTIDIMENSIONAL SCALING (MDS)
cmd = cmdscale(as.dist(1 - cor(exprs(eset), method = "spearman")))
plot(cmd, type = "n")
text(cmd, survival_time, col = ourcolors[batch_days[as.character(batch),]], cex = 0.9)
legend("topleft", paste("Batch", unique(batch)), fill = ourcolors[batch_days[as.character(unique(batch)),]], inset = 0.01)

# maskNormalSamples = treatment == "Normal"
# normalSampleClustering = hclust(as.dist(1 - cor(exprs(eset[, maskNormalSamples]), 
#     method = "spearman")))
# normalSampleDendrogram = as.dendrogram(normalSampleClustering, hang = 0.1)
# normalBatch = batch[maskNormalSamples]
# normalOutcome = treatment[maskNormalSamples]
# normalSampleDendrogram = dendrapply(normalSampleDendrogram, function(x, batch, labels) {
#     ## for every node in the dendrogram if it is a leaf node
#     if (is.leaf(x)) {
#         attr(x, "nodePar") = list(lab.col = as.vector(batch[attr(x, "label")]))  ## color by batch
#         attr(x, "label") = as.vector(labels[attr(x, "label")])  ## label by outcome
#     }
#     x
# }, normalBatch, normaloutcome)

# QUANTIFY CONFOUNDING
s = fast.svd(t(scale(t(exprs(eset)), center=TRUE, scale=TRUE)))
plot(s$d^2/sum(s$d^2), type="b", lwd=2, las=1, xlab="Principal component", ylab="Proportion of variance")

par(mfrow=c(2, 3))
for (i in 1:6) {
  boxplot(split(s$v[, i], batch), main=sprintf("PC%d %.0f%%", i, 100 * s$d[i]^2/sum(s$d^2)))
}
