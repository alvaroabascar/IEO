library(affyPLM)
library(Biobase)
library(corpcor)
library(sva)

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
sort(days)

batch_days = data.frame(row.names=sort(unique(days)), batch=c(1:length(unique(days))))

# Just group samples by day
ourcolors = rainbow(length(unique(days)))

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
dim(cmd)
head(cmd)

# Ploting the results
plot(cmd, type = "n")
text(cmd, survival_time, col = ourcolors[batch_days[as.character(batch),]], cex = 0.9)
legend("topleft", paste("Batch", unique(batch)), fill = ourcolors[batch_days[as.character(unique(batch)),]], inset = 0.01)


# # Protocol to perform the Batch effect detection with only a concrete group, 
# # in the example we select a range of survival times smaller than 20 months.
# # This would be useful in case of needing to check in detail the batch effect 
# # in a determined set of samples. 
# 
# masktwentySamples = survival_time <= 20
# twentySampleClustering = hclust(as.dist(1 - cor(exprs(eset[, masktwentySamples]), 
#     method = "spearman")))
# twentySampleDendrogram = as.dendrogram(twentySampleClustering, hang = 0.1)
# twentyBatch = batch[masktwentySamples]
# twentyOutcome = survival_time[masktwentySamples]
# twentySampleDendrogram = dendrapply(twentySampleDendrogram, function(x, batch, labels) {
#     ## for every node in the dendrogram if it is a leaf node
#     if (is.leaf(x)) {
#     print(x)
#     attr(x, "nodePar") = list(lab.col = as.vector(ourcolors[batch_days[as.character(batch[attr(x, "label")]),]]))
#     attr(x, "label") = as.vector(labels[attr(x, "label")])  ## label by survival_time
#     }
#   x
# }, batch, survival_time)
# 
#     x
# }, twentyBatch, twentyOutcome)
# 
# plot(twentySampleDendrogram)
# legend("topright", paste("Batch", unique(batch)), fill = ourcolors[batch_days[as.character(unique(batch)),]], inset = 0.01)


# QUANTIFY CONFOUNDING by means of Principal Component Analysis

# Plot of proportion of variance explained perprincipal component
s = fast.svd(t(scale(t(exprs(eset)), center=TRUE, scale=TRUE)))
plot(s$d^2/sum(s$d^2), type="b", lwd=2, las=1, xlab="Principal component", ylab="Proportion of variance")


# Boxplots:
par(mfrow=c(2, 3))
for (i in 1:6) {
  boxplot(split(s$v[, i], batch), main=sprintf("PC%d %.0f%%", i, 100 * s$d[i]^2/sum(s$d^2)))
}


# SURROGATE VARIABLE ANALYSIS
st_bin <- ifelse(survival_time >= 20, 1, 0)
# Generating a full model
mod <- model.matrix(~st_bin, data = pData(eset))
head(mod)

# Generating the null model
mod0 <- model.matrix(~1, data = pData(eset))

# surrogate variables calling
sv <- sva(exprs(eset), mod, mod0)

# Plot the correlations
par(mfrow = c(2, 3))
for (i in 1:sv$n.sv) boxplot(sv$sv[, i] ~ batch, main = sprintf("SV %d", i), xlab = "Batch")

# Numbero of genes changing
pValues <- f.pvalue(exprs(eset), mod, mod0)
sum(p.adjust(pValues, method = "BH") < 0.05)

dim(eset)

# Changes after adjustment
modSv <- cbind(mod, sv$sv)
mod0Sv <- cbind(mod0, sv$sv)
pValuesSv <- f.pvalue(exprs(eset), modSv, mod0Sv)
sum(p.adjust(pValuesSv, method = "BH") < 0.05)