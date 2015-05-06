library(affyPLM)
library(Biobase)

glioData = readRDS('glioData.rds')
eset = readRDS('eset.rds')
sampleNames(glioData) = letters[1:21]

# RMA normalization (Robust Multi Array)

rmaEset = rma(glioData)

# MA
png('figures/normalization/MA.png', height=8, width=14, units='in', res=300)
par(mfrow=c(3, 7), mar=c(1, 1, 1, 1))
MAplot(rmaEset[-grep('^AFFX', featureNames(rmaEset))], plot.method='smoothScatter', cex=0.75, ref.title='')
dev.off()

