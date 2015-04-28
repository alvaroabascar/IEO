library(affyPLM)
library(Biobase)

MLL.B = readRDS('MLLB.rds')
eset = readRDS('eset.rds')
sampleNames(MLL.B) = letters[1:21]

# RMA normalization (Robust Multi Array)

rmaEset = rma(MLL.B)

# MA
png('figures/normalization/MA.png', height=8, width=14, units='in', res=300)
par(mfrow=c(3, 7), mar=c(1, 1, 1, 1))
MAplot(rmaEset[-grep('^AFFX', featureNames(rmaEset))], plot.method='smoothScatter', cex=0.75, ref.title='')
dev.off()

