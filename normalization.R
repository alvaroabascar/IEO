library(affyPLM)
library(Biobase)

dir.create(file.path('figures_all/normalization'), showWarnings = FALSE)

MLL.B = readRDS('MLLB.rds')
eset = readRDS('eset.rds')
sampleNames(MLL.B) = as.character(c(1:length(sampleNames(MLL.B))))

# RMA normalization (Robust Multi Array)

rmaEset = rma(MLL.B)

# MA
png('figures_all/normalization/MA.png', height=8, width=14, units='in', res=300)
par(mfrow=c(3, 7), mar=c(1, 1, 1, 1))
MAplot(rmaEset[-grep('^AFFX', featureNames(rmaEset))], plot.method='smoothScatter', cex=0.75, ref.title='')
dev.off()

