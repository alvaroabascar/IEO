library(affyPLM)
library(Biobase)

if (!file.exists('figures/normalization')) {
  dir.create('figures/normalization')
}

glioData = readRDS('glioData.rds')
eset = readRDS('eset.rds')

nsamples = length(sampleNames(glioData))
sampleNames(glioData) = as.character(c(1:nsamples))

# RMA normalization (Robust Multi Array)

rma_glioData = rma(glioData)

# MA
png('figures/normalization/MA.png', height=8, width=14, units='in', res=300)
par(mfrow=c(3, 7), mar=c(1, 1, 1, 1))
MAplot(rma_glioData[-grep('^AFFX', featureNames(rma_glioData))], plot.method='smoothScatter', cex=0.75, ref.title='')
dev.off()

saveRDS(rma_glioData, 'rma_glioData.rds')

