library(affyPLM)

glioData = readRDS('glioData.rds')
nsamples = length(sampleNames(glioData))
sampleNames(glioData) = as.character(c(1:nsamples))

# Scanner plot analysis
png('figures/qa/scannerplot.png', width=4, height=4, units='in', res=500)
par(mfrow = c(3, round(nsamples/3)), mar = c(1, 1, 3, 1))
image(glioData)
dev.off()

# Raw intensity boxplot
png('figures/qa/intensity_boxplot.png', width=8, height=4, units='in', res=500)
par(mar = c(4, 5, 0, 1))
boxplot(glioData)
dev.off()

# Density plot
png('figures/qa/density_plot.png')
par(mar = c(3, 3, 4, 1))
plotDensity.AffyBatch(glioData, lwd=2, col=1:21, lty=1:21)
legend('topright', LETTERS[1:21], col=1:21, lty=1:21, lwd=2, inset=0.1)
dev.off()

Pset = fitPLM(glioData)
for (i in 1:21) {
  png(paste('figures/qa/PLM', i, '.png', sep=''), width=9, height=3, units='in', res=300)
  par(mfrow=c(1, 3))
  image(glioData[, i]) # raw intensities for sample A
  image(Pset, type="resids", which=i) # PLM resids for sample A
  image(Pset, type="weights", which=i) # PLM resids for sample A
  dev.off()
}

# Probe level modeling
png('figures/qa/NUSE.png', width=6, height=4, units='in', res=300)
NUSE(Pset)
dev.off()

# RLE
png('figures/qa/RLE.png', width=6, height=4, units='in', res=300)
RLE(Pset)
dev.off()

nuseDiag = NUSE(Pset, type='stats')
rleDiag = RLE(Pset, type='stats')
