library(affyPLM)

dir.create(file.path('figures_all/qa'), showWarnings = FALSE)

MLL.B = readRDS('MLLB.rds')
names = as.character(c(1:length(sampleNames(MLL.B))))
sampleNames(MLL.B) = names

# Scanner plot analysis
png('figures_all/qa/scannerplot.png', width=4, height=4, units='in', res=500)
par(mfrow = c(3, 7), mar = c(1, 1, 3, 1))
image(MLL.B)
dev.off()

# Raw intensity boxplot
png('figures_all/qa/intensity_boxplot.png', width=8, height=4, units='in', res=500)
par(mar = c(4, 5, 0, 1))
boxplot(MLL.B)
dev.off()

# Density plot
png('figures_all/qa/density_plot.png')
par(mar = c(3, 3, 4, 1))
plotDensity.AffyBatch(MLL.B, lwd=2, col=1:length(names), lty=1:length(names))
legend('topright', names, col=1:length(names) lty=1:length(names) lwd=2, inset=0.1)
dev.off()

Pset = fitPLM(MLL.B)
for (i in 1:length(names)) {
  png(paste('figures_all/qa/PLM', i, '.png', sep=''), width=9, height=3, units='in', res=300)
  par(mfrow=c(1, 3))
  image(MLL.B[, i]) # raw intensities for sample A
  image(Pset, type="resids", which=i) # PLM resids for sample A
  image(Pset, type="weights", which=i) # PLM resids for sample A
  dev.off()
}

# Probe level modeling
png('figures_all/qa/NUSE.png', width=6, height=4, units='in', res=300)
NUSE(Pset)
dev.off()

# RLE
png('figures_all/qa/RLE.png', width=6, height=4, units='in', res=300)
RLE(Pset)
dev.off()

nuseDiag = NUSE(Pset, type='stats')
rleDiag = RLE(Pset, type='stats')
