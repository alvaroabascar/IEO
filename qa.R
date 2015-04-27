library(affyPLM)

MLL.B = readRDS('MLLB.rds')
# Scanner plot analysis
sampleNames(MLL.B) <- letters[1:21]
png('figures/qa/scannerplot.png', width=4, height=4, units='in', res=500)
par(mfrow = c(3, 7), mar = c(1, 1, 3, 1))
image(MLL.B)
dev.off()

# Raw intensity boxplot
png('figures/qa/intensity_boxplot.png', width=8, height=4, units='in', res=500)
par(mar = c(4, 5, 0, 1))
boxplot(MLL.B)
dev.off()

# Density plot
png('figures/qa/density_plot.png')
par(mar = c(3, 3, 4, 1))
plotDensity.AffyBatch(MLL.B, lwd=2, col=1:21, lty=1:21)
legend('topright', LETTERS[1:21], col=1:21, lty=1:21, lwd=2, inset=0.1)
dev.off()

# # Probe level modeling
# Pset = fitPLM(MLL.B)
# png('figures/qa/NUSE.png', width=4, height=4, units='in', res=500)
# NUSE(Pset)
# dev.off()
#
