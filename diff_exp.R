library(affyPLM)
library(Biobase)
library(genefilter)
library(limma)

MLL.B = readRDS('MLLB.rds')
eset = readRDS('eset.rds')
sampleNames(MLL.B) = letters[1:21]

# Load the two different survival groups & plot the log fold-changes between means:
zeroExp <- rowMeans(exprs(eset[, eset$characteristics_ch1.6 == "survival status: 0"]))
oneExp <- rowMeans(exprs(eset[, eset$characteristics_ch1.6 == "survival status: 1"]))

png('figures/diff_exp/fold_change_survival.png', width=10, height=6, units='in', res=700)
par(mfrow = c(1, 2))
plot(zeroExp, oneExp, xlab = "Zero", ylab = "One", pch = ".", cex = 4, las = 1)
plot((oneExp + zeroExp)/2, oneExp - zeroExp, pch = ".", cex = 4, las = 1)
dev.off()


# Multiple testing (we need two classes so the first step is just to remove the NA class)
eset_glio <- eset[,eset$characteristics_ch1.6 != "survival status: NA"]
eset_glio$characteristics_ch1.6 = factor(eset_glio$characteristics_ch1.6)
allTests <- rowttests(eset_glio, eset_glio$characteristics_ch1.6)

# FDR adjustment
padjFDR <- p.adjust(allTests$p.value, method = "BH")
#sum(padjFDR < 0.01)

# After doing this we realize there are too many genes (so too many t-tests) and the
# adjustment returns 0 useful samples. 
# SOLUTION: remove genes without biological value. 

IQRs <- esApply(eset_glio, 1, IQR)
png('figures/diff_exp/CDF_IQR.png', width=6, height=6, units='in', res=700)
plot.ecdf(IQRs, xlab = "IQR", main = "Empirical CDF of IQR values")
abline(v = quantile(IQRs, prob = 0.3), col = "red", lwd = 2)
dev.off()

maskHighVariability <- IQRs > quantile(IQRs, prob = 0.3)
eset_glio_filtered <- eset_glio[maskHighVariability, ]
#dim(eset_glio_filtered)


# Again: 

fewerTests <- rowttests(eset_glio_filtered, eset_glio_filtered$characteristics_ch1.6)
padjBonf <- p.adjust(fewerTests$p.value, method = "bonferroni")
#sum(padjBonf < 0.01)

padjFDR <- p.adjust(fewerTests$p.value, method = "BH")
#sum(padjFDR < 0.01)


######
# MODERATED t-test using limma
######

eset_glio$characteristics_ch1.6 <- relevel(eset_glio$characteristics_ch1.6, ref = "survival status: 0")
design <- model.matrix(~characteristics_ch1.6, data = eset_glio)

fit <- lmFit(eset_glio, design)
fit <- eBayes(fit)

res <- decideTests(fit)
# summary(res)

res <- decideTests(fit, p.value = 0.1)
# summary(res)

tt <- topTable(fit, coef = 2, n = Inf)

png('figures/diff_exp/p_values_distr.png', width=14, height=6, units='in', res=700)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt$P.Value, xlab = "Raw P-values", main = "")
hist(tt$P.Value, xlab = "Raw P-values", breaks = 1000, main = "")
dev.off()

png('figures/diff_exp/volcano.png', width=14, height=6, units='in', res=700)
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75), cex.axis=1.2, las=1, cex.lab=1.5, xlab=expression(paste(log[2], " Fold change")), ylab=expression(paste(-log[10], " P-value")))
points(tt[tt$adj.P.Val < 0.1, "logFC"], -log10(tt[tt$adj.P.Val < 0.1, "P.Value"]), pch=".", cex=4, col="red")
abline(h=-log10(max(tt[tt$adj.P.Val < 0.1, "P.Value"])), col=grey(0.5), lty=2)
dev.off()

# SLIDE 48
