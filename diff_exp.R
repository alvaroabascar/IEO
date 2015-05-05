library(affyPLM)
library(Biobase)
library(genefilter)
library(limma)
library(sva)

MLL.B = readRDS('MLLB.rds')
eset = readRDS('eset.rds')
sampleNames(MLL.B) = letters[1:21]

eset = eset[, eset$characteristics_ch1.6 != "survival status: NA"]
eset$characteristics_ch1.1 = factor(eset$characteristics_ch1.1)
eset$characteristics_ch1.2 = factor(eset$characteristics_ch1.2)
eset$characteristics_ch1.3 = factor(eset$characteristics_ch1.3)
eset$characteristics_ch1.4 = factor(eset$characteristics_ch1.4)
eset$characteristics_ch1.5 = factor(eset$characteristics_ch1.5)
eset$characteristics_ch1.6 = factor(eset$characteristics_ch1.6)
eset$characteristics_ch1.7 = factor(eset$characteristics_ch1.7)
eset$characteristics_ch1.8 = factor(eset$characteristics_ch1.8)
eset$mgmt = eset$characteristics_ch1.9

threshold = 20
eset$survival_time = as.double(gsub('.*: ', '', eset$characteristics_ch1.7)) > threshold
eset$survival_time = factor(eset$survival_time)

# Load the two different survival groups & plot the log fold-changes between means:
# zeroExp <- rowMeans(exprs(eset[, eset$survival_time == "survival status: 0"]))
# oneExp <- rowMeans(exprs(eset[, eset$survival_time == "survival status: 1"]))

zeroExp <- rowMeans(exprs(eset[, eset$survival_time]))
oneExp <- rowMeans(exprs(eset[, eset$survival_time != TRUE]))

png('figures/diff_exp/fold_change_survival.png', width=10, height=6, units='in', res=700)
par(mfrow = c(1, 2))
plot(zeroExp, oneExp, xlab = "Zero", ylab = "One", pch = ".", cex = 4, las = 1)
plot((oneExp + zeroExp)/2, oneExp - zeroExp, pch = ".", cex = 4, las = 1)
dev.off()

allTests <- rowttests(eset, eset$survival_time)

# FDR adjustment
padjFDR <- p.adjust(allTests$p.value, method = "BH")
#sum(padjFDR < 0.01)

# After doing this we realize there are too many genes (so too many t-tests) and the
# adjustment returns 0 useful samples. 
# SOLUTION: remove genes without biological value. 

IQRs <- esApply(eset, 1, IQR)
png('figures/diff_exp/CDF_IQR.png', width=6, height=6, units='in', res=700)
plot.ecdf(IQRs, xlab = "IQR", main = "Empirical CDF of IQR values")
abline(v = quantile(IQRs, prob = 0.3), col = "red", lwd = 2)
dev.off()

maskHighVariability <- IQRs > quantile(IQRs, prob = 0.3)
eset_filtered <- eset[maskHighVariability, ]
#dim(eset_filtered)


# Again: 

fewerTests <- rowttests(eset_filtered, eset_filtered$survival_time)
padjBonf <- p.adjust(fewerTests$p.value, method = "bonferroni")
#sum(padjBonf < 0.01)

padjFDR <- p.adjust(fewerTests$p.value, method = "BH")
#sum(padjFDR < 0.01)


######
# MODERATED t-test using limma
######

design <- model.matrix(~survival_time, data = eset)

fit <- lmFit(eset, design)
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
# abline(h=-log10(max(tt[tt$adj.P.Val < 0.1, "P.Value"])), col=grey(0.5), lty=2)
dev.off()

# Adjust for covariates

design <- model.matrix(~survival_time + characteristics_ch1.4, data = eset)
fit <- lmFit(eset, design)
fit <- eBayes(fit)
# summary(decideTests(fit, p.value = 0.1))

# Surrogate Variables
mod <- model.matrix(~survival_time + characteristics_ch1.4, data = eset)
mod0 <- model.matrix(~characteristics_ch1.4, data = eset)

svaobj <- sva(exprs(eset), mod, mod0)

modSVs <- cbind(mod, svaobj$sv)

fit <- lmFit(eset, modSVs)
fit <- eBayes(fit)
ttadj <- topTable(fit, coef = 2, n = Inf)

png('figures/diff_exp/p_values_distr_after_sva.png', width=14, height=6, units='in', res=700)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(ttadj$P.Value, xlab = "Raw P-values", main = "")
hist(ttadj$P.Value, xlab = "Raw P-values", breaks = 1000, main = "")
dev.off()

# png('figures/diff_exp/volcano_after_sva.png', width=14, height=6, units='in', res=700)
# plot(ttadj$logfc, -log10(ttadj$P.Value), pch=".", cex=4, col=grey(0.75), cex.axis=1.2, las=1, cex.lab=1.5, xlab=expression(paste(log[2], " Fold change")), ylab=expression(paste(-log[10], " P-value")))
# points(ttadj[tt$adj.p.val < 0.1, "logFC"], -log10(ttadj[ttadj$adj.P.Val < 0.1, "P.Value"]), pch=".", cex=4, col="red")
# abline(h=-log10(max(ttadj[ttadj$adj.P.Val < 0.1, "P.Value"])), col=grey(0.5), lty=2)
# dev.off()
