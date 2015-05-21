options(error=utils::recover) 

library(affyPLM)
library(Biobase)
library(genefilter)
library(limma)
library(sva)

glioData = readRDS('rma_glioData.rds')
eset = readRDS('eset.rds')

nsamples = length(sampleNames(glioData))
sampleNames(glioData) = as.character(c(1:nsamples))

eset = eset[, eset$characteristics_ch1.6 != "survival status: NA"]
eset$characteristics_ch1.1 = factor(eset$characteristics_ch1.1) # disease status: re-recurrent GBM
eset$characteristics_ch1.2 = factor(eset$characteristics_ch1.2) # patient id
eset$characteristics_ch1.3 = factor(eset$characteristics_ch1.3) # age: dd.d
eset$characteristics_ch1.4 = factor(eset$characteristics_ch1.4) # gender: MALE/FEMALE
eset$characteristics_ch1.5 = factor(eset$characteristics_ch1.5) # treatment: radiotherapy or TMZ/radiotherapy
eset$characteristics_ch1.6 = factor(eset$characteristics_ch1.6) # survival status: 0 or 1
eset$characteristics_ch1.7 = factor(eset$characteristics_ch1.7) # survival time in months: 25.13
eset$characteristics_ch1.8 = factor(eset$characteristics_ch1.8) # mgmt status: M or U
eset$mgmt = eset$characteristics_ch1.8
eset$treatment = eset$characteristics_ch1.5

eset$survival_status = as.double(gsub('.*: ', '', eset$characteristics_ch1.6)) == 1
eset$survival_status = factor(eset$survival_status)

survival_threshold = 18
eset$survival_time = as.double(gsub('.*: ', '', eset$characteristics_ch1.7)) > survival_threshold
eset$survival_time = factor(eset$survival_time)

age_threshold = 50
eset$aged = as.double(gsub('.*: ', '', eset$characteristics_ch1.3)) > age_threshold
eset$aged = factor(eset$aged)

# Load the two different survival groups & plot the log fold-changes between means:
 zeroExp <- rowMeans(exprs(eset[, eset$survival_status]))
 oneExp <- rowMeans(exprs(eset[, eset$survival_status != TRUE]))

#zeroExp <- rowMeans(exprs(eset[, eset$survival_time]))
#oneExp <- rowMeans(exprs(eset[, eset$survival_time != TRUE]))

png('figures/diff_exp/fold_change_survival.png', width=10, height=6, units='in', res=700)
par(mfrow = c(1, 2))
plot(zeroExp, oneExp, xlab = "Zero", ylab = "One", pch = ".", cex = 4, las = 1)
plot((oneExp + zeroExp)/2, oneExp - zeroExp, pch = ".", cex = 4, las = 1)
dev.off()

allTests <- rowttests(eset, eset$survival_status)

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
#eset_filtered <- eset[maskHighVariability, ]
#dim(eset_filtered)

# Remove low variability genes by selecting 
# the most variable probe sets (standard deviation, >0.75),
# following the proposal from the paper
maskHighVariability_sd <- apply(exprs(eset), 1, sd)>0.75
#eset_filtered <- eset[maskHighVariability_sd, ]

mask = maskHighVariability & maskHighVariability_sd
eset_filtered <- eset[mask, ]

# Again:

fewerTests <- rowttests(eset_filtered, eset_filtered$survival_status)
padjBonf <- p.adjust(fewerTests$p.value, method = "bonferroni")
#sum(padjBonf < 0.01)

padjFDR <- p.adjust(fewerTests$p.value, method = "BH")
#sum(padjFDR < 0.01)

# To reduce even more the number of analysed probe sets 
# by using the nsFilter function. From the nsFilter help page:
# "Filtering features exhibiting little variation, or
#     a consistently low signal, across samples can be advantageous for
#     the subsequent data analysis (Bourgon et al.).  Furthermore, one
#     may decide that there is little value in considering features with
#     insufficient annotation."

library('hgu133plus2.db')
annotation(eset_filtered)<-'hgu133plus2.db'

filtered <- nsFilter(eset_filtered, require.entrez=TRUE,
         require.GOBP=FALSE, require.GOCC=FALSE,
         require.GOMF=FALSE, require.CytoBand=FALSE,
         remove.dupEntrez=TRUE, var.func=IQR,
         var.cutoff=0.5, var.filter=TRUE,
         filterByQuantile=TRUE, feature.exclude="^AFFX")
eset_filtered <- filtered$eset
fewerTests <- rowttests(eset_filtered, eset_filtered$survival_status)
padjFDR <- p.adjust(fewerTests$p.value, method = "BH")


######
# MODERATED t-test using limma
######

eset_filtered <- eset_filtered[, eset_filtered$treatment == "treatment: TMZ/radiotherapy"] ###########hack
eset_filtered$mgmtM = eset_filtered$mgmt == "mgmt status: M"

design <- model.matrix(~survival_status, data = eset_filtered)

fit <- lmFit(eset_filtered, design)
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

#design <- model.matrix(~survival_time + characteristics_ch1.5, data = eset_filtered)
#fit <- lmFit(eset_filtered, design)
#fit <- eBayes(fit)
# summary(decideTests(fit, p.value = 0.1))

# Surrogate Variables
# adjusted for age (> 50 years) and MGMT methylation status
mod <- model.matrix(~survival_status + mgmtM + aged, data = eset_filtered)
mod0 <- model.matrix(~ mgmtM + aged, data = eset_filtered)

svaobj <- sva(exprs(eset_filtered), mod, mod0)

modSVs <- cbind(mod, svaobj$sv)

fit <- lmFit(eset_filtered, modSVs)
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
