library(GSEABase)
library(affyPLM)
library(Biobase)
library(genefilter)
library(hgu133plus2.db)
library(GSVA)
library(sva)
library(limma)
library(org.Hs.eg.db)

glioData = readRDS('rma_glioData.rds')
eset = readRDS('eset.rds')
annotation(eset)<-'hgu133plus2.db'

####
age_threshold = 50
survival_threshold = 20

eset$mgmtM = eset$characteristics_ch1.8 == "mgmt status: M"
eset$mgmtM = factor(eset$mgmtM)
eset$survival = as.double(gsub('.*: ', '', eset$characteristics_ch1.7)) > survival_threshold
eset$survival = factor(eset$survival)
eset$aged = as.double(gsub('.*: ', '', eset$characteristics_ch1.3)) > age_threshold
eset$aged = factor(eset$aged)
eset$recurrent = factor(grepl("recurrent", eset$characteristics_ch1.1))
###

# Analyse for recurrence
mod <- model.matrix(~recurrence + mgmtM + aged, data=eset)
mod0 <- model.matrix(~ aged + mgmtM, data=eset)
svaobj <- sva(exprs(eset), mod, mod0)
modSVs <- cbind(mod, svaobj$sv)

fit <- lmFit(eset, modSVs)
fit <- eBayes(fit)
ttadj <- topTable(fit, coef = 2, n = Inf, adjust.method="bonferroni")

#png('figures/diff_exp/p_values_distr_after_sva.png', width=14, height=6, units='in', res=700)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(ttadj$P.Value, xlab = "Raw P-values", main = "")
hist(ttadj$P.Value, xlab = "Raw P-values", breaks = 1000, main = "")
DEgenes = rownames(tt)

## GO ONTOLOGY ##
microarrayGO = select(hgu133plus2.db, columns=c("SYMBOL", "GO"), keys=featureNames(glioData))
dim(microarrayGO)
length(microarrayGO$GO)
length(unique(microarrayGO$GO))
head(microarrayGO)
