library(GSEABase)
library(affyPLM)
library(Biobase)
library(genefilter)
library(hgu133plus2.db)
library(GSVA)
library(sva)
library(limma)

glioData = readRDS('rma_glioData.rds')
eset = readRDS('eset.rds')
annotation(eset)<-'hgu133plus2.db'

clusters = read.csv('scrap_clusters/clusters.csv', header=TRUE)

all_genesets = c()
for (cluster_name in levels(clusters$cluster)) {
    gene_ids = as.character(clusters[clusters$cluster == cluster_name,]$geneId)
    cluster = GeneSet(AnnotationIdentifier("hgu133plus2.db"), geneIds=gene_ids, setName=cluster_name)
    all_genesets = c(all_genesets, cluster)
}

cluster_collection = GeneSetCollection(all_genesets)
cluster_collection
length(cluster_collection)
head(names(cluster_collection))

cluster_collection = mapIdentifiers(cluster_collection, AnnoOrEntrezIdentifier(annotation(eset)))
Im = incidence(cluster_collection)

GSeset = gsva(eset, cluster_collection, min.sz=10, max.sz=500, verbose=FALSE)$es.obs
class(GSeset)
dim(GSeset)

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

# mod <- model.matrix(~survival + mgmtM + aged, data=eset)
# mod0 <- model.matrix(~ mgmtM + aged, data=eset)
mod <- model.matrix(~recurrent + mgmtM + aged, data=eset)
mod0 <- model.matrix(~ aged + mgmtM, data=eset)

svaobj <- sva(exprs(eset), mod, mod0)

modSVs <- cbind(mod, svaobj$sv)

fit <- lmFit(GSeset, modSVs)
fit <- eBayes(fit)
ttadj <- topTable(fit, coef = 2, n = Inf, adjust.method="bonferroni")

png('figures/diff_exp/p_values_distr_after_sva.png', width=14, height=6, units='in', res=700)
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(ttadj$P.Value, xlab = "Raw P-values", main = "")
hist(ttadj$P.Value, xlab = "Raw P-values", breaks = 1000, main = "")

###### COX
pvalues = c()
nm = c()
for (rn in rownames(GSeset)) {
  g_eset = GSeset[rn,]
  survival_time = as.double(gsub('.*: ', '', g_eset$characteristics_ch1.7))
  survival_status = as.double(gsub('.*: ', '', g_eset$characteristics_ch1.6)) == 1
  mgmt = g_eset$characteristics_ch1.8 == "mgmt status: M"
  age_threshold = 50
  aged = as.double(gsub('.*: ', '', g_eset$characteristics_ch1.3)) > age_threshold
  exprs_vals=as.vector(exprs(g_eset))
  cox_df = data.frame(exprs_vals=exprs_vals, mgmt=mgmt, aged=aged, survival_time=survival_time, survival_status=survival_status)
  cox_result=coxph(Surv(survival_time, survival_status)~exprs_vals+mgmt+aged, cox_df)
  #print(cox_result)
  #pvalues = c(pvalues, summary(cox_result)$waldtest[3])
  pvalues = c(pvalues, summary(cox_result)$coefficients[1,5])
  nm = c(nm, rn)
}

stable_clusters = c("G2", "G7", "G9", "G11", "G12", "G13", "G14", "G16", "G18", "G12", "G21", "G22", "G23", "G24", "G25", "G27", "G28", "G29")
