library(affy)
library(Biobase)
library(GEOquery)

if (file.exists('full_eset.rds')) {
  eset = readRDS('full_eset.rds')
} else {
  eset = getGEO(filename='GSE7696_series_matrix.txt.gz')
  saveRDS(eset, 'full_eset.rds')
}

pData = pData(eset)

# extract sample names
samples = colnames(eset)

# take total 21 samples, with the same proportion across treatments as the
# original dataset
no_treated = samples[eset$characteristics_ch1.5 == "treatment: NA"]
tmz_rt = samples[eset$characteristics_ch1.5 == "treatment: TMZ/radiotherapy"]
rt = samples[eset$characteristics_ch1.5 == "treatment: radiotherapy"]
no_treated = sample(no_treated, size=1)
tmz_rt = sample(tmz_rt, size=13)
rt = sample(rt, size=7)

filtered_samples = c(no_treated, tmz_rt, rt)
filtered_samples_fullnames = paste(filtered_samples, '.cel.gz', sep='')

MLL.B = ReadAffy(celfile.path = "./data", filenames=filtered_samples_fullnames)
filteredeset = eset[, filtered_samples]

saveRDS(filteredeset, 'eset.rds')
saveRDS(MLL.B, 'MLLB.rds')

