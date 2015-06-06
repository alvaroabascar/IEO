library(affy)
library(Biobase)
library(GEOquery)

# Cache the data to avoid having to load it again next time from the .txt
if (file.exists('full_eset.rds')) {
  eset = readRDS('full_eset.rds')
} else {
  eset = getGEO(filename='data/GSE7696_series_matrix.txt.gz')
  saveRDS(eset, 'full_eset.rds')
}

# extract sample names
samples = colnames(eset)

# separate the samples by treatment group
#no_treated = samples[eset$characteristics_ch1.5 == "treatment: NA"]
tmz_rt = samples[eset$characteristics_ch1.5 == "treatment: TMZ/radiotherapy"]
#rt = samples[eset$characteristics_ch1.5 == "treatment: radiotherapy"]

# take 21 samples, keeping the proportions of the original dataset
set.seed(123456)
#no_treated = sample(no_treated, size=1)
#tmz_rt = sample(tmz_rt, size=13)
#rt = sample(rt, size=7)

#subsample_names = c(no_treated, tmz_rt, rt)
subsample_names = c(tmz_rt)
# turn sample names into filenames
subsample_filenames = paste(subsample_names, '.cel.gz', sep='')

# read only the necessary the expression data
glioData = ReadAffy(celfile.path = "./data/raw", filenames=subsample_filenames)
subsampled_eset = eset[, subsample_names]

# save the subsampled data for posterior use
saveRDS(subsampled_eset, 'eset.rds')
saveRDS(glioData, 'glioData.rds')
