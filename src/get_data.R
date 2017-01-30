
# Table S1 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2850316/#s3title
download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2850316/bin/pntd.0000653.s002.xls',
              'data/Westenberger_et_al.xls')

# Table S1 from Pelle et al
download.file('https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-015-0133-7/MediaObjects/13073_2015_133_MOESM1_ESM.xlsx',
              'data/network_project/Pelle_et_al_TableS1.xlsx')

# GEO platform and series data from Boopathia et al.
library(GEOquery)
tmp <- getGEO("GSE55644", destdir = 'data/')