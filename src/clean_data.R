rm(list=ls())
library(xlsx)
library(magrittr)
library(dplyr)
library(siyuan)

# clean up data for network project-------------------------------
# Table for asex-vs-gam/stage annotation. Gam clusters are visually inspected
# from Pelle et al., Fig. 6B. Asex clusters are the rest, and their
# sequestration info obtained from Pelle et al., Table S1.
df.clst.TableS1 <- read.xlsx2('data/network_project/Pelle_et_al_TableS1.xlsx',
                         sheetIndex = 1,
                         # row.names = 2,
                         check.names = F,
                         stringsAsFactors = F) # Pelle et al. Table S1
df.clst.TableS1 <- df.clst.TableS1 %>% 
  mutate(mean_peak_time = `Mean Peak Time (hr) per cluster based on Bozdech et al time course (2003)(Figure 2A)` %>% as.numeric, 
         cluster = `Cluster index` %>% as.numeric)
# sanity check, as clusters should be ordered by peak time
plot( df.clst.TableS1$cluster, df.clst.TableS1$mean_peak_time )

df.clst.Figure6B <- read.table('data/network_project/pelleetall_fig6_clsts.txt',
                               header = F,
                               stringsAsFactors = F,
                               sep = '\t')

df.clst <- merge(df.clst.TableS1, df.clst.Figure6B, 
                 by.x='cluster', by.y=1, all = T)
df.clst <- df.clst %>% 
  mutate(stage=ifelse(is.na(V2), 'asex', V2),
         sequestration=ifelse(stage == 'asex', 
                              ifelse(mean_peak_time<=22, 
                                     'circulation',
                                     'sequestration'),
                              V3))
df.clst %>% 
  select(cluster, mean_peak_time, stage, sequestration) %>% 
  write.table(file='data/network_project/cluster_table_clean.txt',
              row.names=F,
              quote=F,
              sep='\t')

# clean up per-clsuter gene names (into old version)
l.clst <- read.table('data/network_project/cluster_genes.txt',
                     header = T, stringsAsFactors = F, 
                     sep = '\t')[, 2] %>% 
  lapply(function(i.id) {
    i.id %>% strsplit(split = ', ') %>% '[['(1) %>% 
      gsub('W$', 'w', ., perl = T) %>% gsub('C$', 'c', ., perl = T) %>% 
      return
  })

dir.create('data/network_project/cluster_genes_old_id_clean/', showWarnings = F)
(1:length(l.clst)) %>% lapply(function(i.clst) {
  i.id <- l.clst[[i.clst]]
  data.frame(i.id) %>% 
    write.table(paste0('data/network_project/cluster_genes_old_id_clean/clst_',
                       i.clst, '.txt'),
                col.names = F,
                row.names = F,
                quote = F
                )
})

# Wetenberger expression data
df.exprs.Westenberger <- read.xlsx2('data/Westenberger_et_al.xls',
                                    sheetIndex = 1,
                                    startRow = 3,
                                    # row.names = 1,
                                    check.names = F,
                                    stringsAsFactors = F)
# It seems gene "PVX_123355" has two identical rows
# df.exprs.Westenberger[df.exprs.Westenberger[, 1]=='PVX_123355', ]
df.exprs.Westenberger <- df.exprs.Westenberger[!duplicated(df.exprs.Westenberger[, 1]), ]

row.names(df.exprs.Westenberger) <- df.exprs.Westenberger[, 1]
df.exprs.Westenberger.clean <- df.exprs.Westenberger[, -(1:19)]
mat.exprs.Westenberger.clean <- df.exprs.Westenberger.clean %>% data.matrix

# data is already normalized
dir.create2('result/clean_data/') %>% 
  paste0('boxplot_Westenberger.pdf') %>% 
  pdf(width=8, height=6)
par(mar=c(14,3,1,1))
boxplot(log(mat.exprs.Westenberger.clean), las=2) 
dev.off()

# cleaner column names
colnames(mat.exprs.Westenberger.clean) <- colnames(mat.exprs.Westenberger.clean) %>% 
  gsub('^.*\\: ', '', ., perl=T)
# shifting columns to move the asexual blood samples first
mat.exprs.Westenberger.clean <- mat.exprs.Westenberger.clean[, c(3:14, 1:2)]
write.csv(mat.exprs.Westenberger.clean, file = 'data/Westenberger_et_al_clean.csv')
