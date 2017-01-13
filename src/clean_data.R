rm(list=ls())
library(xlsx)
library(magrittr)

# clean up script for network project
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
l_clstGenes <- lapply( read.table( 'data/network_project/cluster_genes.txt', header=T, stringsAsFactors = F, sep='\t' )[, 2], 
                       function( genes ) {
                         genes <- strsplit( genes, split=', ' )[[1]] 
                         return( gsub( 'W$', 'w', gsub( 'C$', 'c', genes, perl=T ) ) )
                       } )

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
pdf('result/clean_data/boxplot_Westenberger.pdf', width=8, height=6)
par(mar=c(14,3,1,1))
boxplot(log(mat.exprs.Westenberger.clean), las=2) 
dev.off()

# cleaner column names
colnames(mat.exprs.Westenberger.clean) <- colnames(mat.exprs.Westenberger.clean) %>% 
  gsub('^.*\\: ', '', ., perl=T)
# shifting columns to move the asexual blood samples first
mat.exprs.Westenberger.clean <- mat.exprs.Westenberger.clean[, c(3:14, 1:2)]
write.csv(mat.exprs.Westenberger.clean, file = 'data/Westenberger_et_al_clean.csv')
