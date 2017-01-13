rm(list=ls())
library(ggplot2)
library(tidyr)
source('src/source.R')
source('src/generate_pv_annotation.R') # needed for PV annotations

# log transform data first
df.exprs.Westenberger <- read.csv('data/Westenberger_et_al_clean.csv', 
                                  header=T, row.names = 1,
                                  check.names = F) %>% log
# #genes available in Westenberger et al.
l.anno.pv %>% lapply(intersect, rownames(df.exprs.Westenberger)) %>% 
  lapply(unique) %>% sapply(length)

df.exprs.avg <- l.anno.pv %>% sapply(function(i.group){
  return( df.exprs.Westenberger[rownames(df.exprs.Westenberger) %in% i.group,] %>% 
    apply(2, mean) 
  )
}) %>% data.frame
df.exprs.avg$sample_type <- 'blood_asexual'
df.exprs.avg$sample_type[11] <- 'blood_gamete/zygote'
df.exprs.avg$sample_type[12] <- 'blood_ookinete'
df.exprs.avg$sample_type[13:14] <- 'saliva_sporozoite'
df.exprs.avg$sample_type <- df.exprs.avg$sample_type %>%
  factor(c('blood_asexual', 'blood_gamete/zygote', 'blood_ookinete', 'saliva_sporozoite'))

df.exprs.avg.long <- gather(df.exprs.avg, 
                            key = stage,
                            value = mean_expression,
                            asex:rings)
df.exprs.avg.long$asex_or_gam <- 'gam'
df.exprs.avg.long$asex_or_gam[df.exprs.avg.long$stage %in% 
                                  c('asex', 'asexCirc', 'asexSeq')] <- 'asex'
df.exprs.avg.long$stage <- df.exprs.avg.long$stage %>% factor(
  levels=c('asex', 'asexCirc', 'asexSeq', 'rings', 'IG', 'IGYG', 'MG', 'commit'))

pdf('result/boxplot_averageByStage_Westenberger.pdf', width=8, height=6)
ggplot(subset(df.exprs.avg.long, !(stage %in% c('asex', 'IGYG'))),
       aes(x=stage, y=mean_expression)) +
  geom_boxplot(aes(color=sample_type), alpha=0.5, outlier.shape = NULL) +
  geom_point(aes(color=sample_type), position = position_jitterdodge()) +
  facet_grid(.~asex_or_gam, scale='free_x', space='free_x') +
  coord_fixed(ratio=1) +
  theme_bw()
dev.off()

#####################################################################
# by cluster analysis?
#####################################################################
library(gplots)

df.clst <- read.table('data/network_project/pelleetall_fig6_clsts.txt', 
                      header = F,
                      stringsAsFactors = F,
                      row.names = 1,
                      sep = '\t')
colnames(df.clst) <- c('stage', 'circulation')

my.breaks <- seq(min(df.exprs.Westenberger), max(df.exprs.Westenberger),
                    length.out = 300)
my.pallete <- greenred(299)
mat.avg.clst.exprs <- rownames(df.clst) %>% sapply(function(i.clst) {
  i.dir.output <- df.clst[i.clst, ]$stage %>% 
    paste0('result/Westenberger/by_clst_heatmap/', ., '/')
  dir.create(i.dir.output, recursive = T, showWarnings = F)
  
  i.gene <- i.clst %>% 
    paste0('data/network_project/cluster_genes_old_id_clean/clst_', 
           ., '.txt') %>% 
    read.table(header = F, stringsAsFactors = F) %>% '['(, 1) %>% 
    pf.to.pv(old.ID = T, print.n.match = F) %>% 
    intersect(row.names(df.exprs.Westenberger))
  i.mat.exprs <- as.matrix(df.exprs.Westenberger[i.gene, ])
  
  if(length(i.gene) > 1) {
    pdf(paste0(i.dir.output, 'clst_', i.clst, '.pdf'), 
        width = 8, 
        height = 3/ncol(i.mat.exprs)*nrow(i.mat.exprs) + 3)
    heatmap.2(i.mat.exprs, Rowv = NA, Colv = NA, symm = F, 
              dendrogram = 'none', scale = 'none', 
              trace = 'none', symkey = F,
              breaks = my.breaks, col=my.pallete)
    dev.off()
  } 
  else write.table(data.frame(i.gene), 
                   paste0(i.dir.output, 'clst_', i.clst, 'txt'),
                   col.names = F, row.names = F)
  
  i.mat.exprs %>% apply(2, mean) %>% return
})

pdf('result/Westenberger/heatmap_clst_avg.pdf',
    width = 4/nrow(mat.avg.clst.exprs)*ncol(mat.avg.clst.exprs), 
    height = 8)
heatmap.2(mat.avg.clst.exprs, Rowv = NA, Colv = NA, symm = F, 
          dendrogram = 'none', scale = 'none', 
          trace = 'none', symkey = F,
          col = greenred(300))
dev.off()