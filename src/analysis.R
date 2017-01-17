rm(list=ls())
library(ggplot2)
library(tidyr)
library(dplyr)
source('src/source.R')

# read in data------------------------------------------------------
# expression data
df.exprs.Westenberger <- read.csv('data/Westenberger_et_al_clean.csv', 
                                  header=T, row.names = 1,
                                  check.names = F) %>% log # log-transformation is needed
# disard the non-blood samples, average the techinical replicates (as indicated
# in the manuscript)
df.exprs.Westenberger <- data.frame(df.exprs.Westenberger[, 1:6],
                                    CM12 = apply(df.exprs.Westenberger[, 7:8], 1, mean),
                                    CM13 = apply(df.exprs.Westenberger[, 9:10], 1, mean))

# read in (clustered) PF genes from Pelle et al. and conver to PV
l.clst.pv <- 1:284 %>% lapply(function(i.clst) {
  i.clst %>% 
    paste0('data/network_project/cluster_genes_old_id_clean/clst_', ., '.txt') %>% 
    read.table(header = F, stringsAsFactors = F) %>% '['(, 1) %>% 
    pf.to.pv(old.ID = T, print.n.match = F) %>% 
    return
})
# number of genes with cluster annotation available in Westenberger et al.
# ~60% has annotation
l.clst.pv %>% unlist %>% intersect(rownames(df.exprs.Westenberger)) %>% 
  length
l.clst.pv <- l.clst.pv %>% 
  lapply(intersect, rownames(df.exprs.Westenberger))


# by cluster----------------------------------------------------------
library(gplots)
df.clst <- read.table('data/network_project/cluster_table_clean.txt',
                      header=T,
                      sep='\t')
df.clst <- df.clst %>% 
  mutate(stage = factor(stage, levels=c('asex', 'rings', 'IG', 'MG', 'commit')),
         sequestration = factor(sequestration, levels=c('circulation', 'sequestration', '')))

# per-cluster mean expression matrix
my.breaks <- seq(min(df.exprs.Westenberger), max(df.exprs.Westenberger),
                 length.out = 300)
my.pallete <- greenred(299) # unified color scheme for per-cluster heatmaps
mat.avg.clst.exprs <- df.clst$cluster %>% sapply(function(i.clst) {
  i.gene <- l.clst.pv[[i.clst]]
  i.mat.exprs <- as.matrix(df.exprs.Westenberger[i.gene, ])
  
  # if gam cluster, output per-cluster heatmap
  if(df.clst[i.clst, ]$stage != 'asex') {
    i.dir.output <- df.clst[i.clst, ]$stage %>% 
      paste0('result/Westenberger/by_clst_heatmap/', ., '/') %>% 
      dir.create2
    
    if(length(i.gene) > 1 ) {
      i.dir.output %>% paste0('clst_', i.clst, '.pdf') %>% 
        pdf(width = 8, 
            height = 3/ncol(i.mat.exprs)*nrow(i.mat.exprs) + 5)
      heatmap.2(i.mat.exprs, Rowv = NA, Colv = NA, symm = F, 
                dendrogram = 'none', scale = 'none', 
                trace = 'none', symkey = F,
                breaks = my.breaks, col=my.pallete)
      dev.off()
    } 
    else write.table(i.mat.exprs, 
                     paste0(i.dir.output, 'clst_', i.clst, 'txt'),
                     col.names = F, row.names = F)
  }
  
  i.mat.exprs %>% apply(2, mean) %>% return
}) # per-cluster mean expression matrix

# gametocyte clusters heatmap
clst.gam <- read.table('data/network_project/pelleetall_fig6_clsts.txt',
                       header = F,
                       sep = '\t')[, 1]
mat.gam <- mat.avg.clst.exprs[, match(clst.gam, df.clst$cluster)]
mat.gam.avg <- mat.gam %>% apply(2, mean) %>% 
  matrix(nrow = 1, dimnames = list('mean', colnames(mat.gam))) # additional row of mean expression across patients
colnames(mat.gam) <- clst.gam
my.breaks <- c(seq(min(mat.gam), median(mat.gam), length.out = 500),
               seq(median(mat.gam), quantile(mat.gam, 0.9), length.out = 401)[-1],
               seq(quantile(mat.gam, 0.9), max(mat.gam), length.out = 101)[-1])
my.pallete <- greenred(999) # unified color scheme for per-cluster heatmaps
pdf('result/Westenberger/heatmap_gam_clst_avg.pdf',
    width = 4/nrow(mat.gam)*ncol(mat.gam), 
    height = 8)
heatmap.2(rbind(mat.gam.avg,
                matrix(NA, nrow = 1, ncol = ncol(mat.gam.avg),
                       dimnames = list('', colnames(mat.gam.avg))),
                mat.gam), 
          Rowv = NA, Colv = NA, symm = F, 
          dendrogram = 'none', scale = 'none', 
          trace = 'none', symkey = F,
          breaks = my.breaks,
          col = my.pallete)
dev.off()

# asexual clusters heatmap
clst.asex <- setdiff(df.clst$cluster, clst.gam)
mat.asex <- mat.avg.clst.exprs[, match(clst.asex, df.clst$cluster)]
colnames(mat.asex) <- clst.asex
mat.asex.avg <- mat.asex %>% apply(2, mean) %>% 
  matrix(nrow = 1, dimnames = list('mean', colnames(mat.asex))) # additional row of mean expression across patients
my.breaks <- c(seq(min(mat.asex, na.rm = T), median(mat.asex, na.rm = T), length.out = 500),
               seq(median(mat.asex, na.rm = T), quantile(mat.asex, 0.9, na.rm = T), length.out = 401)[-1],
               seq(quantile(mat.asex, 0.9, na.rm = T), max(mat.asex, na.rm = T), length.out = 101)[-1])
my.pallete <- greenred(999) # unified color scheme for per-cluster heatmaps
pdf('result/Westenberger/heatmap_asex_clst_avg.pdf',
    width = 1/nrow(mat.asex)*ncol(mat.asex), 
    height = 8)
heatmap.2(rbind(mat.asex.avg,
                matrix(NA, nrow = 1, ncol = ncol(mat.asex.avg),
                       dimnames = list('', colnames(mat.asex.avg))),
                mat.asex), 
          Rowv = NA, Colv = NA, symm = F, 
          dendrogram = 'none', scale = 'none', 
          trace = 'none', symkey = F,
          breaks = my.breaks,
          col = my.pallete)
dev.off()


# identify breakpoint in asex data ----------------------------------------
library(fpc)
mat.asex.complete <- mat.asex %>% 
  apply(2, function(x) !any(is.na(x))) %>% 
  '['(mat.asex, ,.)
dist.asex <-  mat.asex.complete %>% t %>% dist
mat.summary <- 2:ncol(mat.asex.complete) %>% sapply(function(i.col) {
  c(rep(1, i.col - 1), rep(2, ncol(mat.asex.complete) - i.col + 1)) %>% 
    cluster.stats(d = dist.asex, clustering = ., 
                  silhouette = TRUE,
                  wgap=FALSE, 
                  sepindex=FALSE,
                  sepwithnoise=FALSE) %>% 
    '['(c('within.cluster.ss', 'avg.silwidth', 'ch')) %>% 
    unlist %>% 
    return
})
df.summary <- data.frame(t(mat.summary), 
                         cluster = factor(colnames(mat.asex.complete)[-1],
                                          levels = colnames(mat.asex.complete)[-1])
                         )
pdf('result/Westenberger/asex_break_summary.pdf', width = 12, height = 10)
df.summary %>% 
  gather(key = statistic, value = value, within.cluster.ss:ch) %>% 
  ggplot(aes(x = as.numeric(cluster), y = value)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name = 'cluster', 
                     breaks = as.numeric(df.summary$cluster), 
                     labels = levels(df.summary$cluster)) +
  facet_wrap(~statistic, ncol = 1, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=4),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

dist.asex.avg <-  mat.asex.avg %>% 
  apply(2, function(x) !any(is.na(x))) %>% 
  '['(mat.asex.avg, ,.) %>% dist
mat.summary.avg <- 2:ncol(mat.asex.complete) %>% sapply(function(i.col) {
  c(rep(1, i.col - 1), rep(2, ncol(mat.asex.complete) - i.col + 1)) %>% 
    cluster.stats(d = dist.asex.avg, clustering = ., 
                  silhouette = TRUE,
                  wgap=FALSE, 
                  sepindex=FALSE,
                  sepwithnoise=FALSE) %>% 
    '['(c('within.cluster.ss', 'avg.silwidth', 'ch')) %>% 
    unlist %>% 
    return
})
df.summary.avg <- data.frame(t(mat.summary.avg), 
                         cluster = factor(colnames(mat.asex.complete)[-1],
                                          levels = colnames(mat.asex.complete)[-1])
)
pdf('result/Westenberger/asex_break_summary_avg.pdf', width = 12, height = 10)
df.summary.avg %>% 
  gather(key = statistic, value = value, within.cluster.ss:ch) %>% 
  ggplot(aes(x = as.numeric(cluster), y = value)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name = 'cluster', 
                     breaks = as.numeric(df.summary.avg$cluster), 
                     labels = levels(df.summary.avg$cluster)) +
  facet_wrap(~statistic, ncol = 1, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=4),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# by stage averaging----------------------------------------------------
df.exprs.avg <- cbind( 
  c('rings', 'IG', 'MG', 'commit') %>% 
    sapply(function(i.stage) {
      df.clst %>% filter(stage == i.stage) %>% '['(, 'cluster') %>% 
        '['(l.clst.pv, .) %>% unlist %>% unique %>% 
        '['(df.exprs.Westenberger, ., ) %>% 
        apply(2, mean) %>% return
    }) %>% data.frame,
  c('circulation', 'sequestration') %>% 
    sapply(function(i.seq) {
      df.clst %>% filter(stage == 'asex', sequestration == i.seq) %>% '['(, 'cluster') %>% 
        '['(l.clst.pv, .) %>% unlist %>% unique %>% 
        '['(df.exprs.Westenberger, ., ) %>% 
        apply(2, mean) %>% return
    }) %>% data.frame
)

df.exprs.avg.long <- gather(df.exprs.avg, 
                            key = stage,
                            value = mean_expression,
                            rings:sequestration) %>% 
  mutate(stage = factor(stage, 
                        levels=c('circulation', 'sequestration', 'rings', 
                                 'IG', 'MG', 'commit')),
         asex_vs_gam = ifelse(stage %in% c('circulation', 'sequestration'),
                              'asex', 'gam'))

pdf('result/Westenberger/boxplot_stage_avg.pdf', width=8, height=6)
ggplot(df.exprs.avg.long,
       aes(x=stage, y=mean_expression)) +
  geom_boxplot() +
  # geom_point(aes(color=sample_type), position = position_jitterdodge()) +
  facet_grid(.~asex_vs_gam, scale='free_x', space='free_x') +
  # coord_fixed(ratio=1) +
  theme_bw()
dev.off()


