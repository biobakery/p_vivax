library(magrittr)
source('src/source.R')
# mapping of pf gametocyte stages to pv genes

group.anno <- c( 'asex', 'asexCirc', 'asexSeq',
                 'commit', 'IG', 'IGYG', 'MG',
                 'rings' )
l.anno.pf <- group.anno %>% lapply(function(i.group) {
  i.df.anno <- i.group %>% 
    paste0('data/network_project/genes_', ., '_clean.txt') %>% 
    read.table(header = F, stringsAsFactors = F)
  return(i.df.anno[, 1])
})
names(l.anno.pf) <- group.anno

##########################################################################
# some sanity check
# There are some overlapping between different annotation groups,
# possibly due to the conversion from old identifier to new ones.
# The proportions of overlapped genes seem okay to neglect (less than 10%).
# Because commit, IG, MG, and rings genes are derived as they are depicted
# in Figure 6, Pelle et al., use these groups as the annotations.

# mat.nIntersect <- group.anno %>% 
#   sapply(function(i.group) {
#     group.anno %>% sapply(function(j.group) {
#       intersect(l.anno.pf[[i.group]], l.anno.pf[[j.group]]) %>% 
#         unique %>% length %>% return
#     })
#   })
# mat.nIntersect <- sapply(l.anno.pf, length) %>% rbind(mat.nIntersect)
# mat.nIntersect[1:4, ]
# mat.nIntersect[-(2:4), -(1:3)]
##########################################################################

l.anno.pv <- l.anno.pf %>% lapply(pf.to.pv, old.ID = FALSE, print.n.match = FALSE)
names(l.anno.pv) <- group.anno

##########################################################################
# run the same sanity check as before with pv annotations

# mat.nIntersect <- group.anno %>%
#   sapply(function(i.group) {
#     group.anno %>% sapply(function(j.group) {
#       intersect(l.anno.pv[[i.group]], l.anno.pv[[j.group]]) %>%
#         unique %>% length %>% return
#     })
#   })
# mat.nIntersect <- sapply(l.anno.pv, length) %>% rbind(mat.nIntersect)
# mat.nIntersect[1:4, ]
# mat.nIntersect[-(2:4), -(1:3)]
##########################################################################