pf.to.pv <- function(gene, old.ID = TRUE, print.n.match = FALSE) {
  # mapping function for mapping pf genes to pv genes
  
  # pv (new identifier) to pf (old identifier) genes mapping
  rawOrtho <- read.csv( 'data/p_vivax.csv', header=T, stringsAsFactors = F )[, 1]
  l.genes.pf <- lapply( rawOrtho, function( str ) {
    str <- strsplit( str, ' ', fixed=T )[[1]]
    str_PF <- str[grepl( 'falciparum', str, fixed=T )]
    str_PF <- gsub( '\\(.*\\)', '', str_PF, perl = T )
    return( str_PF )
  } )
  l.genes.pv <- lapply( rawOrtho, function( str ) {
    str <- strsplit( str, ' ', fixed=T )[[1]]
    str_PV <- str[grepl( 'vivax', str, fixed=T )]
    str_PV <- gsub( '\\(.*\\)', '', str_PV, perl = T )
    return( str_PV )
  } )
  
  if (!old.ID) {
    # pf mapping, old identifier to new identifier
    genes.pf.mapping <- read.table( 'data/network_project/genes_PF_mapping.txt', 
                                    header=T, stringsAsFactors = F, 
                                    check.names = F,
                                    sep='\t' )
    gene <- genes.pf.mapping[, 2] %>% 
      '['(genes.pf.mapping[, 1] %in% gene)
  }
  
  ind.match <- sapply(l.genes.pf, function(i.gene) {
    any(i.gene %in% gene)
  })
  if (print.n.match) {
    l.genes.pf[ind.match] %>% unlist %>% 
      intersect(gene) %>% print
  }
  
  return( l.genes.pv[ind.match] %>% unlist %>% unique )
}