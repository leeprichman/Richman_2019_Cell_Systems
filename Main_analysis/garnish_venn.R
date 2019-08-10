## ---- garnish_venn
#' Return venn diagram ready table for garnish_predictions output.
#'
#' @param dt Data table of garnish_predictions output to summarize overlap.
#'
#' @import stringr
#' @import data.table
#' @import magrittr
#'
#' @return Data table of neoantigen types to make venn diagram figs.
#'
#' @export garnish_venn
#' @md

garnish_venn <- function(dt){

  dtl <- lapply(dt[, sample_id %>% unique], function(si){

    c <- dt[sample_id  == si]

    c <- c[pep_type != "wt" & Consensus_scores < 500]

   `all MHC binders` <-  c[, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   CDNs <- c[Consensus_scores < 50, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ADNs <- c[min_DAI > 10, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   `IEDB high` <- c[iedb_score > 0.9, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   dissimilar <- c[dissimilarity > 0.75, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ad <- c[min_DAI > 10 & Consensus_scores < 50, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ie <- c[iedb_score > 0.9 & Consensus_scores < 50, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_st <- c[dissimilarity > 0.75 & Consensus_scores < 50, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ad_st <- c[dissimilarity > 0.75 & min_DAI > 10, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ad_ie <- c[iedb_score > 0.9 & min_DAI > 10, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ag_st <- c[dissimilarity < 0.5 & min_DAI > 10, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ag_ie <- c[iedb_score > 0.5 & min_DAI > 10, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ie_st <- c[iedb_score > 0.9 & dissimilarity > 0.75, paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ad_ie <- c[iedb_score > 0.9 & min_DAI > 10 & Consensus_scores < 50,
                    paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ad_st <- c[dissimilarity > 0.75 & min_DAI > 10 & Consensus_scores < 50,
                    paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ie_st <- c[iedb_score > 0.9 & dissimilarity > 0.75 & Consensus_scores < 50,
                    paste(nmer, MHC, sep = "_") %>% unique %>% length]

   ad_st_ie <- c[iedb_score > 0.9 & dissimilarity > 0.75 & min_DAI > 10,
                    paste(nmer, MHC, sep = "_") %>% unique %>% length]

   cd_ad_st_ie <- c[iedb_score > 0.9 & dissimilarity > 0.75 & min_DAI > 10 & Consensus_scores < 50,
                    paste(nmer, MHC, sep = "_") %>% unique %>% length]

   outdt <- data.table::data.table(sample_id = si, CDNs = CDNs, ADNs = ADNs,
                                                  `IEDB high neoantigens` = `IEDB high`, `high dissimilarity neoantigens` = dissimilar,
                                                  ad_st, ad_ie, ad_st_ie, cd_ad,
                                                  ag_st, ag_ie,
                                                  cd_ad_st, cd_ad_ie, cd_ad_st_ie, cd_st, cd_ie, cd_ie_st, ie_st,
                                                  `all MHC binders` = `all MHC binders`)

  return(outdt)

 }) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

 return(dtl)

}
