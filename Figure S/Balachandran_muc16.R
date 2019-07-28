library(antigen.garnish)
library(magrittr)

# Balachandran et al. Nature 2017
# Supplementary table 1
# MUC16 R15C is the reported mutant
# They did not provide the peptides in the paper so I found all R to C mutations
# then found the match to MUC16 sequence
# Only one hit so safe to assume these are the right peptides
# sssptRslm (WT) and sssptCslm  (Mut)

c("SSSPTRSLM", "SSSPTCSLM") %>% garnish_dissimilarity(db = "human") %>% dput

# structure(list(nmer = c("SSSPTRSLM", "SSSPTCSLM"), dissimilarity = c(0,
# 8.285462e-06)), class = c("data.table", "data.frame"), row.names = c(NA,
# -2L))
# Mutant peptide has dissimilarity <  1 and WT does not.
