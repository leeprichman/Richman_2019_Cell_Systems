library(antigen.garnish)
library(magrittr)
library(data.table)

setDTthreads(parallel::detectCores())

# read in data
# taken  from Keskin et al. Nature 2018, see original excel file in this folder
# Supplementary Data 5 41586_2018_792_MOESM5_ESM.xlsx
dt <-  "Keskin_vaccines.txt" %>% data.table::fread

dt <-  dt[!is.na(MTseq)]

# include only patients who did not recieve dexamethasone per manuscript
dt <- dt[`Patient ID` == 7]

dt[, MTseq %>% unique] %>% length
# 17 unique mutant peptides tested

st <-  dt[, MTseq %>% unique] %>% garnish_dissimilarity(db = "human")

# use the dissimilarity cutoff that worked in Chowell dataset  (Figure 2C)
st[dissimilarity > 0] %>% dput
# structure(list(nmer = c("NVTNIPLLR", "TTAATHREK", "TFTTAATHR",
# "RVSIYDYKR", "MVNTVAGAM"), dissimilarity = c(1.73790004609486e-09,
# 2.04947170345804e-13, 1.99840144432528e-15, 6.79456491070596e-14,
# 1.02029495963052e-13)), .Names = c("nmer", "dissimilarity"), class = c("data.table",
# "data.frame"), row.names = c(NA, -5L))

# check the CD8 immunogenic peptides that were reported in bold in original excel file
st[nmer %chin% c("MVNTVAGAMK", "TTAATHREK")]  %>% dput
# structure(list(nmer = c("MVNTVAGAMK", "TTAATHREK"), dissimilarity = c(0,
# 2.049472e-13)), class = c("data.table", "data.frame"), row.names = c(NA,
# -2L))
# TTAATHREK has dissimilarity greater than 0
# MVNTVAGAMK does not but...
# according to Extended data figure 2, the 9mer MVNTVAGAM was also immunogenic

st[nmer %chin% c("MVNTVAGAM"), dissimilarity > 0]
# TRUE
# the two immunogenic peptides both are dissimilar
