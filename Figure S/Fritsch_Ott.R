# load dependencies
library(magrittr)
library(data.table)
library(antigen.garnish)

# set max cores for dissimilarity
setDTthreads(parallel::detectCores())

# get baseline rates
dt <- data.table::fread("../Figure 3/all_MHC_binders_w_id.txt")

# 8 to 10mers only
dt %<>% .[nchar(nmer) %in% 8:10]

d <- data.table::fread("../Main_analysis/combined_output.txt") %>%
.[, .SD, .SDcols = c("sample_id", "source")]

dt <- merge(dt, d, by = "sample_id", all.x = TRUE)

rt <- lapply(dt[, source %>% unique], function(s){

  d <- dt[source == s]

  rate_all <- nrow(d[dissimilarity > 0])/ nrow(d)
  rate_CDN <- nrow(d[dissimilarity > 0 &
                    Ensemble_score < 50]) /
                    nrow(d[Ensemble_score < 50])

  return(data.table::data.table(source = s, rate_all, rate_CDN))

}) %>% data.table::rbindlist

# get average rate as baseline
bl <- dt[, rate_CDN %>% mean]

# read CIR table of validated neoantigens
# http://cancerimmunolres.aacrjournals.org/content/2/6/522
# Fritsch et al. "HLA-Binding Properties of Tumor Neoepitopes in Humans"
#   Cancer Immunology Research 2014
nmers <- data.table::fread("Fritsch_CIR_2014.txt",  header = FALSE) %>%
  unlist

st <- garnish_dissimilarity(nmers, db = "human")

# N.B. none met 0.75 threshold
st[, Dissimilarity_0 := dissimilarity > 0]

rate_Fritsch <- nrow(st[Dissimilarity_0 == TRUE]) /  nrow(st)

#   Dissimilarity_0  N
# 1:            TRUE 11
# 2:           FALSE 20

# read table of melanoma antigens
dt <- data.table::fread("Wu_mel.csv")

s <- dt[, Sequence  %>% unique]  %>%
  garnish_dissimilarity(db = "human")

s %>% data.table::setnames("nmer", "Sequence")

dt <- merge(dt, s, by = "Sequence", all.x  = TRUE)

dt[, Dissimilarity_0 := dissimilarity > 0]

# N.B. none met 0.25 threshold
dt[, .N, by = "Dissimilarity_0"]

# Dissimilarity_0   N
# 1:            TRUE  59
# 2:           FALSE 108

rate_Ott_non <- nrow(dt[Dissimilarity_0 == TRUE &
  `Peptide pulsed autologous APC` == 0]) /
  nrow(dt[!is.na(Sequence) & `Peptide pulsed autologous APC` == 0])

rate_Ott_imm <- nrow(dt[Dissimilarity_0 == TRUE &
  `Peptide pulsed autologous APC` == 1]) /
  nrow(dt[!is.na(Sequence) & `Peptide pulsed autologous APC` == 1])

dt[, .N, by = c("Dissimilarity_0", "Peptide pulsed autologous APC")] %>%
.[order(`Peptide pulsed autologous APC`, Dissimilarity_0)]

#     Dissimilarity_0 Peptide pulsed autologous APC  N
# 1:           FALSE                             0 98
# 2:            TRUE                             0 51
# 3:           FALSE                             1 10
# 4:            TRUE                             1  8

# rate of dissim in <50nM binders 10074 / 39843
# rate in all binders 51877 out of 186580
# Wu list is enriched for dissimilar neos over baseline

# read in data
# taken  from Keskin et al. Nature 2019, see original excel file in this folder
# Supplementary Data 5 41586_2018_792_MOESM5_ESM.xlsx
dt <-  "Keskin_vaccines.txt" %>% data.table::fread

dt <-  dt[!is.na(MTseq)]

# include only patients who did not recieve dexamethasone per manuscript
dt <- dt[`Patient ID` == 7]

dt[, MTseq %>% unique] %>% length
# 17 unique mutant peptides tested

st <-  dt[, MTseq %>% unique] %>% garnish_dissimilarity(db = "human")

# check the CD8 immunogenic peptides that were reported in Ext data fig 2.
st[, imm := FALSE]
st[nmer %chin% c("MVNTVAGAMK", "TTAATHREK", "MVNTVAGAM"), imm := TRUE]

# use the dissimilarity cutoff that worked in Chowell dataset  (Figure 2C)
rate_Keskin_non <-  nrow(st[imm == FALSE & dissimilarity > 0]) /
                      nrow(st[imm == FALSE])

rate_Keskin_imm <-  nrow(st[imm == TRUE & dissimilarity > 0]) /
                      nrow(st[imm == TRUE])
# Le et al. Science 2017 (see Le_analysis.R)
Le_imm <- 5 / 7

m <- data.table::data.table(Study =  c("Baseline",
                                        "MMR",
                                        "GBM_non",
                                        "GBM_imm",
                                        "Mel_non",
                                        "Mel_imm",
                                        "Fritsch"),
                              rate = c(bl,
                                        Le_imm,
                                        rate_Keskin_non,
                                        rate_Keskin_imm,
                                      rate_Ott_non,
                                      rate_Ott_imm,
                                      rate_Fritsch))

m[, rate := rate %>% round(digits = 2)]

# combined p value
# MMR 5 2
# MEL 8 10
# GBM 2 3
# > binom.test(x = 15, n  = 30, p = 0.25, alternative = "two.sided")
#
#         Exact binomial test
#
# data:  15 and 30
# number of successes = 15, number of trials = 30, p-value = 0.004714
# alternative hypothesis: true probability of success is not equal to 0.25
# 95 percent confidence interval:
#  0.3129703 0.6870297
# sample estimates:
# probability of success
#                    0.5
# Balachandran et al. Nature 2017
