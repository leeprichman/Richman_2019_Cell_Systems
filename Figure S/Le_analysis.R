library(antigen.garnish)
library(magrittr)

# Le et al Science 2017 MANA sequences
# Fig 2c table of MANAs with positive elispots
seqs  <-  c("NSDITLYVY",
            "LSSVSFFLY",
            "RSFPSWWSR",
            "SLMSKKFLPL",
            "MALSYSPEY",
            "MAMTLHPF",
            "MPSAVSCF")

st <- garnish_dissimilarity(seqs, db = "human")

dput(st)

# 5/7 have dissimilarity > 0
# structure(list(nmer = c("NSDITLYVY", "LSSVSFFLY", "RSFPSWWSR",
# "SLMSKKFLPL", "MALSYSPEY", "MAMTLHPF", "MPSAVSCF"), dissimilarity = c(1.73759295840625e-09,
# 1.71161795847752e-09, 0, 0.052270523176709, 0, 0.977252049602365,
# 0.00733844074870105)), .Names = c("nmer", "dissimilarity"), class = c("data.table",
# "data.frame"), row.names = c(NA, -7L))
