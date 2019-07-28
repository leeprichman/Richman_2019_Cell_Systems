library(antigen.garnish)
library(magrittr)
library(data.table)

# Chheda et al. J Ex Med 2018
# mutant H3.3 K27M RMSAPSTGGV, immunogenic
# wt H3.3 RKSAPSTGGV, non-immunogenic

# input for antigen.garnish
dt <- data.table::data.table(pep_mut = c("RMSAPSTGGV", "RKSAPSTGGV"),
                      sample_id =  c("mut", "wt"),
                      mutant_index = 2,
                      MHC = "HLA-A*02:01")

# none of these actually bind  < 500 oddly enough, set binding cutoff to let analysis continue
# remove_wt == FALSE to  prevent wt seq from getting tossed
dto <- dt %>% garnish_affinity(remove_wt = FALSE, blast = TRUE,
                                fitness = FALSE, binding_cutoff = 50000) %>%
# drop the wt peptides passed to output for DAI
.[pep_type != "wt"]

# dissimilar wild-type subpeptides?
dto[sample_id  ==  "wt", dissimilarity > 0] %>% any
# FALSE

# dissimilar mutant subpeptides
dto[sample_id  ==  "mut", dissimilarity > 0] %>% any
# TRUE
dto[dissimilarity > 0, nmer %>%  unique] %>% print
# MSAPSTGG RMSAPSTG have > 0 dissimilarity

dto[nmer %chin% c("MSAPSTGG", "RMSAPSTG"), dissimilarity] %>% dput
# corresponding wildtypes KSAPSTGG RKSAPSTG have dissimilarity ==  0
dto[nmer %chin% c("KSAPSTGG", "RKSAPSTG"), dissimilarity] %>% dput

# we use dput so that the R interpreter doesn't clean  up near 0 values as 0 like with print
# you can also check this with == 0
