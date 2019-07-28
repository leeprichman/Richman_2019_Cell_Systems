library(antigen.garnish)
library(data.table)
library(magrittr)


options(mc.cores  = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

# NB you need to decompress using tar +/- xz the input files before running this
dt <- lapply(file.path("../Main_analyis",
c("Riaz_input.txt","Hellman_input.txt", "VASR_input.txt")),
            data.table::fread) %>%
            data.table::rbindlist(use.names = TRUE, fill = TRUE)

dtl <- dt %>% split(by = "sample_id")

# we run the first half of garnish_affinity to generate all predicted
# neopeptides from variants in order to compute total dissimilar peptide burden
ot <- lapply(dtl %>% seq_along, function(i){

  print(i)
  print(length(dtl))

  d <- garnish_affinity(dtl[[i]], save = FALSE,
  predict = FALSE, fitness = FALSE, blast = FALSE)

  d[, mutID := paste(CHR, POS, REF, ALT, sep = "_")]

  d %<>% .[, .SD, .SDcols = c("sample_id", "nmer", "mutID")]

  return(d)

}) %>% data.table::rbindlist(use.names = TRUE)

# 9mers only for consistency and to dramatically reduce the computational workload
ps <- ot[nchar(nmer) == 9, nmer %>% unique]

# chunk it so we don't memory crash
# write to disk to be safe too
l <- ps %>% data.table::as.data.table %>% split(1:1000)

lapply(l %>% seq_along, function(i){

  print(i)
  print(length(l))

  p <- l[[i]] %>% unlist

  st <- garnish_dissimilarity(p, db = "human")

  list.files(pattern = "blastdt_[0-9]+\\.txt$") %>% file.remove

  st %>% data.table::fwrite(paste("strangeness_", i, ".txt", sep = ""), sep = "\t")

  return(NULL)

})  %>% data.table::rbindlist(use.names = TRUE)

sto <- lapply(list.files(pattern = "strangeness"), fread) %>% rbindlist(use.names = TRUE)

sto %>%  data.table::fwrite("9mer_peptides_dissimilarity.txt", sep = "\t")

dt <- "9mer_peptides_dissimilarity.txt" %>% data.table::fread

dt <- merge(sto, dt, by = "nmer")

dt[,  Dissimilar := dissimilarity > 0]

dto <- dt[Dissimilar == TRUE, nmer %>% length, by = "sample_id"]

dto %>%
data.table::setnames("V1", "all dissimilar neopeptides")

# write table of total number of dissimilar 9mers per patient
dto %>% data.table::fwrite("all_dissimilar_neopeptides_patient_level.txt", sep = "\t")
