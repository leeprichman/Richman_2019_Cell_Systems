library(magrittr)
library(data.table)
library(antigen.garnish)

# set core usage
options(mc.cores = parallel::detectCores())

# source our modified summary function
# identical to antigen.garnish::garnish_summary but returns overlap for venn diagrams and excludes metrics irrelevant to this analysis.
source("garnish_venn.R")

# this is provided as a maximally compressed file using pixz and tar
# download and install pixz on linux from shell with: sudo apt install -y pixz
# decompress the file from the shell with:
# pixz -d VASR_input.tpxz
# tar -xvf VASR_input.tar

dt <- "VASR_input.tsv" %>% data.table::fread
# this is a SnpEff annotated table or variants compiled from the original data
# adapted from Rech et al. Cancer Immunol. Res. 2018
# this could similarly be done with direct annotated VCF input to antigen.garnish::garnish_variants
# get TMB here from unique CHROM, POS, REF, ALT per sample_id here
# example: dt[, paste(CHROM, POS, REF, ALT, sep = "_") %>% unique, by = "sample_id"]

# chunk and run.  See caveats Get_Hellmann_files.R script about balancing parallelizing with memory limits
dt %<>% split(by = "sample_id")

# this is a lot so lets write to disk
# you probably want to chunk this even further (10-20 rows per job) to maximize core utilization without segmentation faults
lapply(dt, function(i){

  i %>% data.table::fwrite(paste(i[, sample_id %>% unique], "_in.txt", sep = ""),
                            sep = "\t", quote = FALSE, row.names = FALSE)

})

f <- list.files("_in.txt$", full.names = TRUE)

lapply(f %>% seq_along, function(i){

  dt <- data.table::fread(f[i])

  dto <- garnish_affinity(dt)

  # to save nearly two terabytes of disk space, feel free to only save dto[Ensemble_score < 5000]
  dto %>% data.table::fwrite(f[i] %>% stringr::str_replace("_in\\.txt$", "_out.txt"), sep = "\t",
                              quote = FALSE, row.names = FALSE)

  return(NULL)

})

# now lets read them in and combine
# if these are full files this may not be possible with the contiguous memory limit, so consider reading and then dropping rows that bind  > 500nM
l <- list.files(pattern = "_out\\.txt$", full.names = TRUE)

dtl <- lapply(l, data.table::fread) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

dtl[Ensemble_score < 500] %>% data.table::fwrite("VASR_MHC_binders.txt")

dtl[Ensemble_score < 5000 & Ensemble_score > 1000] %>%
  data.table::fwrite("VASR_non_binders.txt")

dt <- garnish_venn(dtl)

dt <- merge(dt,
            data.table::fread("TMB_table.txt"),
            by = c("sample_id", "source"),
            all.x = TRUE)

dt %>% data.table::fwrite("VASR_output.txt", sep = "\t", quote = FALSE, row.names = FALSE)
