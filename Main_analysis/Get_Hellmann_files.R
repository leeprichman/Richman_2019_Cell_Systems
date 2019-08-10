library(magrittr)
library(data.table)
library(antigen.garnish)
library(stringr)

# source our modified summary function
# identical to antigen.garnish::garnish_summary but returns overlap for venn diagrams and excludes metrics irrelevant to this analysis.
source("garnish_venn.R")

# metadata file with get URLs from European Variation Archive
# deposited by Hellmann et al. Cancer Cell 2018:
# https://www.sciencedirect.com/science/article/pii/S1535610818301235?via%3Dihub#sec5
# retrieved from: https://www.ebi.ac.uk/eva/?eva-study=PRJEB24995
dt <- data.table::fread("PRJEB24995.txt")

# get our vcf files list
v <- dt[, submitted_ftp] %>% strsplit(split = ";") %>% unlist()

v <- v[which(v %like% "vcf\\.gz$")]

options(mc.cores = parallel::detectCores())

# download in parallel and iterate over files list
# if not on linux replace wget with curl piped to file
parallel::mclapply(v, function(i){

  paste("wget ", "ftp://", i, sep = "") %>% system

})

# make sure this is all 43 published files before you proceed
v <- list.files(pattern = "vcf.gz")
length(v)

# unzip and  annotate with SnpEff
# you can install SnpEff from:
# https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip?r=http%3A%2F%2Fsnpeff.sourceforge.net%2F&ts=1501375106&use_mirror=cfhcable
#  requires Java install
parallel::mclapply(v, function(i){

  paste("gzip -d ", i, sep = "") %>% system

})

v <- list.files(pattern = "vcf")

# annotate the VCFs
lapply(v, function(i){

  print(which(i == v))
  print(length(v))

  system(paste("java -jar ",
              file.path(getwd(), "snpEff/snpEff.jar"),
              " ",
              "hg19",
              " ",
              i,
              " > ",
              i %>% stringr::str_replace("\\.vcf$", ".ann.vcf"),
              sep = ""))

})

v <- list.files(pattern = "ann\\.vcf")

# Correct the genotype field headings to work with antigen.garnish
l <- lapply(v, function(i){

  vcf <- vcfR::read.vcfR(i)

  # make sure that you  have peaked in on the files to make sure the headings are in this order
  colnames(vcf@gt) <- c("FORMAT", "TUMOR", "NORMAL")

  o <- paste("corrected_", i, ".gz", sep  = "")

  vcfR::write.vcf(vcf, file = o)

  return(o)

}) %>% unlist

# get input variants table
dt <- antigen.garnish::garnish_variants(l)

# only keep 15 alt allele counts or greater, per Hellmann et al. methods
dt <- dt[TUMOR_AD_alt >= 15]

# now we pair the patient HLA to the table of variants
met <- data.table::fread("Hellmann_meta.txt")

hla <- met[, .SD %>% unique, .SDcols = c("Patient ID", names(met)[which(names(met) %like% "allele")])]

hla %<>% melt(id.vars = "Patient ID")

# reformat alleles to work with antigen.garnish syntax, see `antigen.garnish::list_mhc()`
hla[, value := paste("HLA-", value, sep = "")]

# convert to space separated string for antigen.garnish input
hla[, MHC := paste(value %>% unique, collapse = " "), by = "Patient ID"]

hla <- hla[, .SD %>% unique, .SDcols = c("Patient ID", "MHC")] %>%
  data.table::setnames("Patient ID", "sample_id")

hla[, sample_id := sample_id %>% as.character]

# extract patient name from vcf name for matching
dt[, sample_id := sample_id %>% stringr::str_extract("(?<=corrected_)[0-9]+")]

dt <- merge(dt, hla, by = "sample_id")

# numeric only sample ids mess up final merge in garnish_affinity function so make it alphanumeric
dt[, sample_id := paste("Patient_", sample_id, sep = "")]

dt %>% data.table::fwrite("Hellmann_input.txt", sep  = "\t")

# TMB was taken determined after this filter and SnpEff annotation
# ensures that TMB refers to variants that we are confident are coding variants
# correlates with but does not match exact published data in Hellmann et al.
# get by all unique combinations of CHROM, POS, REF, ALT in table per sample_id
# example: tmb_tab <- dt[, paste(CHROM, POS, REF, ALT, sep = "_") %>% unique %>% length, by = "sample_id"]

# chunking this by sample_id so we don't run out of memory
# warning: if you don't some steps may seg fault but not error, so disregard at your own risk if you intend to parallelize at all
# it is also likely that splitting by sample_id will not be sufficient on most servers so
# consider further chunking into 10-20 row increments to allow max core usage
options(mc.cores = parallel::detectCores())

dt %<>% split(by = "sample_id")

dtl <- lapply(dt %>% seq_along, function(i){

  print(paste(i, "of", length(dt)))

  dti <- dt[[i]] %>% unique

  dto <- garnish_affinity(dti)

  # to save HDD space, can choose to only return binders here for further analysis with return(dto[Ensemble_score < 5000])
  return(dto)

})

dtl %<>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

dtl[Ensemble_score < 500] %>% data.table::fwrite("Hellmann_MHC_binders.txt")

dtl[Ensemble_score < 5000 & Ensemble_score > 1000] %>%
  data.table::fwrite("Hellmann_non_binders.txt")

dt <- dtl %>% garnish_venn

dt[, source := "Hellmann"]

dt <- merge(dt,
            data.table::fread("TMB_table.txt"),
            by = c("sample_id", "source"),
            all.x = TRUE)

dt %>% data.table::fwrite("Hellmann_output.txt", sep = "\t", quote = FALSE, row.names = FALSE)
