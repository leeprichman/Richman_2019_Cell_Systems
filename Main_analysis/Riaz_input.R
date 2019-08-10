library(antigen.garnish)
library(magrittr)
library(data.table)

# need biomaRt, install via bioconductor: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
library(biomaRt)

# source our modified summary function
# identical to antigen.garnish::garnish_summary but returns overlap for venn diagrams and excludes metrics irrelevant to this analysis.
source("garnish_venn.R")

# Supplementary Table S5 from Riaz et al. Cell 2017
# https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311224-mmc5.xlsx
dt <- "Riaz_Alleles.xlsx" %>% rio::import_list(which = 1) %>% data.table::as.data.table

# format alleles for antigen.garnish
mhc <- dt[, c(`WT Allele`, `MT Allele`) %>%
        unique %>%
        paste("HLA-", ., sep = ""), by = "Patient"]

mhc[, V1 := paste(substr(V1, 1, 5),
                  "*",
                  substr(V1, 6,7),
                  ":",
                  substr(V1, 8, nchar(V1)),
                  sep = "")]

mhc <- mhc[, V1 := V1 %>% unique %>% paste(collapse = " "), by = "Patient"] %>% unique %>%
        data.table::setnames(c("Patient", "V1"), c("sample_id", "MHC"))

# Supplementary Table S3 from Riaz et al. Cell 2017
# https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311224-mmc3.xlsx
dt <- "Riaz_Mutations.xlsx" %>% rio::import_list(which = 1) %>%
      data.table::as.data.table

vdt <- dt[, .SD %>% unique, .SDcols = c("Patient", "Chromosome", "Start", "HGVS_c", "Hugo Symbol")] %>%
        data.table::setnames("Hugo Symbol", "external_gene_name")


# convert HUGO names to ensemble ID for compatibility with antigen.garnish transcript level prediction method
HUGO2ensembl <- function (ensembl.genes){

    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    getBM(values = ensembl.genes, attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"), filters = "external_gene_name", mart = ensembl,
        bmHeader = FALSE)

}

a <- vdt[, external_gene_name %>% unique] %>% HUGO2ensembl %>% data.table::as.data.table

vdt <- merge(vdt, a[, .SD %>% unique, .SDcols = c("external_gene_name", "ensembl_transcript_id")],
              by = "external_gene_name", allow.cartesian = TRUE)

vdt %>% data.table::setnames(c("Patient", "HGVS_c"),
                            c("sample_id", "cDNA_change"))
# Take TMB at this point, unique ensembl_gene_id, cDNA_change by sample_id
# example: vdt[, paste("ensembl_gene_id", "cDNA_change", sep = "_") %>% unique, by = "sample_id"]

vdt <- merge(vdt, mhc, by = "sample_id")

vdt %>% data.table::fwrite("Riaz_input.txt", sep = "\t")

# loop over samples to allow higher parallelization without memory limit
# unless you are using an extremely large server, you will likely need to further chunk this by sample
# see warnings in annotated code for Hellmann dataset (Get_Hellmann_files.R) about caution related to
# quiet seg faults returning incomplete or erroneous results without error throwing.
vdt %<>% split(by = "sample_id")

# set cores
options(mc.cores = parallel::detectCores())

dtl <- lapply(vdt %>% seq_along, function(i){

  dt <- vdt[[i]]

  dto <- garnish_affinity(dt)

  # safety check to make sure that cDNA change matches input.  Drop transcript IDs that don't match this cDNA_change.
  # transcript level input for garnish_affinity does not yet do this by default
  dt[, c_ind := cDNA_change %>% stringr::str_extract("(?<=c\\.)[0-9]+(?=[AGCT])") %>% as.numeric] %>%
  .[, prior := cDNA_change %>% stringr::str_extract("[AGCT](?=\\>)")] %>%
  .[, keep := prior == substr(coding, c_ind, c_ind)]

  dt <- dt[keep == TRUE | pep_type == "wt"]

  # delete the dummy columns from the safety check.
  dt[, c("c_ind", "prior", "keep") := NULL]

  # to save HDD space, can choose to only return binders here for further analysis with return(dto[Ensemble_score < 5000])
  return(dt)

}) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

dtl[Ensemble_score < 500] %>% data.table::fwrite("Riaz_MHC_binders.txt")

dtl[Ensemble_score < 5000 & Ensemble_score > 1000] %>%
  data.table::fwrite("Riaz_non_binders.txt")

dt <- dtl %>% garnish_venn

dt[, source := "Riaz"]

dt <- merge(dt,
            data.table::fread("TMB_table.txt"),
            by = c("sample_id", "source"),
            all.x = TRUE)

dt %>% data.table::fwrite("Riaz_output.txt", sep = "\t", quote = FALSE, row.names = FALSE)
