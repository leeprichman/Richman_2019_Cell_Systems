library(antigen.garnish)
library(magrittr)
library(microbenchmark)
library(ggplot2)
library(data.table)

# on 8 cores, 32gb ram
# amazon t3.2xlarge
options(mc.cores = 8)
data.table::setDTthreads(8)

set.seed(42)

vcf <- "../Main_analysis/Hellmann_input.txt" %>%
  data.table::fread %>%
  .[, sample_id := "Patient_XX"]

# only benchmark on one HLA allele
vcf[, MHC := "HLA-A*24:02"]

expr <- function(v, n){

  v <- v[sample(1:nrow(v), size = n)]

  try(garnish_affinity(v, fitness = FALSE))

}

expr2 <- function(v, n){

  v <- v[sample(1:nrow(v), size = n)]

  # binding_cutoff of 0 means IEDB_score and dissimilarity will not be run.
  try(garnish_affinity(v, blast = FALSE, fitness = FALSE, binding_cutoff = 0))

}

# benchmark with quality analysis
mb1 <- microbenchmark(expr(v = vcf, 10), times = 100)

npep1  <- lapply(list.files(pattern = "ag_out"), function(f){

  n <- data.table::fread(f) %>%
  .[pep_type != "wt", paste(nmer, MHC) %>% unique %>% length]

  file.remove(f)

}) %>% unlist

# benchmark without quality analysis
mb2 <- microbenchmark(expr2(v = vcf, 10), times = 100)

npep2  <- lapply(list.files(pattern = "ag_out"), function(f){

  n <- data.table::fread(f) %>%
  .[pep_type != "wt", paste(nmer, MHC) %>% unique %>% length]

  file.remove(f)

}) %>% unlist

mb  <- list(
  mb1 %>% data.table::as.data.table %>%
  .[, x := "With quality analysis"],
  mb2 %>% data.table::as.data.table %>%
  .[, x :=  "Without quality analysis"]) %>%
  data.table::rbindlist

mb[, seconds :=  time * 1e-9]

mb[, nvar := 10]

mb[x == "With quality analysis", neos := npep1]
mb[x == "Without quality analysis", neos := npep2]

mb %>% data.table::fwrite("antigen.garnish_affinity_benchmark.txt", sep = "\t")

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            plot.title = ggplot2::element_text(face = "bold",
                hjust = 0.5, size = ggplot2::rel(0.5)),
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.9),
                color = "black", face = "bold"),
            axis.title = ggplot2::element_text(size = ggplot2::rel(1.0),
                face = "bold"),
                axis.line = element_line(size = 2),
                strip.text = element_text(size = ggplot2::rel(0.7), face = "bold"),
                legend.title = element_blank(),
            legend.position = "bottom")

          return(pt)

}

dtb <- data.table::fread("antigen.garnish_affinity_benchmark.txt")

dtb[, neosec := neos / seconds]

dtb[, x := x %>% stringr::str_replace_all("\\ ", "\n")]

# boxplot of peptide output
pl <- ggplot(dtb, aes(x = x, y = neosec)) +
  geom_boxplot(fill = "dodgerblue", size = 2, outlier.size = 2, col = "black") +
  plot_theme() +
  xlab(NULL) +
  ylab("Peptides\noutput / sec\n") +
  ylim(c(NA, 11.5))

ggsave(plot = pl, filename = "FigS1.pdf")
