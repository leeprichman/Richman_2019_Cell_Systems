library(magrittr)
library(data.table)
library(antigen.garnish)
library(ggplot2)

# input and output from Luksza et al. reported in supplementary data files
# see Key resources table for original source
dti <- "Lukza_input.txt" %>% data.table::fread
dto <- "Lukza_output.txt" %>% data.table::fread

# take only columns we  care about.
a <- dto[, .SD %>% unique, .SDcols = c("MutantPeptide", "R")]

# we only want to compare R to iedb_score
b <- dti[, MT.Peptide %>% unique] %>% antigen.garnish::iedb_score(db = "human")

a %>% data.table::setnames("MutantPeptide", "nmer")

dt <- merge(a, b, by = "nmer", all.x = TRUE)

# our methodology for iedb_score returns NA when there are no BLAST
# alignments to IEDB, Luksza scripts return 0.  Fix that discrepancy here
# and circumvent issues with non-numeric values downstream.
dt[is.na(iedb_score), iedb_score := 0]

# save intermediate output
dt %>%  data.table::fwrite("parameters_comp.txt", sep = "\t")

# compute correlation coeff
ct <- cor.test(dt$R, dt$iedb_score, method = "spearman")

# reformat p-val for annotation
pval <- paste("p = ", ct$p.value, sep = "")

if (ct$p.value == 0)
 pval <- "p < 2.2e-16"

# plot annotation
anno <- paste("rho = ", ct$estimate %>% signif(digits = 4), "\n", pval , sep = "")

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            axis.line = element_line(size = 2),
            axis.ticks = element_line(size = 2),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom")

          return(pt)

}

pl <- ggplot(dt, aes(x = R, y = iedb_score)) +
geom_point(fill = "dodgerblue", size = 4, shape = 23, col = "black", stroke = 2) +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = 1e-5, y = 1e-31, label = anno, size = 12, fontface = "bold") + plot_theme() +
  xlab("\nLuksza et al. 2017\nTCR recognition probability") +
  ylab("antigen.garnish\nIEDB_score\n") +
  theme(
        axis.text = element_text(size = rel(1), face = "bold", color = "black"),
        axis.title = element_text(size = rel(1.3), face = "bold", color = "black")
      )

ggsave(plot = pl, filename = "Fig2A.pdf", height = 9, width = 12)
