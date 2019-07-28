library(ggplot2)
library(antigen.garnish)
library(magrittr)
library(data.table)

CDNcol <- "#ff6d00"
ADNcol <- "#2962ff"
IEcol <- "#00c853"
stcol <- "#DE7AA7"

# set up our classification colors
colpal <- c(CDNcol, ADNcol, IEcol, stcol, "#FFFFFF", "#BEBEBE")
names(colpal) <- c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens", "TMB", "all MHC binders")

# read in the Main analysis output
dtl <- data.table::fread("../Main_analysis/combined_output.txt")

n <- c("TMB", "CDNs", "ADNs",
      "all MHC binders", "high dissimilarity neoantigens",
      "IEDB high neoantigens")

# keep only the columns we care about here
dt <- dtl[, .SD, .SDcols = c("source", n)]

# convert to long table for plotting correlation with TMB
dtm2 <- dt %>% melt(id.vars = c("source", "TMB"))

# exclude one outlier for linear correlation plot (Figure 3A),
# all MHC binders for this patient was > 30k, removed for plot scale.
# stats are rank based so this has virtually no effect on the data
dtm2 <- dtm2[TMB != 7358]

# factor stuff to get correct orders on graphs.
dtm2[, variable_f :=
   factor(variable, levels = c("all MHC binders", "CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens"), ordered = TRUE)]

# theming
plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            plot.title = ggplot2::element_text(face = "bold",
                hjust = 0.5, size = ggplot2::rel(0.5)), axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7),
                color = "black", face = "bold"),
            axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.8),
                angle = 30, color = "black", hjust = 0.96, face = "bold"),
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(1),
                face = "bold"), axis.title.y = ggplot2::element_text(size = ggplot2::rel(1),
                face = "bold"), legend.title = ggplot2::element_blank(),
                axis.line = element_line(size = 2),
                strip.text = element_text(size = ggplot2::rel(0.7), face = "bold"),
                axis.ticks = element_line(size = 2),
            legend.position = "bottom")

          return(pt)

}

# make a table of SRCC and p-vals
tab <- lapply(dtm2[, variable %>% unique], function(i){

  b <- stats::cor.test(x = dtm2[variable == i, TMB],  y = dtm2[variable == i, value],
                  method = "spearman", exact = TRUE)

  a <- data.table::data.table(srcc = b$estimate,
                              p = b$p.value,
                              variable = as.character(i))

  return(a)

}) %>% data.table::rbindlist(use.names = TRUE)

# correct for multiple comparisons,  get sigfigs for plot
tab[, p := p %>% p.adjust(method = "bonferroni")] %>%
    .[, p := p %>% signif(digits = 3)] %>%
    .[, srcc := srcc %>% signif(digits = 4)]

# now decide label positioning based on values for aesthetics
lims <- dtm2[, max(value), by = "variable"]

tab <- merge(tab, lims, by = "variable")

# use 5% and 20% over y max for label position
tab[, srcc_pos := 1.2 * V1] %>%
  .[, p_pos := 1.05 * V1]

# factor for correct order
tab[, variable_f := variable %>% stringr::str_replace("\\ neoantigens", "")]
tab[, variable_f :=
     factor(variable_f, levels = c("all MHC binders", "CDNs", "ADNs", "IEDB high", "high dissimilarity"), ordered = TRUE)]

dtm2[, variable_f := variable %>% stringr::str_replace("\\ neoantigens", "")]
dtm2[, variable_f :=
          factor(variable_f, levels = c("all MHC binders", "CDNs", "ADNs", "IEDB high", "high dissimilarity"), ordered = TRUE)]

# figure 4A, TMB correlation with classifications and all MHC binders
pl2 <- ggplot(dtm2, aes(x = TMB, y = value)) + geom_point(aes(col = variable)) +
        geom_smooth(col = "black", method  = "lm") +
  facet_wrap(~variable_f, scales = "free") +
  scale_color_manual(values = colpal) +
  plot_theme() +
  geom_text(data = tab, inherit.aes = FALSE, aes(x = 3000, y = srcc_pos, label = paste("rho =", srcc)), fontface = "bold", size = rel(6)) +
  geom_text(data = tab, inherit.aes = FALSE, aes(x = 3000, y = p_pos, label = paste("p =", p)), fontface = "bold", size = rel(6)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.8)),
      axis.text.y = element_text(size = rel(0.9)),
      strip.background = element_blank(),
      panel.spacing = unit(2, "lines")
      ) +
  xlab("\nTumor mutational burden") +
  ylab("Peptides\n")

ggsave(plot = pl2, filename = "Fig4A.pdf", width = 12, height = 8)
