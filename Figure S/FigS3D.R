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

# convert to long table for plotting distribution of all metrics
dtm <- dt %>% melt(id.vars = c("source"))

# factor stuff to get correct orders on graphs.
dtm[, source := factor(source, levels = c("Van Allen", "Snyder", "Riaz", "Rizvi", "Hellmann"), ordered = TRUE)]
dtm[, variable_f :=
  factor(variable, levels = c("TMB", "all MHC binders", "CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens"), ordered = TRUE)]

# set 0 to 1 so log scale doesn't exclude
dtm[value == 0, value := 1]

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

# Figure S3D boxplot distribution of all MHC binders, tmb, and neoantigen classifications
pl1box <- ggplot(dtm, aes(x = variable_f, y = value)) +
  geom_boxplot(aes(fill = variable), size = 2, outlier.size = 2, col = "black") +
  facet_wrap(~source, scales = "free_x") +
  plot_theme() +
  scale_fill_manual(values = colpal) +
  scale_x_discrete(labels = c("TMB", "all MHC binders", "CDNs", "ADNs", "IEDB high", "high dissimilarity")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.7), face = "bold"),
        axis.text.x = element_text(size = rel(0.6), angle = 45),
      axis.text.y = element_text(size = rel(0.7)),
    strip.background = element_blank()) +
  xlab(NULL) +
  ylab(NULL)

pl1box <- pl1box + theme(panel.spacing = unit(2, "lines"))

ggsave(plot = pl1box, filename = "FigS3D.pdf", width = 9, height = 8)
