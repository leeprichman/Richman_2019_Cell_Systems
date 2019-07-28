library(data.table)
library(magrittr)
library(ggplot2)

# read in table of all neos with binding < 500nM by ensemble method, with UUID for each to pass to blast
dt <- "all_MHC_binders_w_id.txt" %>% data.table::fread

# only 9mers for clarity
dt <- dt[nchar(nmer) == 9]

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32, color = "black"),
            axis.line = element_line(size = 2),
            axis.ticks = element_line(size = 2),
            axis.text = element_text(face = "bold", color = "black", size = rel(0.6)),
            axis.text.y = element_text(angle =  30),
            axis.title = element_text(face = "bold", color = "black"),
            legend.title = ggplot2::element_blank(),
            legend.text  = element_text(size = rel(0.6)),
            legend.position = c(0.7, 0.25),
            legend.key.width = unit(1, "cm"), legend.key.height = unit(1, "cm")
          )

          return(pt)

}

# figure of dissimilarity value densities
pl <- ggplot(dt, aes(x = dissimilarity)) +
  annotate("rect", xmin = 0.75, xmax = 1.05, ymin = 0, ymax = 1e05,
  alpha = 0.6, fill = "grey", color = "transparent") +
  annotate("text", label = "high\ndissimilarity", x = 0.9, y = 1e04,
  fontface = "bold", size = 5) +
  geom_histogram(aes(y=..count..),
    color = "black", fill = "#DE7AA7",
   bins = ceiling(log(length(dt$dissimilarity), base = 2) + 1)) +
  geom_vline(xintercept = 0.75, color = "black", linetype = "dashed", size = 1.25) +
  scale_y_log10(breaks = c(0, 10, 100, 1000, 10000, 100000)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  plot_theme() +
  xlab("\nDissimilarity") +
  ylab("Peptides\n")

ggsave(plot = pl, "Fig3A.pdf", device = "pdf")

# figure of iedb score value densities
pl <- ggplot(dt, aes(x = iedb_score)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 0, ymax = 1e05,
  alpha = 0.6, fill = "grey", color = "transparent") +
  annotate("text", label = "IEDB\nhigh", x = 1.0, y = 2e04,
  fontface = "bold", size = 5) +
  geom_histogram(aes(y=..count..),
    color = "black", fill = "#00c853",
   bins = ceiling(log(length(dt$dissimilarity), base = 2) + 1)) +
  geom_vline(xintercept = 0.9, color = "black", linetype = "dashed", size = 1.25) +
  scale_y_log10(breaks = c(0, 10, 100, 1000, 10000, 100000)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  plot_theme() +
  xlab("\nIEDB score") +
  ylab("Peptides\n")

ggsave(plot = pl, "FigS2B.pdf", device = "pdf")
