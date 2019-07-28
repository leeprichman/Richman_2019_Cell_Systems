library(data.table)
library(antigen.garnish)
library(magrittr)
library(ggplot2)

mdt <- data.table::fread("ensemble_pred_mtable.txt")

mdt[, rel := abs(prediction / meas - 1) * 100]

# change variable names to look pretty on plots
mdt[, tool := tool %>% stringr::str_replace("affinity\\(nM\\)_", "")]
mdt[, tool := tool %>% stringr::str_replace("_prediction", "")]
mdt[, tool := tool %>% stringr::str_replace("_pred", "")]
mdt[tool == "Ensemble_score", tool := "a.g_ensemble"]

# remove the duplicate predictions from best MHC choice
mdt <- mdt[tool != "best_netMHC"]

# take predictons with top 10% of variation between tools
mdt[!tool %chin% c("a.g_ensemble"),
  sigma := var(prediction, na.rm =  TRUE),
  by = c("nmer", "MHC")]

# change argument passed ot quantile if you want more or less than the top decile
co <- mdt[, sigma  %>% quantile(0.9, na.rm =  TRUE)]

l <- mdt[sigma > co, .SD %>% unique, .SDcols = c("nmer", "MHC")]

mdt <- merge(mdt, l, by = c("nmer", "MHC"))

mdt <- mdt[!is.na(rel)]

plot_col <- c("#c51162",
              "#aa00ff",
              "#0091ea",
              "#64dd17",
              "#ffab00",
              "#00b8d4")

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            plot.title = ggplot2::element_text(face = "bold",
                hjust = 0.5, size = ggplot2::rel(0.5)), axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7),
                color = "black", face = "bold"),
            axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.8),
                angle = 45, color = "black", hjust = 0.96, face = "bold"),
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.6),
                face = "bold"), axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8),
                face = "bold"), legend.title = ggplot2::element_blank(),
            axis.line = element_line(size = 2),
            legend.position = "bottom")

          return(pt)

}

kruskal.test(mdt[, rel], mdt[, tool] %>% as.factor, alternative = "two.sided")
# Kruskal-Wallis rank sum test
#
# data:  mdt[, rel] and mdt[, tool] %>% as.factor
# Kruskal-Wallis chi-squared = 9313.9, df = 5, p-value < 2.2e-16

pairwise.wilcox.test(mdt$rel, mdt$tool,
  p.adjust.method = "bonferroni",
  alternative = "two.sided")
#   Pairwise comparisons using Wilcoxon rank sum test
#
# data:  mdt$rel and mdt$tool
#
# a.g_ensemble mhcflurry mhcnuggets_gru mhcnuggets_lstm netMHC
# mhcflurry       < 2e-16      -         -              -               -
# mhcnuggets_gru  < 2e-16      1.0000    -              -               -
# mhcnuggets_lstm < 2e-16      1.0000    1.0000         -               -
# netMHC          5.4e-10      1.2e-09   1.7e-06        1.3e-06         -
# netMHCpan       2.9e-13      8.7e-06   0.0021         0.0035          1.0000

pl <- ggplot(mdt, aes(x = tool, y = rel)) +
  geom_boxplot(aes(fill = tool), col = "black", outlier.shape = NA, size = 2) +
  scale_fill_manual(values = plot_col) +
  annotate("text", label = "***", x = "a.g_ensemble", y = 150, size = 10) +
  plot_theme() +
  ylab("Prediction error\n") +
  xlab(NULL) +
  ggtitle("High variance pMHC") +
  theme(legend.position = "none", plot.title = element_text(size = rel(1))) +
  scale_y_continuous(limits = c(0, 170),
  breaks =  c(0, 50, 100, 150),
  labels = c("0%", "50%", "100%", "150%"))

# boxplot of prediction error on top 10% variable pMHC
ggsave(plot = pl, filename = "Fig1D.pdf")
