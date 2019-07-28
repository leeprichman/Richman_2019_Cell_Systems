library(data.table)
library(antigen.garnish)
library(magrittr)
library(ggplot2)

# set max cores, you will need to balance this with your memory maximum.
# best practice will be to chunk the upcoming `garnish_affinity` runs to allow max core usage
options(mc.cores = parallel::detectCores())

# Retrieve the input data, complete data without any blind/withholding (also provided in this repo)
# from http://tools.iedb.org/main/datasets/
# Kim et al. 2014, retrieval date in manuscript
system("curl -fssL http://tools.iedb.org/static/main/binding_data_2013.zip > binding_data_2013.zip;
unzip binding_data_2013.zip")

dt <- "bdata.20130222.mhci.txt" %>% data.table::fread

# keep only human and colums we need, change names
dt <- dt[species %chin% c("human", "mouse")]

dt <- dt[, .SD %>% unique, .SDcols = c("mhc", "sequence")]

dt %>% data.table::setnames(c("mhc", "sequence"),
                            c("MHC", "pep_mut"))

# antigen.garnish will break into all possible 8:15mers, to ensure that covers
# the actual input peptides, force it to iterate over all AAs
dt[, mutant_index := "all"]

# mock sample_id so it runs
dt[, sample_id := "ensemble_pred"]

# dirty check for supported HLA class I formatting, downstream prediction programs
# will throw warnings and return NULL but do this to save time
# use this opportunity to also take only human input first
dth <- dt[MHC %like% "HLA-[A-Z]\\*[0-9]{2}:[0-9]{2}"]

# generate all nmer inputs, dont filter WT since we want 100% of sequences returned
# ignore fitness model because we only want predictions here
# you will probably want to chunk this and run it to maximize parallelization
# without hitting memory limits. netMHC will seg fault if you overdo it.
dtho <- garnish_affinity(dth, fitness = FALSE, remove_wt = FALSE)

# keep only the nmers that are in the original data with binding measurements (in the pep_mut column)
dtho <- dtho[nmer == pep_mut]

# same but now get mouse class I
dtm <- dt[MHC %like% "H-2-[A-Z][a-z]"]

dtmo <- garnish_affinity(dtm, fitness = FALSE, remove_wt = FALSE)

# keep only the nmers that are in the original data with binding measurements (in the pep_mut column)
dtmo <- dtmo[nmer == pep_mut]

dtl <- data.table::rbindlist(list(dtho, dtmo), use.names = TRUE, fill = TRUE)

# keep only the columns we need for our analysis
dtl <- dtl[, .SD %>% unique, .SDcols = c("nmer", "MHC", "Ensemble_score",
                                      "mhcflurry_prediction", "mhcnuggets_pred_gru",
                                      "mhcnuggets_pred_lstm", "best_netMHC", "affinity(nM)_netMHC",
                                      "affinity(nM)_netMHCpan")]

# melt it for ease of use with ggplot2
mdt <- dtl %>% melt(id.vars = c("nmer", "MHC"), variable.name = "tool", value.name = "prediction")

# add back in our actual measurements to compare to predictions
a <- "bdata.20130222.mhci.txt" %>% data.table::fread %>%
        data.table::setnames(c("mhc",  "sequence"),  c("MHC", "nmer"))

mdt <- merge(mdt, a[, .SD %>% unique, .SDcols = c("MHC", "nmer", "meas")], by = c("MHC", "nmer"), all.x = TRUE)

# save  the output.  it is also in this repo but compressed
# see the README for decompression instructions
mdt %>% data.table::fwrite("ensemble_pred_mtable.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# start here if you are not repeating the prediction portion of the analysis.
mdt <- data.table::fread("ensemble_pred_mtable.txt")

# get ratio of predicted to measured
mdt[, rel := prediction / meas]

# change variable names to look pretty on plots
mdt[, tool := tool %>% stringr::str_replace("affinity\\(nM\\)_", "")]
mdt[, tool := tool %>% stringr::str_replace("_prediction", "")]
mdt[, tool := tool %>% stringr::str_replace("_pred", "")]
mdt[tool == "Ensemble_score", tool := "a.g_ensemble"]

# will stratify by this later
mdt[, nmer_length := nchar(nmer)]
mdt[, nmer_length := as.factor(nmer_length)]

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
                angle = 30, color = "black", hjust = 0.96, face = "bold"),
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.6),
                face = "bold"), axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.6),
                face = "bold"), legend.title = ggplot2::element_blank(),
            axis.line = element_line(size = 2),
            legend.position = "bottom")

          return(pt)

}


# 9mers only
dt9 <- mdt[tool != "best_netMHC" & nmer_length == 9]

dt9[, uid := paste(nmer, MHC, sep =  "_")]

dt9n <- dt9 %>% data.table::copy %>% na.omit

# stats for the IQR on 9mer predictions
bs <- parallel::mclapply(1:2000, function(s){

  print(s)

  set.seed(s)

  wt <- dt9n[uid %chin% sample(dt9n[, uid], size = 1000)]

  l <-  wt[, IQR(rel, na.rm = TRUE), by = "tool"]

  return(l)

})  %>% data.table::rbindlist(use.names = TRUE)

bs %>% data.table::setnames("V1", "IQR")

# the boostraps for each tool should be normally distributed
anova(lm(IQR  ~  tool, data = bs))
# Analysis of Variance Table
#
# Response: IQR
#              Df  Sum Sq Mean Sq F value    Pr(>F)
# tool          5 29.0751  5.8150  9218.1 < 2.2e-16 ***
# Residuals 11994  7.5661  0.0006
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.wilcox.test(bs$IQR, bs$tool,
  p.adjust.method = "bonferroni",
  alternative = "two.sided")

# just an aside,  signif by both methods, t.test and rank-based

# Pairwise comparisons using t tests with pooled SD
#
# data:  bs$IQR and bs$tool
#
#           a.g_ensemble mhcflurry mhcnuggets_gru mhcnuggets_lstm netMHC
# mhcflurry       <2e-16       -         -              -               -
# mhcnuggets_gru  <2e-16       <2e-16    -              -               -
# mhcnuggets_lstm <2e-16       <2e-16    <2e-16         -               -
# netMHC          <2e-16       0.06      <2e-16         <2e-16          -
# netMHCpan       <2e-16       <2e-16    <2e-16         <2e-16          <2e-16
#
# P value adjustment method: bonferroni

# Pairwise comparisons using Wilcoxon rank sum test
#
# data:  bs$IQR and bs$tool
#
#         a.g_ensemble mhcflurry mhcnuggets_gru mhcnuggets_lstm netMHC
# mhcflurry       1.7e-15      -         -              -               -
# mhcnuggets_gru  < 2e-16      < 2e-16   -              -               -
# mhcnuggets_lstm < 2e-16      < 2e-16   < 2e-16        -               -
# netMHC          < 2e-16      < 2e-16   < 2e-16        < 2e-16         -
# netMHCpan       < 2e-16      < 2e-16   < 2e-16        < 2e-16         < 2e-16
#
# P value adjustment method: bonferroni

# build table for our geoms for crossbars and vlines
st <- dt9n[, list(median(rel), IQR(rel)), by  = "tool"]

st[, low := V1 - (V2 / 2)]
st[, high := V1 + (V2 / 2)]

st[, lab := V2 %>% round(digits = 2) %>% sprintf(fmt = "%.2f") %>% paste("IQR:", .)]
st[tool == "a.g_ensemble", signif := "***"]

pl <- ggplot(dt9n, aes(x = rel)) +
  geom_histogram(bins = 1000, aes(col = tool, fill = tool)) +
  geom_vline(data = st, aes(xintercept = low), color = "gray", linetype = "dashed") +
  geom_vline(data = st, aes(xintercept = high), color = "gray", linetype = "dashed") +
  geom_vline(data = st, aes(xintercept = V1), color = "black", linetype = "longdash") +
  geom_segment(data = st, aes(x = low, xend = high, y = 250, yend = 250),
    lineend = "square", color = "black") +
  geom_segment(data = st, aes(x = low, xend = low, y = 150, yend = 350)) +
  geom_segment(data = st, aes(x = high, xend = high, y = 150, yend = 350)) +
  geom_text(data = st, aes(label = lab, x = 0.20, y = 300), fontface = "bold", size = 5) +
  geom_text(data = st, aes(label = signif, x = 0.20, y = 400), fontface = "bold", size = 8) +
  facet_wrap(~tool, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = plot_col) +
  scale_color_manual(values = plot_col) +
  plot_theme() +
  scale_x_continuous(limits  = c(0.1, 3), trans = "log2", labels = function(x) x %>%round(digits = 2) %>% sprintf(fmt = "%.2f")) +
  # don't care about looking at the tails
  ylab("Predictions\n") +
  xlab("\nPredicted/measured affinity") +
  theme(legend.position = "none",
        axis.line = element_line(size = 1),
        axis.text.x = element_text(angle = 45),
        panel.spacing.x = unit(2, "lines"),
        strip.text = element_text(size = ggplot2::rel(0.7), face = "bold"),
        strip.background = element_blank(),
        axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.8)),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8)))

# plot of prediction error distributions
ggsave(plot = pl, filename = "Fig1C.pdf")

# SRCC analysis  by peptide length, get table of SRCC values for each tool
# downsample and bootstrap 2,000 times
tsl <- parallel::mclapply(1:2000, function(s){

  print(s)

  set.seed(s)

  ts <- lapply(8:15, function(n){

    adt <- mdt[nchar(nmer) == n] %>% stats::na.omit

    if (nrow(adt) == 0) return(NULL)

    adt %<>% split(by = "tool")

    adt <- lapply(adt %>% seq_along, function(i){

      v <- adt[[i]] %>% nrow

      # least abundant nmer_length in dataset  has 227 measurements
      # 100 measurements per tool is a good downsample
      return(adt[[i]][sample(1:v, size = 100)])

    }) %>% data.table::rbindlist(use.names = TRUE)


    cls <- lapply(adt[, tool %>% unique], function(i){

      # this function returns rho
      try(cor.test(x = adt[tool == i, prediction],  y = adt[tool == i, meas],
                  method = "spearman"))

    })

    names(cls) <- adt[, tool %>% unique]

    tab <- lapply(cls %>% seq_along, function(i){

      obj <- cls[[i]]

      if (class(obj) == "try-error") return(NULL)

      nm <- names(cls[i])

      # extract rho
      srcc <- obj$estimate

      return(data.table::data.table(nm, srcc))

    }) %>% data.table::rbindlist(use.names = TRUE)

    return(tab[, nmer_length := n])

  }) %>% data.table::rbindlist(use.names = TRUE)

  return(ts[, seed := s])

}) %>% data.table::rbindlist(use.names = TRUE)

#  set up some of our theming
plot_col <- c("#c51162",
    "#aa00ff",
    "#0091ea",
    "#64dd17",
    "#ffab00",
    "#00b8d4")

# manage aesthetics for SRCC plot, take single algorithms only
tsl <- tsl[nm != "best_netMHC"]

tsl[, nm_f := factor(nm, levels = c("a.g_ensemble", "mhcflurry", "mhcnuggets_gru", "mhcnuggets_lstm", "netMHC", "netMHCpan"))]

# add and nudge spacing for legend
levels(tsl$nm_f) <- paste(" ", levels(tsl[, nm_f]), "  ", sep = "")

names(plot_col) <- levels(tsl$nm_f)

# factor x so x scale is discrete on graph
tsl[, nmer_length := factor(nmer_length)]

# figure 1B: SRCC by tool by peptide length bootstraps boxplot
pl <- ggplot(tsl, aes(x = nmer_length, y = srcc)) +
    geom_boxplot(aes(fill = nm_f, color = nm_f), size = 1, outlier.size = 2, color = "black",
  position = position_dodge2(padding = 0.5)) +
    scale_y_continuous(limits = c(NA, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_fill_manual(values = plot_col) +
    scale_color_manual(values = plot_col) +
    plot_theme() +
    xlab("\nPeptide length") +
    ylab("SRCC\n") +
    theme(legend.key.width  = unit(1, "inches"),
          legend.key.height = unit(1, "inches"),
          legend.text = element_text(size = rel(0.7), face = "bold"),
          axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(1.3), angle = 0, hjust = 0.5),
          axis.title.x = element_text(size = rel(1.0)),
          axis.title.y = element_text(size = rel(1.3)),
          plot.title = element_text(size = rel(1))
        )

# plot of SRCC bootstraps
ggsave(plot = pl, filename = "Fig1B.pdf", width = 12, height = 9)
