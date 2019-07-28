library(gglogo)
library(ggplot2)
library(magrittr)
library(data.table)

#  set your parallelization here
options(mc.cores = parallel::detectCores())

# color scheme
CDNcol <- "#ff6d00"
ADNcol <- "#2962ff"
IEcol <- "#00c853"
stcol <- "#DE7AA7"


# function to return median + CI

median_95 <- function(v){

  if (length(v) == 1)
    return(data.table::data.table(y = 1, ymin = 1, ymax = 1))

  wt <- wilcox.test(v, alternative =  "two.sided", conf.level = 0.95, conf.int = TRUE, paired = FALSE)

  y <- wt$estimate

  ci <- wt$conf.int

  ymin <- min(ci)

  ymax <- max(ci)

  return(data.table::data.table(y, ymin, ymax))

}

# set up our classification colors
colpal <- c(CDNcol, ADNcol, IEcol, stcol)
names(colpal) <- c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens")

# read in table of all neos with binding < 500nM by ensemble method, with UUID for each to pass to blast
dt <- "all_MHC_binders_w_id.txt" %>% data.table::fread

# first lets classify this list
dt[, c("CDN", "ADN", "dissimilar", "IEDB") := 0]

dt[Ensemble_score < 50, CDN := 1]
dt[min_DAI > 10, ADN := 1]
dt[dissimilarity > 0.75, dissimilar := 1]
dt[iedb_score > 0.9, IEDB := 1]

# analyze only neos exclusive to one class
dt[CDN == 1 & ADN == 0 & IEDB == 0 & dissimilar == 0, classification := "CDNs"]
dt[CDN == 0 & ADN == 1 & IEDB == 0 & dissimilar == 0, classification := "ADNs"]
dt[CDN == 0 & ADN == 0 & IEDB == 1 & dissimilar == 0, classification := "IEDB high neoantigens"]
dt[CDN == 0 & ADN == 0 & IEDB == 0 & dissimilar == 1, classification := "high dissimilarity neoantigens"]

# Logo plots
plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            plot.title = ggplot2::element_text(face = "bold",
                hjust = 0.5, size = ggplot2::rel(0.5)),
              #  axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7), color = "black", face = "bold"),
            axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.8),
                color = "black", face = "bold"),
            axis.title = ggplot2::element_text(size = ggplot2::rel(1.4),
                face = "bold"),
                axis.line.x = element_line(size = 2),
                strip.text = element_text(size = ggplot2::rel(0.7), face = "bold"),
                legend.title = element_blank(),
            legend.position = "bottom")

          return(pt)

}


# all MHC
dt2 <- data.table::copy(dt[!is.na(classification) & nchar(nmer) == 9])

dt2[, classification := paste(classification, "all MHC")]

pl <- ggplot(data = ggfortify(dt2[!is.na(classification) & nchar(nmer) == 9], "nmer", treatment = "classification")) +
       geom_logo(aes(x=position, y=bits, group=element,
          label=element, fill=interaction(Polarity, Water)),
          alpha = 0.6)  +
          facet_wrap(~classification) +
          plot_theme() +
       ggplot2::scale_fill_manual("Amino-acids properties",
      values = c("#d50000",
    "#ff8a80",
    "#ff6d00",
    "#ffd180",
    "#304ffe",
    "#8c9eff",
    "#64dd17",
    "#ccff90")) +
    ylab("bitscore") +
    xlab("\nAmino acid position") +
    theme(axis.text.y = element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4, byrow = TRUE, keywidth = 1, keyheight = 1))

ggsave(plot = pl, filename = "Logo_all.pdf", width = 9, height = 8)

# common  alleles in data
lapply(c("HLA-A*02", "HLA-B*15", "HLA-B*07", "HLA-A*03", "HLA-A*68", "HLA-C*12"), function(p){

  # replace literal asterisk with escaped asterisk
  patt <-  p %>% stringr::str_replace("\\*", stringr::fixed("\\\\*"))

  dt2 <- data.table::copy(dt[!is.na(classification) &
                            nchar(nmer) == 9 &
                            MHC %like% patt])

  dt2[, classification_f := paste(classification, "\n", p, sep = "")]

  # append number of peptides  to strip label
  nt <- dt2[, .N, by = "classification_f"]

  dt2  <- merge(dt2, nt, by = "classification_f")

  dt2[, classification_f := paste(classification_f, ", n =", " ", N, sep = "")]

  lvls <- dt2[, classification_f %>% unique]

  # put it in order
  lvls <- c(lvls[which(lvls %like% "CDN")],
            lvls[which(lvls %like% "ADN")],
            lvls[which(lvls %like% "IEDB")],
            lvls[which(lvls %like% "similar")])

  # order factor levels
  dt2[, classification_f := factor(classification_f,
                                  levels = lvls,
                                  ordered = TRUE)]

  pl <- ggplot(data = ggfortify(dt2, "nmer", treatment = "classification_f")) +
         geom_logo(aes(x=position, y=bits, group=element,
            label=element, fill=interaction(Polarity, Water)),
            alpha = 0.6)  +
            facet_wrap(~classification_f) +
            plot_theme() +
         ggplot2::scale_fill_manual("Amino-acids properties",
        values = c("#d50000",
      "#ff8a80",
      "#ff6d00",
      "#ffd180",
      "#64dd17",
      "#ccff90",
      "#304ffe",
      "#8c9eff")) +
      ylab("bitscore") +
      xlab("\nAmino acid position") +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(size = rel(1.2))) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4, byrow = TRUE, keywidth = 2, keyheight = 2)) +
      theme(legend.key = element_rect(size = 0))

  # add spacing to the legend for legibility
  pldt <- pl[[1]] %>% data.table::as.data.table

  pldt[, Polarity := paste(" ", Polarity, sep = "")]
  pldt[, Water := paste(Water, " ", sep = "")]

  pl[[1]] <- pldt

  fn <-  paste("Logo_", p, ".pdf", sep = "")

  ggsave(plot = pl, filename = fn, width = 12, height = 10)

})

file.copy("Logo_HLA-A*02.pdf", "Fig3D.pdf")

# hydropathy plots
# load Kyte-Doolittle hydropathy indices
kd <- "Kyte_Doolittle.txt" %>% data.table::fread

hp <- dt[!is.na(classification) & nchar(nmer) == 9, .SD %>% unique, .SDcols = c("nmer", "classification")]

# this is all 9mer non-binders (1000-5000) from the data, nmer, MHC allele, Ensemble_score, and sample_id
nb <- "all_non_binders_w_names.txt" %>% data.table::fread %>%
  .[, .SD %>% unique, .SDcols = c("nmer", "classification")]

hp <- data.table::rbindlist(list(hp, nb))

# bootstrap 2,000 times to ensure hydropathy result is not downsampling artifact
aal <- parallel::mclapply(1:2000, function(s){

  set.seed(s)

  print(s)

  # sample at half the size of the smallest group (dissimilarity neos)
  n <- nrow(hp[classification == "high dissimilarity neoantigens"]) / 2 %>% floor

  hi <- hp %>% data.table::copy

  hi %<>% split(by = "classification")

  # split by classification then sample v rows
  di <- lapply(hi %>% seq_along, function(i){

    v <- hi[[i]] %>% nrow

    return(hi[[i]][sample(1:v, size = n)])

  }) %>% data.table::rbindlist(use.names = TRUE)

  # melt nmers into single AA and position index
  l <- lapply(di[, nmer %>% unique], function(n){

    nl <- n %>% strsplit(split = "") %>% unlist

    i <- seq_along(nl)

    return(data.table::data.table(pos = i, AA = nl, nmer = n))

  }) %>% data.table::rbindlist

  di <- merge(l, di, by = "nmer")

  di <- merge(di, kd, by = "AA")

  # compute our median hyropathy  for the iteration per class and position
  di[, med := median(Hydropathy_index), by = c("classification", "pos")]

  di <- di[, .SD %>% unique, .SDcols = c("med", "classification", "pos")]

  di[, seed := s]

  return(di)

}) %>% data.table::rbindlist(use.names = TRUE)

# run Kruskal at each AA position (these are independent so appropriate)
pt <- lapply(aal[, pos  %>% unique],  function(n){

  a <- aal[pos == n]

  a[, classification := factor(classification)]

  av <- kruskal.test(med ~ classification,  a)$p.value

    return(data.table::data.table(pos = n, pval = av))

}) %>% data.table::rbindlist(use.names = TRUE)

# all global tests reject null, now run local test for dissimilarity neos
lapply(aal[, pos  %>% unique],  function(n){

  a <- aal[pos == n]

  a[, classification := factor(classification)]

  av <- pairwise.wilcox.test(a[, med], a[, classification],
        p.adjust.method = "bonferroni")

    return(av)

})
# for all comparisons at all positions, dissimilarity p-value is < 0.001

# factor position for plotting
aal[, pos := factor(pos)]

#  get our 95% CI for medians
aal[, c("med_med", "ymin", "ymax") := median_95(med), by = c("pos", "classification")]

# add signif where KW passes global hypothesis test
drv <- pt[pval  < 0.05, pos]

# replace NA CIs, meaning all AAs had same hydropathy, with median value
aal[is.na(ymin), ymin := med_med]
aal[is.na(ymax), ymax := med_med]

# set up to annotate signif comparisons
aal[pos %in% drv, lab := "***"]
aal[pos %in% drv, lab_height := max(ymax) + 0.2, by = "pos"]

# make plot table, factor variable, fix aesthetics
hyp <- aal[, .SD %>% unique, .SDcols = c("pos", "med_med", "ymin", "ymax", "classification", "lab", "lab_height")]
hyp[, classification := factor(classification,
    levels = c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens", "Non-binders"),
    ordered = TRUE)]

levels(hyp$classification) <- paste(" ", levels(hyp[, classification]), " ", sep = "")

colpal2 <- c(colpal, "grey")
names(colpal2) <- c(" CDNs ", " ADNs ", " IEDB high neoantigens ", " high dissimilarity neoantigens ", " Non-binders ")

hyplot <- ggplot(hyp, aes(x = pos, y = med_med)) +
    geom_point(aes(fill =  classification), color = "black", shape = 23, size = 5, stroke = 2) +
    geom_line(aes(col = classification, group = classification), linetype = "solid", size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, group = classification), width  = 0.2, size = 1, color = "black") +
    geom_text(aes(label = lab, y = lab_height), color = "black", size = 12) +
    scale_color_manual(values = colpal2) +
    scale_fill_manual(values = colpal2) +
    scale_y_continuous(limits = c(NA, 4.2)) +
    ylab("Hydropathy\n") +
    xlab("\nAmino acid position") +
    plot_theme() +
    theme(axis.line = element_line(size = 2),
          axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.4), color = "black", face = "bold"),
          axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.4), color = "black", face = "bold"),
          axis.title = element_text(size = rel(1.2), face = "bold"),
          legend.position = "bottom",
        legend.text = element_text(size = rel(0.8), face = "bold"))  +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4, byrow = TRUE, keywidth = 2, keyheight = 2),
  color = ggplot2::guide_legend(ncol = 2, byrow = TRUE, keywidth = 2, keyheight = 2))

ggsave(plot = hyplot, filename = "Fig3E.pdf", width  = 12, height = 8)

# that hydropathy plot made us concerned high dissimilarity neos just bind better because of hydrophobic bias
dt[, median(Ensemble_score), by = "classification"]
# they don't, lets graph it

# factoring, adding new lines
plot_dt <- dt[!is.na(classification)] %>% data.table::copy
plot_dt[, classification := classification  %>% stringr::str_replace("neoantigens", "\nneoantigens")]
plot_dt[, classification := factor(classification,
                            c("CDNs", "ADNs",
                            "IEDB high \nneoantigens",
                            "high dissimilarity \nneoantigens"),
                            ordered = TRUE)]

# colors
colpal2 <- colpal
names(colpal2) <- c("CDNs", "ADNs", "IEDB high \nneoantigens", "high dissimilarity \nneoantigens")

# pairwise tests and format text labels for significance bars
statd <- pairwise.wilcox.test(plot_dt[, Ensemble_score], plot_dt[, classification],
  p.adjust.method = "bonferroni", paired = FALSE)$p.value %>%
  melt %>%
  data.table::as.data.table %>%
  stats::na.omit %>%
  .[Var1 == "high dissimilarity \nneoantigens"] %>%
  .[, value := value %>% signif(digits = 2) %>% as.character]

statd[, value := paste("p =", value)]

statd[value == "p = 0", value := "p < 2.2e-16"]

# get heights for signif bars and text
h <- plot_dt[, max(Ensemble_score), by = "classification"] %>%
  data.table::setnames("classification", "Var1")

statd <- merge(statd, h, by = "Var1")

statd[, lheight := V1 * c(1.05, 1.2, 1.35)]
statd[, theight := V1 * c(1.12, 1.27, 1.42)]

statd[, txpos := c(2.5, 3, 3.5)]

# plot the median binding boxplot
plot <- ggplot(plot_dt, aes(x = classification, y = Ensemble_score)) +
  geom_boxplot(aes(fill = classification), size = 2, col = "black")  +
  geom_segment(data = statd, aes(x = Var1, xend = Var2, y =lheight, yend = lheight), size = 2) +
  geom_text(data = statd, aes(x = txpos, y = theight, label = value),
            hjust = 0.5, size = 8, fontface = "bold") +
  plot_theme() +
  xlab(NULL) +
  ylab("MHC affinity (nM)") +
  scale_fill_manual(values = colpal2) +
  theme(axis.line = element_line(size = 2),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.4), color = "black", face = "bold"),
        axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.4), hjust = 1,
                                          color = "black", face = "bold", angle = 30),
        axis.title = element_text(size = rel(1.2), face = "bold"),
        legend.position = "none")

ggsave(plot = plot, "../Figure S/binding_boxplot.pdf", width = 10, height = 8)
# now set up to run blast to generate alignments for all neos.
# below is source code adapted from the internals of antigen.garnish::garnish_dissimilarity
# output is also provided here in all_align_lengths.txt so skip ahead to after 478

# generate fastas to query
nt <- dt[, .SD %>% unique, .SDcols = c("nmer", "uuid")]
v <- nt[, nmer %>% unique]
names(v) <- nt[, uuid %>% unique]

AA <- Biostrings::AAStringSet(v, use.names = TRUE)
Biostrings::writeXStringSet(AA, file = "all_MHC_binders.fa", format = "fasta")

# run blastp for self matches
# https://www.ncbi.nlm.nih.gov/books/NBK279684/
# flags here taken from Luksza et al. Nature 2017:
# -matrix use BLOSUM62 substitution matrix
# -evalue expect value for saving hits
# -gapopen, -gapextend, numeric cost to a gapped alignment and
# -outfmt, out put a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
# length of overlap, number of mismatches, percent identical, expected value, bitscore
# you will need to add blastp to path for this.  If you install antigen.garnish correctly, this should already be the case.

# use the human db provided with antigen.garnish install
db <- "-db antigen.garnish/human.bdb"

# blast call with the parameters from Luksza model, E value, matrix, and gap costs
system(paste0(
      "blastp -query all_MHC_binders.fa ", db, " -evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -out blastp_self.csv -num_threads ", parallel::detectCores(),
      " -outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
      ))

# after broad blast, compute local alignment scores
SW_align <- function(col1,
                      col2,
                      gap_open = -11,
                      gap_extend = -1){

scores <-  parallel::mclapply(col1 %>% seq_along, function(i){

        if (i %in% seq(1, length(col1), by = (length(col1) / 10) %>% floor))
          print(paste(i, "of", length(col1), "alignments"))

        aa1 <- Biostrings::AAString(col1[i])
        aa2 <- Biostrings::AAString(col2[i])

        al <- Biostrings::pairwiseAlignment(aa1, aa2,
                                  substitutionMatrix = "BLOSUM62",
                                  gapOpening = gap_open,
                                  gapExtension = gap_extend,
                                  type = "local",
                                  scoreOnly = TRUE)

        if (length(al) == 0) al <- as.numeric(NA)

        return(al)

    }) %>% unlist

    return(scores)

  }

# read in blast output
blastdt <- list.files(pattern = "blastp_self\\.csv")

# 15gb table, don't do this on your local machine probably
# a compressed processed table is provided in the repo
blastdt <- blastdt %>% data.table::fread

blastdt %>% data.table::setnames(names(.),
                                        c("uuid",
                                        "self_anno",
                                        "nmer",
                                        "q_start",
                                        "q_stop",
                                        "WT.peptide",
                                        "s_start",
                                        "s_end",
                                        "overlap_length",
                                        "mismatch_length",
                                        "pident",
                                        "evalue",
                                        "bitscore"))
# shrink table for memory constraints
blastdt <- blastdt[, .SD %>% unique, .SDcols = c("uuid", "self_anno", "WT.peptide", "nmer")]

# remove gapped alignments, necessary given permissive blast params
blastdt <- blastdt[, WT.peptide := WT.peptide %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
                    .[, nmer := nmer %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
                    .[!is.na(WT.peptide) & !is.na(nmer)] %>%
                    .[nchar(nmer) == nchar(WT.peptide)]

# get table of alignment lengths per nmer
blastdt[, align_l := nchar(WT.peptide)]

al <- blastdt[, .N, by = c("uuid", "align_l")]

al <- merge(al, dt[, .SD %>% unique, .SDcols = c("uuid", "classification")],
            all = TRUE, by = "uuid")

# compressed copy of this is present in this repository
# uncompress with tar and pixz and then start from next line in script with with:
# al <- data.table::fread("all_align_lengths.txt")
al %>% data.table::fwrite("all_align_lengths.txt", sep = "\t")

# again only focus on single classification neos, take median # of alignments per neoantigen by length
lin <- al[classification != "", median_95(N), by = c("classification", "align_l")]

# run KW test at each alignment length (these are independent so appropriate)
# see paper methods for rational for alignment length choices
ptg <- lapply(6:14,  function(n){

  a <- al[align_l == n & classification != ""]

  a[, classification := factor(classification)]

  av <- kruskal.test(N ~ classification,  a)$p.value

    return(data.table::data.table(align_l = n, pval = av))

}) %>% data.table::rbindlist(use.names = TRUE)

# all global tests reject null, now run local test for dissimilarity compared to others
ptl <- lapply(6:14,  function(n){

  a <- al[align_l == n & classification != ""]

  a[, classification := factor(classification,
                              c("CDNs", "ADNs",
                              "IEDB high neoantigens",
                              "high dissimilarity neoantigens"),
                              ordered = TRUE)]

  av <- pairwise.wilcox.test(a[, N], a[, classification],
        p.adjust.method = "bonferroni")$p.value

  # take the largest (least significant) p-value from any comparison the plot
  av <- av["high dissimilarity neoantigens", ] %>% max

  return(data.table::data.table(align_l = n, pval = av))

}) %>% data.table::rbindlist(use.names = TRUE)

# signif labels
drv <- ptl[pval > 0.05, align_l]
drk3 <- ptl[pval < 0.001, align_l]
drk2 <- ptl[pval < 0.01 & pval > 0.001, align_l]
drk1 <- ptl[pval < 0.05 & pval > 0.01, align_l]

# set up to annotate n.s. comparisons
lin[align_l %in% drv, lab := "n.s."]
lin[align_l %in% drk1, lab := "*"]
lin[align_l %in% drk2, lab := "**"]
lin[align_l %in% drk3, lab := "***"]
lin[, lab_height := max(y) + 10, by = "align_l"]

ant <-  lin[, .SD %>% unique, .SDcols =  c("align_l", "lab", "lab_height")]

lin[, classification := factor(classification,
    levels = c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens"),
    ordered = TRUE)]

# add spaces for legend aesthetics
levels(lin$classification) <- paste(" ", levels(lin[, classification]), " ", sep = "")

colpal2 <- colpal
names(colpal2) <- c(" CDNs ", " ADNs ", " IEDB high neoantigens ", " high dissimilarity neoantigens ")
# alignments per nmer for each length
ap <- ggplot(lin, aes(x = align_l, y = y)) +
    geom_line(aes(col = classification), linetype = "solid", size = 2) +
    geom_point(aes(fill =  classification), color = "black", shape = 23, size = 5, stroke = 2) +
    geom_text(data = ant, aes(label = lab, y = lab_height), size  = 12) +
    geom_errorbar(data = lin,
      aes(ymin = ymin, ymax = ymax), width  = 0.2, size = 1, color = "black") +
    scale_color_manual(values = colpal2) +
    scale_fill_manual(values = colpal2) +
    scale_x_continuous(limits = c(6, 14), breaks = 6:14) +
    ylab("Alignments\n") +
    xlab("\nLength of alignment") +
    plot_theme() +
    theme(axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.1), color = "black", face = "bold"),
          axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.1), color = "black", face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(face = "bold", size = rel(1.1)),
          axis.line.y = element_line(size = 2),
        legend.text = element_text(size = rel(0.8), face = "bold"))  +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, keywidth = 2, keyheight = 2),
  color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, keywidth = 2, keyheight = 2))

ggsave(plot = ap, filename = "Fig3C.pdf", width  = 12, height = 8)
