library(pROC)
library(magrittr)
library(data.table)
library(ggplot2)

# table of dissimilarity and IEDB values for all entries
dt <- "Chowell_st_ie.txt" %>% data.table::fread

# table of predictions for those Chowell entries with 4digit HLA
b <- "Chowell_preds.txt" %>% data.table::fread  %>%
  .[, .SD, .SDcols = c("nmer", "MHC", "Ensemble_score")]

dt2 <- merge(dt, b, by = "nmer")

# function to compute mean value  across all positions in a peptide
# all 5 atchley factors and  KD hydropathy
Atchleyize <- function(nmers){

  if (class(nmers) != "character") stop("Input is character vector of peptides.")

  # in case of malformation
  nmers %<>% toupper

  # break into letters
  nl <- lapply(nmers, function(n){strsplit(n, split = "")}) %>%
    unlist(recursive  = FALSE)

  names(nl) <- nmers

  nl <- lapply(nl %>% seq_along, function(i){

    nv <- nl[[i]]

    names(nv) <- paste("Pos_", 1:length(nv), sep = "")

    nv %<>% as.matrix %>% t

    rownames(nv) <- names(nl[i])

    return(nv %>% data.table::as.data.table(keep.rownames = TRUE))

  })

  nl %<>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

  # create long t able of amino acids per position
  nl %<>% melt(id.cols = "rn",
    measure.vars = names(nl)[which(names(nl) %like% "Pos_")],
    value.name = "AA")

  # Atchley factor table. From Atchley et al. PNAS 2005
  atch <- structure(list(`AA` = c("A", "C", "D", "E", "F", "G",
    "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
    "Y"), `Factor I` = c(-0.591, -1.343, 1.05, 1.357, -1.006, -0.384,
    0.336, -1.239, 1.831, -1.019, -0.663, 0.945, 0.189, 0.931, 1.538,
    -0.228, -0.032, -1.337, -0.595, 0.26), `Factor II` = c(-1.302,
    0.465, 0.302, -1.453, -0.59, 1.652, -0.417, -0.547, -0.561, -0.987,
    -1.524, 0.828, 2.081, -0.179, -0.055, 1.399, 0.326, -0.279, 0.009,
    0.83), `Factor III` = c(-0.733, -0.862, -3.656, 1.477, 1.891,
    1.33, -1.673, 2.131, 0.533, -1.505, 2.219, 1.299, -1.628, -3.005,
    1.502, -4.76, 2.213, -0.544, 0.672, 3.097), `Factor IV` = c(1.57,
    -1.02, -0.259, 0.113, -0.397, 1.045, -1.474, 0.393, -0.277, 1.266,
    -1.005, -0.169, 0.421, -0.503, 0.44, 0.67, 0.908, 1.242, -2.128,
    -0.838), `Factor V` = c(-0.146, -0.255, -3.242, -0.837, 0.412,
    2.064, -0.078, 0.816, 1.648, -0.912, 1.212, 0.933, -1.392, -1.853,
    2.897, -2.647, 1.313, -1.262, -0.184, 1.512)), .Names = c("AA",
    "Factor I", "Factor II", "Factor III", "Factor IV", "Factor V"
    ), row.names = c(NA, -20L), class = c("data.table", "data.frame"
    ))

    # Kyte-Doolittle hydropathy
    KD <- structure(list(AA = c("I", "V", "L", "F", "C", "M", "A", "G",
    "T", "W", "S", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"),
    Factor_KD = c(4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.9, -0.8,
      -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5)),
      .Names = c("AA", "Factor_KD"),
      row.names = c(NA, -20L), class = c("data.table", "data.frame"))

    atch <- merge(KD, atch, by = "AA")

  nl <- merge(nl, atch, all.x = TRUE, by = "AA")

  names(nl) <- names(nl) %>% stringr::str_replace_all("\\ ", "_")

  nl <- nl[, .SD, .SDcols =  c("rn", "variable",
          names(nl)[which(names(nl) %like% "Factor_")])] %>%
    melt(id.vars = c("rn", "variable"),
        variable.name  = "Factor")

  nl %<>% split(by = "Factor")

  # now make wide table, one column per AA position
  nl <- lapply(nl, function(t){

      f <- t[, Factor  %>% unique]

      t %<>% dcast(rn ~ variable)

      t[, Factor := f]

  }) %>% data.table::rbindlist(use.names = TRUE)

  nl %>% data.table::setnames("rn", "nmer")

  return(nl)

}

at <-  dt[, nmer %>% unique] %>% Atchleyize

at %<>% melt(id.vars = c("Factor", "nmer")) %>% na.omit

at[, mean := mean(value), by = c("nmer", "Factor")]

at %<>% .[, .SD %>% unique,
          .SDcols = c("mean", "Factor", "nmer")]  %>%
          dcast(nmer ~  Factor, value.var = "mean")

dt <- merge(dt, at, by = "nmer")

dt %<>% melt(id.vars = c("nmer", "Immunogenicity"))

dtl <- dt %>% split(by = "variable")

# list of ROC objects for each of the factors, dissim, KD hydro
rocl <- lapply(dtl, function(t){

  r <- pROC::roc(formula = Immunogenicity ~ value,
                data = t,
                smooth = FALSE,
                percent = TRUE,
                partial.auc = FALSE,
                ci=TRUE,
                ci.method = "bootstrap",
                boot.n=2000,
                ci.alpha=0.9,
                # nearly 50/50 pos and negative so dont need to stratify
                stratified=TRUE,
                # arguments for plot
                plot=FALSE,
                auc.polygon=FALSE,
                max.auc.polygon=FALSE,
                grid=TRUE,
                # print numeric value of AUC on plot
                # (can be set with print.auc.x/y args)
                print.auc=TRUE,
                show.thres=TRUE)

  return(r)

})

names(rocl) <- lapply(dtl, function(x) x$variable %>% unique) %>% unlist

# Build a ROC object and compute the AUC for affinity too
# separate because its only on 6,050 of the entries in the database
aff <- pROC::roc(formula = Immunogenicity ~ Ensemble_score,
              data = dt2,
              smooth = FALSE,
              percent = TRUE,
              partial.auc = FALSE,
              ci=TRUE,
              ci.method = "bootstrap",
              boot.n=2000,
              ci.alpha=0.9,
              # nearly 50/50 pos and negative so dont need to stratify
              stratified=FALSE,
              # arguments for plot
              plot=FALSE,
              auc.polygon=FALSE,
              max.auc.polygon=FALSE,
              grid=TRUE,
              # print numeric value of AUC on plot
              # (can be set with print.auc.x/y args)
              print.auc=TRUE,
              print.auc.y = 30,
              show.thres=TRUE)

rocl[[9]] <- aff

names(rocl)[9] <- "Affinity (nM)"

## Coordinates of the curves ##
lapply(rocl, function(roc){
  pROC::coords(roc, "best", ret=c("threshold", "specificity", "1-npv", "sensitivity"))
})
# $dissimilarity
#   threshold specificity       1-npv sensitivity
#     1.00000    91.50691    21.07041    76.46475
#
# $iedb_score
#   threshold specificity       1-npv sensitivity
#    0.999941   96.598639   39.308380   39.721946
#
# $`Affinity (nM)`
#   threshold specificity       1-npv sensitivity
#   25.76432    82.29927    34.07810    26.28726


rocl %>% saveRDS("Chowell_ROCs.RDS")

labs <- lapply(rocl, function(i){

  (i$auc * 0.01) %>%
    round(digits = 2) %>%
    sprintf(fmt = "%.2f") %>%
    as.character

}) %>% unlist

# pull raw data so we can replot with ggplot2
rt  <- lapply(rocl %>% seq_along, function(n){

  i  <- rocl[[n]]

  data.table::data.table(y = i$sensitivities,
                        x = 100 - i$specificities,
                        thresholds = i$thresholds,
                        metric = names(rocl)[n])

}) %>% data.table::rbindlist(use.names =  TRUE)

# format labels and legend entries
labs <- paste("(AUC = ", labs, ")", sep  = "")

labs <- paste(names(rocl), labs) %>%
  stringr::str_replace("Factor.KD", "Hydropathy") %>%
  stringr::str_replace("iedb_score", "IEDB score") %>%
  stringr::str_replace("^d", "D") %>%
  stringr::str_replace("_", " ")

rt[, metric := metric %>%
factor(levels = rt[, metric %>% unique],
        labels = labs,
        ordered = TRUE)]

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32, color = "black"),
            axis.line = element_line(size = 2),
            axis.ticks = element_line(size = 2),
            axis.text = element_text(face = "bold", color = "black", size = rel(0.6)),
            axis.title = element_text(face = "bold", color = "black"),
            legend.title = ggplot2::element_blank(),
            legend.text  = element_text(size = rel(0.6)),
            legend.position = "right",
            legend.key.width = unit(1, "cm"), legend.key.height = unit(1, "cm"),
            panel.grid.major = element_line(color = "grey"))

          return(pt)

}

colors <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00",  "#00b8d4",
            "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00",
            "#aeea00", "#dd2c00")

# plot of Factors and KD hydro vs dissimilarity
pl <- ggplot(rt[!metric %like% "Aff|IEDB"], aes(x = x, y = y)) +
  #geom_point(aes(col = metric), shape = 43, size = 4) +
  geom_line(aes(col = metric), size = 2) +
  geom_abline(slope =  1, intercept = 0,
              linetype = "dashed", color = "black", size = 2) +
  scale_color_manual(values = colors[c(1:6,10)]) +
  scale_x_continuous(labels = function(x){paste0(x, "%")}) +
  scale_y_continuous(labels = function(x){paste0(x, "%")}) +
  plot_theme() +
  xlab("False positive rate") +
  ylab("True positive rate")

ggsave(pl, filename = "Fig2E.pdf", height = 7, width = 10.5)

# plot of Atchley Factors and Hydropathy vs Dissimilarity
pl <- ggplot(rt[!metric %like% "Hydro|Factor"], aes(x = x, y = y)) +
  #geom_point(aes(col = metric), shape = 43, size = 4) +
  geom_line(aes(col = metric), size = 2) +
  geom_abline(slope =  1, intercept = 0,
              linetype = "dashed", color = "black", size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(labels = function(x){paste0(x, "%")}) +
  scale_y_continuous(labels = function(x){paste0(x, "%")}) +
  plot_theme() +
  xlab("False positive rate") +
  ylab("True positive rate") +
  theme(legend.position =  c(0.7, 0.2))

ggsave(pl, filename = "Fig2C.pdf", height = 7, width = 7)

# Supplemental figure to show affinity prediction by tool vs immunogenicity
dt <- "Chowell_st_ie.txt" %>% data.table::fread

b <- "Chowell_preds.txt" %>% data.table::fread

dt2 <- merge(dt, b, by = "nmer")

dt2[, MHC := NULL]

dt2 %<>% melt(id.vars = c("nmer", "Immunogenicity"))

dtl <- dt2 %>% split(by = "variable")

rocl <- lapply(dtl, function(t){

  r <- pROC::roc(formula = Immunogenicity ~ value,
                data = t,
                smooth = FALSE,
                percent = TRUE,
                partial.auc = FALSE,
                ci=TRUE,
                ci.method = "bootstrap",
                boot.n=2000,
                ci.alpha=0.9,
                stratified=TRUE,
                # arguments for plot
                plot=FALSE,
                auc.polygon=FALSE,
                max.auc.polygon=FALSE,
                grid=TRUE,
                # print numeric value of AUC on plot
                # (can be set with print.auc.x/y args)
                print.auc=TRUE,
                show.thres=TRUE)

  return(r)

})

names(rocl) <- lapply(dtl, function(x) x$variable %>% unique) %>% unlist

## Coordinates of the curves ##
lapply(rocl, function(roc){
  pROC::coords(roc, "best", ret=c("threshold", "specificity", "1-npv", "sensitivity"))
})

rocl %>% saveRDS("Tools_ROCs.RDS")

labs <- lapply(rocl, function(i){

  (i$auc * 0.01) %>%
    round(digits = 2) %>%
    sprintf(fmt = "%.2f") %>%
    as.character

}) %>% unlist

# pull raw data so we can replot with ggplot2
rt  <- lapply(rocl %>% seq_along, function(n){

  i  <- rocl[[n]]

  data.table::data.table(y = i$sensitivities,
                        x = 100 - i$specificities,
                        thresholds = i$thresholds,
                        metric = names(rocl)[n])

}) %>% data.table::rbindlist(use.names =  TRUE)

# formatting labels and legend
labs <- paste("(AUC = ", labs, ")", sep  = "")

labs <- paste(names(rocl), labs) %>%
  stringr::str_replace("iedb_score", "IEDB score") %>%
  stringr::str_replace("^d", "D") %>%
  stringr::str_replace("affinity\\(nM\\)_", "") %>%
  stringr::str_replace("mhc", "MHC") %>%
  stringr::str_replace("_pred(iction)?", "") %>%
  stringr::str_replace("gru", "GRU") %>%
  stringr::str_replace("lstm", "LSTM") %>%
  stringr::str_replace("_", " ")

rt[, metric := metric %>%
factor(levels = rt[, metric %>% unique],
        labels = labs,
        ordered = TRUE)]

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32, color = "black"),
            axis.line = element_line(size = 2),
            axis.ticks = element_line(size = 2),
            axis.text = element_text(face = "bold", color = "black", size = rel(0.6)),
            axis.title = element_text(face = "bold", color = "black"),
            legend.text  = element_text(size = rel(0.6)),
            legend.title = element_text(size = rel(0.6), face = "bold"),
            legend.position = "right",
            legend.key.width = unit(1, "cm"), legend.key.height = unit(1, "cm"),
            panel.grid.major = element_line(color = "grey"))

          return(pt)

}

colors <- c("#c51162", "#aa00ff", "#0091ea", "#64dd17", "#ffab00",  "#00b8d4",
            "#d50000", "#6200ea", "#2962ff", "#a7ffeb", "#00c853", "#ff6d00",
            "#aeea00", "#dd2c00")

# plot of individual prediction tools vs dissimilarity
pl <- ggplot(rt[!metric %like% "IEDB|Ensemble"], aes(x = x, y = y)) +
  #geom_point(aes(col = metric), shape = 43, size = 4) +
  geom_line(aes(col = metric), size = 2) +
  geom_abline(slope =  1, intercept = 0,
              linetype = "dashed", color = "black", size = 2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(labels = function(x){paste0(x, "%")}) +
  scale_y_continuous(labels = function(x){paste0(x, "%")}) +
  plot_theme() +
  labs(color = "Peptides with\npredicted affinity") +
  xlab("False positive rate") +
  ylab("True positive rate")

ggsave(pl, filename = "FigS2A.pdf", height = 7, width = 12)
