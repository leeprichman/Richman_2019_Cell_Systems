library(grid)
library(gridExtra)
library(antigen.garnish)
library(magrittr)
library(data.table)

# run on fullli dataset, AA sequence quality metrics only.
dt <- "pnas_chowell_2015.txt" %>% data.table::fread

##  Run AG  to get predictions for those with 4 digit HLA
dt <- dt[, .SD %>% unique, .SDcols = c("Epitope", "MHC Allele")]
dt <- dt[`MHC Allele` %like% "HLA-[ABC]\\*[0-9]{2}:[0-9]{2}"]
# now we are down to 6023 from 9888

dt  %>% setnames(names(dt), c("pep_mut", "MHC"))
dt[, MHC := MHC %>% stringr::str_extract("HLA-[ABC]\\*[0-9]{2}:[0-9]{2}")]
dt[, mutant_index := "all"]
dt[, sample_id :=  "placeholder"]

# you might need to chunk this if you are running with limited RAM
dt <- antigen.garnish::garnish_affinity(dt, remove_wt = FALSE, fitness = FALSE)

# write nmer and affinity prediction for limited dataset to file.
dt[, .SD  %>% unique,
.SDcols = c("nmer",
            "MHC",
            "Ensemble_score",
            "mhcflurry_prediction",
            "mhcnuggets_pred_lstm",
            "mhcnuggets_pred_gru",
            "affinity(nM)_netMHCpan",
            "affinity(nM)_netMHC")] %>%
  data.table::fwrite("Chowell_preds.txt", sep = "\t")

# now we run sequence metrics (IEDB score and dissimilarity) on full dataset.
dt <- "pnas_chowell_2015.txt" %>% data.table::fread

dt <- dt[, .SD %>% unique, .SDcols = c("Epitope", "Immunogenicity")] %>%
    data.table::setnames("Epitope", "nmer")

# chunk it just in case of memory limit
# BLAST is ok but the SW-alignment of hits is memory intensive
dtl <- split(dt, 1:10)

stl <- lapply(dtl %>% seq_along, function(i){

  print(i)

  st <-  dtl[[i]][, nmer %>% unique] %>% antigen.garnish::iedb_score(db = "human")
  it <-  dtl[[i]][, nmer %>% unique] %>% antigen.garnish::garnish_dissimilarity(db = "human")

  st <- merge(st, it, by  = "nmer")

  return(st)

}) %>% data.table::rbindlist(use.names = TRUE)

dt <- merge(stl, dt, by = "nmer", all.y = TRUE)

dt %<>% unique(by = c("nmer", "Immunogenicity", "dissimilarity", "iedb_score"))

# if blast didn't return anything, no hits so iedb_score is 0.
# do this manually here since this is normally done within garnish_affinity
dt[is.na(iedb_score), iedb_score := 0]

# save our output just in case
dt %>% data.table::fwrite("Chowell_st_ie.txt", sep  = "\t")

# first sequence metrics table only, no excluded entries

# we are using the absolute most permissive metrics here
# This is based on both an ROC curve analysis, and because
# IEDB high and high dissimilarity neos using the cutoffs of  c(0.25, 0.9) as in other figures
# are such a relatively small percent that this analysis is not meaningful
# change these numbers to have < 0.25 and >0.9, OR and trends remain significant
# see Figure S1A for this exact analysis with those cutoffs
# overall sensitivity though is quite low, so we use the permissive threshold of 1 and 0 instead
dt <- "Chowell_st_ie.txt" %>% data.table::fread

# dissimilarity first
dt[dissimilarity > 0, res := "Dissimilarity > 0"]
dt[dissimilarity == 0, res := "Dissimilarity = 0"]

t <- dt[,  .N, by = c("res", "Immunogenicity")] %>% dcast(Immunogenicity ~ res, value.var  = "N")

tab <- t[, .SD, .SDcols = 2:3] %>% as.data.frame
rownames(tab) <-  t[, Immunogenicity]

tab <- tab[, c("Dissimilarity = 0", "Dissimilarity > 0")]

# perform local tests by checking all columns against all others combined
pdt <- lapply(colnames(tab), function(i){

  t <- tab %>% data.table::as.data.table(keep.rownames =  TRUE)

  # sum all  other columns vs iterant
  v <- names(t)[which(!names(t) %chin% c("rn",  i))]

  t[, comp := sum(.SD %>% unlist), .SDcols = v, by = "rn"]

  tt <- t[, .SD, .SDcols = c("comp", i)] %>%
          as.data.frame

  rownames(tt) <- t[, rn]

  ft <-  fisher.test(tt, simulate.p.value = TRUE, B =  2000)

  ci <-  paste(ft$conf.int %>% signif(digits = 3), collapse = "-")

  return(data.table::data.table(var = i, `Odds Ratio` = ft$estimate %>% signif(digits = 3),
                                `p-value` = ft$p.value, conf.interval = ci))

})  %>% data.table::rbindlist

pdt[, `p-value` := `p-value` %>% signif(digits = 3)]

pdt[, `Odds Ratio` := paste(`Odds Ratio`, " (", conf.interval, ")", sep = "")]

pdt %<>% as.data.frame

rownames(pdt) <- pdt["var"] %>% unlist

pdt <- pdt[, c("Odds Ratio", "p-value")]  %>% t

pdt[, 1] <- ""

tab <- rbind(tab, pdt)

tab  %<>% data.table::as.data.table(keep.rownames = TRUE)
tab %>% data.table::setnames("rn", "Immunogenicity")

tab[, names(tab)[2:3] := lapply(.SD, as.character), .SDcols = names(tab)[2:3]]

tab[Immunogenicity == "Positive", Immunogenicity := "Immunogenic"]
tab[Immunogenicity == "Negative", Immunogenicity := "Non-immunogenic"]

if((tab[4, 3]  %>% as.numeric) < 2.2e-16) tab[4, 3] <- "< 2.20e-16"

tt <- ttheme_default(padding = unit(c(10, 4), "mm"),
                    rowhead=list(fg_params=list(face = "bold")),
                    colhead=list(fg_params=list(col="white"),
                                  bg_params = list(fill="navyblue", alpha = 0.5)
                                ),
                    core=list(

          fg_params=list(fontface=c("plain", "plain", "bold", "bold.italic")),
          bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                      length.out=4), "#6BAED6"),
                           alpha = rep(c(1,0.5), each=5))

                         )
                         )

# pad out table to max characters
gtab_dis <- tableGrob(tab, rows = NULL, theme = tt)

gtab_dis$widths <-  unit(rep(1/ncol(gtab_dis), ncol(gtab_dis)), "npc")

# now IEDB alignment score
dt <- "Chowell_st_ie.txt" %>% data.table::fread

# maximally permissive, try to maximize sensitivity
dt[iedb_score > 0, res := "IEDB score > 0"]
dt[iedb_score == 0, res := "IEDB score = 0"]

t <- dt[,  .N, by = c("res", "Immunogenicity")] %>% dcast(Immunogenicity ~ res, value.var  = "N")

tab <- t[, .SD, .SDcols = 2:3] %>% as.data.frame
rownames(tab) <-  t[, Immunogenicity]

tab <- tab[, c("IEDB score = 0", "IEDB score > 0")]

# perform local tests by checking all columns against all others combined
pdt <- lapply(colnames(tab), function(i){

  t <- tab %>% data.table::as.data.table(keep.rownames =  TRUE)

  # sum all  other columns vs iterant
  v <- names(t)[which(!names(t) %chin% c("rn",  i))]

  t[, comp := sum(.SD %>% unlist), .SDcols = v, by = "rn"]

  tt <- t[, .SD, .SDcols = c("comp", i)] %>%
          as.data.frame

  rownames(tt) <- t[, rn]

  ft <-  fisher.test(tt, simulate.p.value = TRUE, B =  2000)

  ci <-  paste(ft$conf.int %>% signif(digits = 3), collapse = "-")

  return(data.table::data.table(var = i, `Odds Ratio` = ft$estimate %>% signif(digits = 3),
                                `p-value` = ft$p.value, conf.interval = ci))

})  %>% data.table::rbindlist

pdt[, `p-value` := `p-value` %>% signif(digits = 3)]

pdt[, `Odds Ratio` := paste(`Odds Ratio`, " (", conf.interval, ")", sep = "")]

pdt %<>% as.data.frame

rownames(pdt) <- pdt["var"] %>% unlist

pdt <- pdt[, c("Odds Ratio", "p-value")]  %>% t

pdt[, 1] <- ""

tab <- rbind(tab, pdt)

tab  %<>% data.table::as.data.table(keep.rownames = TRUE)
tab %>% data.table::setnames("rn", "Immunogenicity")

tab[, names(tab)[2:3] := lapply(.SD, as.character), .SDcols = names(tab)[2:3]]

tab[Immunogenicity == "Positive", Immunogenicity := "Immunogenic"]
tab[Immunogenicity == "Negative", Immunogenicity := "Non-immunogenic"]

if((tab[4, 3]  %>% as.numeric) < 2.2e-16) tab[4, 3] <- "< 2.20e-16"

tt <- ttheme_default(padding = unit(c(10, 4), "mm"),
                    rowhead=list(fg_params=list(face = "bold")),
                    colhead=list(fg_params=list(col="white"),
                                  bg_params = list(fill="navyblue", alpha = 0.5)
                                ),
                    core=list(

          fg_params=list(fontface=c("plain", "plain", "bold", "bold.italic")),
          bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                      length.out=4), "#6BAED6"),
                           alpha = rep(c(1,0.5), each=5))

                         )
                         )

# pad out table to max characters
gtab_ie <- tableGrob(tab, rows = NULL, theme = tt)

gtab_ie$widths <-  unit(rep(1/ncol(gtab_ie), ncol(gtab_ie)), "npc")

# now only peptides that included HLA data for predictions
dt <- "Chowell_st_ie.txt" %>% data.table::fread
b <- "Chowell_preds.txt" %>% data.table::fread

dt <- merge(dt, b, by = "nmer")

dt[Ensemble_score > 50, res := "Affinity > 50nM"]
dt[Ensemble_score < 50, res := "Affinity < 50nM"]

t <- dt[,  .N, by = c("res", "Immunogenicity")] %>% dcast(Immunogenicity ~ res, value.var  = "N")

tab <- t[, .SD, .SDcols = 2:3] %>% as.data.frame
rownames(tab) <-  t[, Immunogenicity]

tab <- tab[, c("Affinity > 50nM", "Affinity < 50nM")]

# perform local tests by checking all columns against all others combined
pdt <- lapply(colnames(tab), function(i){

  t <- tab %>% data.table::as.data.table(keep.rownames =  TRUE)

  # sum all  other columns vs iterant
  v <- names(t)[which(!names(t) %chin% c("rn",  i))]

  t[, comp := sum(.SD %>% unlist), .SDcols = v, by = "rn"]

  tt <- t[, .SD, .SDcols = c("comp", i)] %>%
          as.data.frame

  rownames(tt) <- t[, rn]

  ft <-  fisher.test(tt, simulate.p.value = TRUE, B =  2000)

  ci <-  paste(ft$conf.int %>% signif(digits = 3), collapse = "-")

  return(data.table::data.table(var = i, `Odds Ratio` = ft$estimate %>% signif(digits = 3),
                                `p-value` = ft$p.value, conf.interval = ci))

})  %>% data.table::rbindlist

pdt[, `p-value` := `p-value` %>% p.adjust(method = "bonferroni")]
pdt[, `p-value` := `p-value` %>% signif(digits = 3)]

pdt[, `Odds Ratio` := paste(`Odds Ratio`, " (", conf.interval, ")", sep = "")]

pdt %<>% as.data.frame

rownames(pdt) <- pdt["var"] %>% unlist

pdt <- pdt[, c("Odds Ratio", "p-value")]  %>% t

pdt[, 1] <- ""

tab <- rbind(tab, pdt)

tab  %<>% data.table::as.data.table(keep.rownames = TRUE)
tab %>% data.table::setnames("rn", "Immunogenicity")

tab[, names(tab)[2:3] := lapply(.SD, as.character), .SDcols = names(tab)[2:3]]

tab[Immunogenicity == "Positive", Immunogenicity := "Immunogenic"]
tab[Immunogenicity == "Negative", Immunogenicity := "Non-immunogenic"]

if((tab[4, 3]  %>% as.numeric) < 2.2e-16) tab[4, 3] <- "< 2.20e-16"

tt <- gridExtra::ttheme_default(padding = grid::unit(c(10, 4), "mm"),
                    rowhead=list(fg_params=list(face = "bold")),
                    colhead=list(fg_params=list(col="white"),
                                  bg_params = list(fill="navyblue", alpha = 0.5)
                                ),
                    core=list(

          fg_params=list(fontface=c("plain", "plain", "bold", "bold.italic")),
          bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                      length.out=4), "#6BAED6"),
                           alpha = rep(c(1,0.5), each=5))

                         )
                         )

# pad out table to max characters
gtab3 <- tableGrob(tab, rows = NULL, theme = tt)

gtab3$widths <-  unit(rep(1/ncol(gtab3), ncol(gtab3)), "npc")

png(filename = "Fig2D.png", type = "quartz", width = 4.5, height = 5, units = "in", res = 300)
grid.newpage()
grid.arrange(gtab_dis, gtab_ie, gtab3, nrow = 3)
dev.off()
