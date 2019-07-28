library(grid)
library(gridExtra)
library(antigen.garnish)
library(magrittr)
library(data.table)

# IEDB high and high dissimilarity neos using the cutoffs of  c(0.25, 0.9) as in other figures
# overall sensitivity though is quite low, so we use the permissive threshold of 1 and 0 instead in 2F
# as well, use optimal thresholds from ROC analysis c(1, 0.999941).
lapply(list(c(0.75, 0.9), c(0, 0.999941)), function(v){

  # make sure you have run Figure 2/Fig2C.R first, if you want to perform your own analysis
  # otherwise, use the data from our analysis provided in  the repository.
  dt <- "../Figure 2/Chowell_st_ie.txt" %>% data.table::fread

  # dissimilarity first
  dt[dissimilarity <= v[1], res := paste("Dissimilarity\n <", v[1] %>% substr(1, min(4, nchar(v[1]))))]
  dt[dissimilarity > v[1], res := paste("Dissimilarity\n >", v[1] %>% substr(1, min(4, nchar(v[1]))))]

  if (v[1] == 0) dt[, res := res %>% stringr::str_replace("<", "=")]

  t <- dt[,  .N, by = c("res", "Immunogenicity")] %>% dcast(Immunogenicity ~ res, value.var  = "N")

  tab <- t[, .SD, .SDcols = 2:3] %>% as.data.frame
  rownames(tab) <-  t[, Immunogenicity]

  tab <- tab[, c(names(tab)[names(tab) %like% "<|="], names(tab)[names(tab) %like% ">"])]

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

  tab %<>% data.table::as.data.table(keep.rownames = TRUE)
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
  dt <- "../Figure 2/Chowell_st_ie.txt" %>% data.table::fread

  dt[iedb_score <= v[2], res := paste("IEDB score <", v[2] %>% substr(1, min(4, nchar(v[2]))))]
  dt[iedb_score > v[2], res := paste("IEDB score >", v[2] %>% substr(1, min(4, nchar(v[2]))))]

  if (v[2] == 0) dt[, res := res %>% stringr::str_replace("<", "=")]

  t <- dt[,  .N, by = c("res", "Immunogenicity")] %>% dcast(Immunogenicity ~ res, value.var  = "N")

  tab <- t[, .SD, .SDcols = 2:3] %>% as.data.frame
  rownames(tab) <-  t[, Immunogenicity]

  tab <- tab[, c(names(tab)[names(tab) %like% "<|="], names(tab)[names(tab) %like% ">"])]

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

  fn <- paste("Contingency_table_dissimilarity_",
              v[1], "_iedb_", v[2], ".png", sep = "")

  png(filename = fn, type = "quartz", width = 4.5, height = 3.33, units = "in", res = 300)
  grid.newpage()
  grid.arrange(gtab_dis, gtab_ie, nrow = 2)
  dev.off()

})
