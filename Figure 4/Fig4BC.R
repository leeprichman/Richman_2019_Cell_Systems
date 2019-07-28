library(antigen.garnish)
library(ggplot2)
library(magrittr)
library(data.table)
library(survival)
library(pheatmap)

# heatmaps with melanoma data included
# function to return Cox Proportional Hazard model HR, confidence interval, and p-value based on a single column at a given cutoff
iterateR <- function(dt, col, cut1 = 0, cut2 = 0){

  a <- dt[, .SD %>% unique, .SDcols = c("time", "event", "sample_id", col)]

  a %<>% stats::na.omit %>% unique

  a[, time := time %>% as.numeric]

  a %>% data.table::setnames(col, "col")

  a[, col := col %>% as.numeric]

  a[col > cut1, group := "high"]
  a[col <= cut2, group := "low"]

  a <- a[!is.na(group)]

  a[, group := as.factor(group)]

  fit <- survival::survfit(survival::Surv(time, event) ~ group,
                             data = a)

  hr <- paste("HR: ")
  pval <- paste("p = ")

  dt <- data.table::data.table(HR = 1, p = 1)

  if (!a[, group %>% levels %>% length] < 2){

    cx <- survival::coxph(formula = survival::Surv(time, event) ~ group,
                                                       data = a,
                                                       iter.max = 200)

    hr <- paste("HR:",
    summary(cx)$coefficients[2] %>%
    formatC(., format = "e", digits = 2))

    pval <- paste("p = ",
    summary(cx)$logtest[3] %>%
    signif(digits = 2))

    dt <- data.table::data.table(HR = summary(cx)$coefficients[2], p = summary(cx)$logtest[3],
    upper = summary(cx)$conf.int[1, 4], lower = summary(cx)$conf.int[1, 3])

  }

    return(dt)

}

# set our FDR
FDR  <- 0.05

# read in the main analysis data
dtl <- data.table::fread("../Main_analysis/combined_output.txt")

# merge with output from the all_dissimilarity_9mers
sans <- "all_dissimilar_neopeptides_patient_level.txt" %>% data.table::fread

dtl <- merge(dtl, sans, by = "sample_id")

n <- c("TMB", "CDNs", "ADNs", "IEDB high neoantigens",
      "high dissimilarity neoantigens", "all MHC binders",
     "all dissimilar neopeptides")

# take only columns we need
dtl <- dtl[, .SD %>% unique, .SDcols = c("sample_id", "event", "time", "source", n)]

# make any NA that slipped through  (none should have though) in TMB or classifications into zero
NA_to_0 <- function(v){

  v[which(is.na(v))] <- 0

  return(v)

}

dtl[, as.character(n) := lapply(.SD, NA_to_0), .SDcols = n]

# set median survival for later normalization per dataset
dtl[, med_surv := median(time), by = "source"]

# normalize to values medians for combination
rescale <- function(v){

  if (median(v) == 0) return(v)

  return(v / median(v))

}

dtl[, as.character(n) := lapply(.SD, rescale), .SDcols = n, by = "source"]

dt2 <- data.table::copy(dtl)

# normalize to time to median survival per dataset for combination
dt2[, time := time / med_surv]

# combine Rizvi and Hellmann in this analysis, both PD-1, PFS, and NSCLC
# we did not combine melanoma sets because of CTLA-4 vs PD-1 and use of treatment resistant and naive patients in Riaz dataset
dt2[source %chin% c("Rizvi", "Hellmann"), source := "Combined\nNSCLC"]

dtl <- data.table::rbindlist(list(dtl, dt2[source == "Combined\nNSCLC"]), use.names = TRUE)

qu <- c(0.50)

names(qu) <- c("median")

# construct table iterating over all cutoffs for all comparisons for each source
mt <- lapply(dtl[, source %>% unique], function(s){

  s2 <- dtl[source == s]

  c <- lapply(qu, function(q){

    q1 <- 1 - q
    q2 <- 0 + q

    b <- lapply(n, function(i){

      a <- iterateR(s2, col = i,
        cut1 = quantile(s2[, get(i)], q1, na.rm = TRUE),
        cut2 = quantile(s2[, get(i)], q2, na.rm = TRUE))

      a[, variable := i]

      return(a)

    }) %>% data.table::rbindlist(use.names = TRUE, fill = TRUE)

    b[, cut := names(qu)[which(qu == q)]]

    return(b)

  }) %>% data.table::rbindlist(use.names = TRUE)

  c[, source := s]

  return(c)

}) %>% data.table::rbindlist(use.names = TRUE)

# adjust p-values to False Discovery Rate (Benjamini-Hochberg)
mt[, p := p %>% p.adjust(method = "BH"), by = c("cut", "source")]

# set names for aesthetics
mt[, variable := variable %>% stringr::str_replace("\\ neoantigens", "")]

# write Supplementary_Data_1, all comparisons table, HR, FDR, and CI for HR
mt %>% data.table::fwrite("Supplementary_Data_1.txt", sep = "\t", quote = TRUE, row.names = FALSE)

# show that no comparisons were significant for melanoma sets
# Riaz, Snyder and Van Allen do not return any comparisons that meet 0.05 family-wise error rate
mt[p < FDR, .N, by = "source"]

# look at only Rizvi, Hellmann, and combined (NSCLC, PD-1, PFS)
mt <- mt[source %chin% c("Rizvi", "Hellmann", "Combined\nNSCLC")]

mat <- data.table::copy(mt)

# for hazard ratio heatmaps, set insignificant comparisons to white tiles this way
mat[p > FDR, HR := as.numeric(NA)]

hrm <- mat[, .SD, .SDcols = c("HR", "variable", "cut", "source")]

# split by quartiles and deciles
hrm_c <- hrm %>% data.table::copy %>%
          split(by = "cut")

# function to save heatmap pdf
save_pheatmap_pdf <- function(hm, filename) {

  pdf(file = filename,
  width = 7, height = 7)

  grid::grid.newpage()

  grid::grid.draw(hm$gtable)

  dev.off()

}

# set column order
orderv1 <- c("Hellmann", "Rizvi",  "Combined\nNSCLC")


# set row order
orderc <- c("TMB", "all MHC binders", "CDNs", "ADNs", "IEDB high",
           "high dissimilarity", "all dissimilar neopeptides")

NA_to_blank <- function(v){

  v[which(is.na(v) | v == "NA")] <- ""

  return(v)

}

# generate HR heatmaps for quartiles and median
lapply(hrm_c %>% seq_along, function(i){

  dt <- hrm_c[[i]] %>% data.table::copy

  q <- dt[, cut %>% unique]

  fname <- paste("HR_", q, sep = "")

  dt[, cut := NULL]

  dt %<>% dcast(source ~ variable, value.var = "HR")

  mat <- dt %>% as.matrix(rownames = "source")

  mat <- mat[orderv1, orderc]

  mat %<>% t

  # format labels for cells
  lmat <- dt[, lapply(.SD, function(i) sprintf(i, fmt = "%.2f")),
          .SDcols = 2:ncol(dt)] %>%
          .[, lapply(.SD, function(i) NA_to_blank(i)),
                  .SDcols = 1:ncol(.)] %>%
          .[, source := dt[, source]] %>%
          as.matrix(rownames = "source")

  lmat <- lmat[orderv1, orderc]

  lmat %<>% t

  if (!all(is.na(mat))){

    cls <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name =
       "Reds"))(100)[1:100]

    bks <- seq(mt[p < 0.05, HR %>% min(na.rm = TRUE) %>% floor],
              mt[p < 0.05, HR %>% max(na.rm = TRUE) %>% ceiling],
             length.out = length(cls) + 1)

    hm <- try(pheatmap::pheatmap(mat, main = paste("Hazard ratios,", q),
              color = cls,
              breaks = bks,
       cluster_rows = FALSE,
       cluster_cols = FALSE,
       na_col = "#FFFFFF",
       display_numbers = lmat,
       number_color = "black",
       fontsize_number = 16,
       fontsize_row = 18,
       fontsize_col = 24))

       # Fontsize for title
       hm$gtable$grobs[[1]]$gp$fontsize <- 21
       # rotation of x axis labels
       hm$gtable$grobs[[3]]$rot <- 65
       hm$gtable$grobs[[3]]$hjust <- 1
       hm$gtable$grobs[[3]]$vjust <- 0.8
       # legend font size
       hm$gtable$grobs[[5]]$children[[2]]$gp$fontsize <- 14

       # get cell fill values, break into one dimension as expected for labels
       cv <- hm$gtable$grobs[[2]]$children[[1]]$gp$fill  %>% as.character

       # function to return binary black or white for labels based on fill of cell
       chromatch <- function(v, cls, na_col){

         dt <- data.table::data.table(orig = v %>% unique)

         dt[orig != na_col, log := which(cls == orig) > 50, by = orig]

         dt[orig != na_col, new := ifelse(log, "#FFFFFF", "#000000")]

         dt[orig == na_col, new := na_col]

         for (i in 1:length(v)){

           v[i] <- dt[orig == v[i], new]

         }

         return(v)

       }


       cv <- chromatch(cv, cls, "#FFFFFF")

       # numbers in cells
       hm$gtable$grobs[[2]]$children[[2]]$gp$col <- cv


    if (class(hm)[1] != "try-error")
       save_pheatmap_pdf(hm, "Fig4B.pdf")

  }

})

# now set up for p-value heatmaps
mat <- data.table::copy(mt)

mat[p > FDR, p := as.numeric(NA)]

pm <- mat[, .SD, .SDcols = c("p", "variable", "cut", "source")]

pm_c <- pm %>% data.table::copy %>%
          split(by = "cut")

lapply(pm_c %>% seq_along, function(i){

  dt <- pm_c[[i]] %>% data.table::copy

  q <- dt[, cut %>% unique]

  fname <- paste("p_", q, sep = "")

  dt[, cut := NULL]

  dt %<>% dcast(source ~ variable, value.var = "p")

  mat <- dt %>% as.matrix(rownames = "source")

  mat <- mat[orderv1, orderc]

  mat %<>% t

  # format labels for cells
  lmat <- dt[, lapply(.SD, function(i) sprintf(i, fmt = "%.3f")),
          .SDcols = 2:ncol(dt)] %>%
          .[, lapply(.SD, function(i) NA_to_blank(i)),
                  .SDcols = 1:ncol(.)] %>%
          .[, source := dt[, source]] %>%
          as.matrix(rownames = "source")

  lmat <- lmat[orderv1, orderc]

  lmat %<>% t

  if (!all(is.na(mat))){

  cls <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name =
       "Blues")))(100)[1:100]

  bks <- seq(mt[, p %>% min(na.rm = TRUE) %>% floor],
             FDR,
             length.out = length(cls) + 1)

  hm <- pheatmap::pheatmap(mat, main = paste("FDR,", q),
  color = cls,
  breaks = bks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "#FFFFFF",
  display_numbers = lmat,
  number_color = "white",
  fontsize_number = 16,
  fontsize_row = 18,
  fontsize_col = 24,
  fontsize  = 18)

  # Fontsize for title
  hm$gtable$grobs[[1]]$gp$fontsize <- 22
  # rotation of x axis labels
  hm$gtable$grobs[[3]]$rot <- 65
  hm$gtable$grobs[[3]]$hjust <- 1
  hm$gtable$grobs[[3]]$vjust <- 0.8
  # legend font size
  hm$gtable$grobs[[5]]$children[[2]]$gp$fontsize <- 14

  # get cell fill values, break into one dimension as expected for labels
  cv <- hm$gtable$grobs[[2]]$children[[1]]$gp$fill  %>% as.character

  # function to return binary black or white for labels based on fill of cell
  chromatch <- function(v, cls, na_col){

    dt <- data.table::data.table(orig = v %>% unique)

    dt[orig != na_col, log := which(cls == orig) < 50, by = orig]

    dt[orig != na_col, new := ifelse(log, "#FFFFFF", "#000000")]

    dt[orig == na_col, new := na_col]

    for (i in 1:length(v)){

      v[i] <- dt[orig == v[i], new]

    }

    return(v)

  }


  cv <- chromatch(cv, cls, "#FFFFFF")

  # numbers in cells
  hm$gtable$grobs[[2]]$children[[2]]$gp$col <- cv


  if (class(hm)[1] != "try-error")
    save_pheatmap_pdf(hm, "Fig4C.pdf")
  }

})

#  unused function to gradient scale colors for heatmap labels
# chromatch <- function(v, cls, na_col){
#
#   l <- grDevices::colorRampPalette(c("white", "black"))(length(cls))
#
#   dt <- data.table::data.table(orig = v %>% unique)
#
#   dt[orig != na_col, new := l[which(cls == orig)], by = orig]
#
#   dt[orig == na_col, new := na_col]
#
#   for (i in 1:length(v)){
#
#     v[i] <- dt[orig == v[i], new]
#
#   }
#
#   return(v)
#
# }
