library(antigen.garnish)
library(magrittr)
library(data.table)
library(VennDiagram)

# ligthen color so  plot can be read
CDNcol <- "#ff6d00"
ADNcol <- "#2962ff"
IEcol <- "#00c853"
stcol <- "#DE7AA7"

lighten <- function(hex){

  hex2 <- colorRampPalette(c("#FFFFFF", hex))(10)[6]

  return(paste(hex2))

}

CDNcol %<>% lighten
ADNcol %<>% lighten
IEcol %<>% lighten
stcol %<>% lighten


# set up our classification colors
colpal <- c(CDNcol, ADNcol, IEcol, stcol, "#FFFFFF", "#BEBEBE")
names(colpal) <- c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens", "TMB", "MHC binders")

# read in combined output
dt <- data.table::fread("../Main_analysis/combined_output.txt")

# sum our metrics of interest for the venn diagram, the first ones are all antigen.garnish::garnish_summary output
n <- c("all MHC binders", "CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens",
     #  these are the overlap metrics produced by garnish_venn function from this analysis
     # see main analysis repo for this function
      "cd_ad", "cd_ie", "cd_st", "ad_ie", "ad_st", "ie_st",
      "cd_ad_ie", "cd_ad_st", "cd_ie_st", "ad_st_ie", "cd_ad_st_ie")

# sum all samples
v <- dt[, lapply(.SD, sum), .SDcols = n]

# get each group as a proportion of total MHC binders
v <- v[, n[2:length(n)] := lapply(.SD, function(i){

   (100 * i / `all MHC binders`) %>% signif(digits = 3)

 }), .SDcols = n[2:length(n)]]

# slot in our metrics from the garnish_venn function
area1 <- v[, `CDNs`]
area2 <- v[, `high dissimilarity neoantigens`]
area3 <- v[, ADNs]
area4 <- v[, `IEDB high neoantigens`]


n13 <- v[, cd_ad]
n14 <-  v[, cd_ie]
n12 <- v[, cd_st]
n34 <- v[, ad_ie]
n23 <- v[, ad_st]
n24 <- v[, ie_st]

n134 <- v[, cd_ad_ie]
n123 <- v[, cd_ad_st]
n124 <- v[, cd_ie_st]
n234 <- v[, ad_st_ie]
n1234 <- v[, cd_ad_st_ie]

# generate initial figure to edit
venn.plot <- VennDiagram::draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
    n34, n123, n124, n134, n234, n1234,
  category = c("CDNs", "high\ndissimilarity\nneoantigens", "ADNs", "IEDB high\nneoantigens"),
  lwd  = rep(1, 4), lty =  rep("dashed",  4),
  fill = colpal[c(1, 4, 2, 3)], fontface = rep("bold", 15), cex = rep(1.2, 15),
  alpha = rep (0.8,  4), cat.fontface = rep("bold", 4),
  cat.fontfamily = rep("sans", 4), fontfamily = rep("sans", 15),
  cat.cex = rep(1.5, 4), cat.dist = c(0.22, 0.25, 0.12, 0.12),
  margin = 0.03, sigdigs = 2, print.mode = c("raw"))

# get total proportion of classified neoantigens to make table later
pcts <- lapply(venn.plot %>% seq_along, function(i){

    return(venn.plot[[i]]$label)

}) %>% unlist %>%
        stringr::str_extract("[0-9,\\.]+") %>%
          na.omit %>% as.numeric

# change the labels to percent and round off
lapply(venn.plot %>% seq_along, function(i){

      l <- venn.plot[[i]]$label

      if (length(l) == 0) return(NULL)

      if (l %like% "^[0-9,\\.]+$")
      # weird behavior requires global env variable use
        venn.plot[[i]]$label <<- l %>% as.numeric %>%
                          signif(digits = 2) %>%
                          as.character %>% paste(., "%", sep = "")

          })

# save Venn diagram
png(filename = "Fig3F.png", type = "quartz", width = 8, height = 7, units = "in", res = 300)
grid::grid.newpage()
grid.draw(venn.plot)
dev.off()

# make accompanying table of total classifications as proportion of MHC binders
vt <- v %>% melt(variable.name = "Neoantigen class", value.name = "Percent of MHC binders")

# number of unclassified AGs
ol <- 100 - sum(pcts)

vt[, `Neoantigen class` := as.character(`Neoantigen class`)]

vt %<>% .[`Neoantigen class` %chin% c("CDNs", "ADNs", "IEDB high neoantigens", "high dissimilarity neoantigens")]

# add unclassified
vt <- data.table::rbindlist(list(vt, data.table(`Neoantigen class` = "unclassified",
                                                `Percent of MHC binders` = ol %>% signif(digits = 3))))
# convert it to a grob
library(gridExtra)

tab <- vt %>% gridExtra::tableGrob(rows = NULL)

# save PNG
png(filename = "Fig3F_tab.png", type = "quartz", width = 5, height = 2, units = "in", res = 300)
grid::grid.newpage()
grid.draw(tab)
dev.off()

# now lets check overlap of high dissimilarity neos with other classes at patient level
# can't get this from the garnish_venn output, need this from original tables
dt <- "../Figure 2/all_MHC_binders_w_id.txt" %>% data.table::fread

# first lets classify this list like we did in Figure 2
dt[, c("CDN", "ADN", "dissimilar", "IEDB") := 0]

dt[Ensemble_score < 50, CDN := 1]
dt[min_DAI > 10, ADN := 1]
dt[dissimilarity > 0.75, dissimilar := 1]
dt[iedb_score > 0.9, IEDB := 1]

dt[CDN == 1 & ADN == 0 & IEDB == 0 & dissimilar == 0, classification := "CDNs"]
dt[CDN == 0 & ADN == 1 & IEDB == 0 & dissimilar == 0, classification := "ADNs"]
dt[CDN == 0 & ADN == 0 & IEDB == 1 & dissimilar == 0, classification := "IEDB high"]
dt[CDN == 0 & ADN == 0 & IEDB == 0 & dissimilar == 1, classification := "high dissimilarity"]

# get our proportions of dissimilar neos that are doubly classified
dtl <- lapply(dt[, sample_id %>% unique], function(s){

  d <- dt[sample_id == s]

  ad <- d[dissimilar == 1 & ADN == 1] %>% nrow
  cd <- d[dissimilar == 1 & CDN == 1] %>% nrow
  ie <- d[dissimilar == 1 & IEDB == 1] %>% nrow

  denom <- d[dissimilar == 1] %>% nrow

  return(data.table::data.table(sample_id = s, `CDNs` = cd / denom, `ADNs` = ad / denom, `IEDB high` = ie /  denom))

}) %>% data.table::rbindlist(use.names = TRUE)

# drop the 73 out of 318 patients with no high dissimilarity neoantigens, can't use zero as a denominator
dtl %<>% na.omit

# long form for ggplot
gdt <- melt(dtl[, .SD, .SDcols = c("sample_id", "CDNs", "ADNs", "IEDB high")],
            id.vars = "sample_id")

# colors
CDNcol <- "#ff6d00"
ADNcol <- "#2962ff"
IEcol <- "#00c853"
stcol <- "#DE7AA7"

# set up our classification colors
colpal <- c(CDNcol, ADNcol, IEcol)
names(colpal) <- c("CDNs", "ADNs", "IEDB high")

plot_theme <- function(){

  pt <- ggplot2::theme_bw() + ggplot2::theme_classic() +
        ggplot2::theme(text = ggplot2::element_text(size = 32),
            plot.title = ggplot2::element_text(face = "bold",
                hjust = 0.5, size = ggplot2::rel(0.5)),
              #  axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7), color = "black", face = "bold"),
            axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.2),
                color = "black", face = "bold", angle = 45, hjust = 1),
            axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.2),
                color = "black", face = "bold"),
            axis.title = ggplot2::element_text(size = ggplot2::rel(1.2),
                face = "bold"),
                axis.line = element_line(size = 2),
                strip.text = element_text(size = ggplot2::rel(0.7), face = "bold"),
            legend.position = "none")

          return(pt)

}

library(ggplot2)

# extra plot of other neo types as proportion  of dissimilar
pl <- ggplot(gdt, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable), size = 2, col = "black", outlier.size = 2) +
  scale_fill_manual(values = colpal) +
  xlab(NULL) +
  ylab("Proportion of high\ndissimilarity neoantigens\n") +
  plot_theme()

ggsave(plot = pl, "../Figure S/high_dissimilarity_prop_boxplot.pdf", width = 9, height = 8.25)
