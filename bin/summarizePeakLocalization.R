#!/usr/bin/env Rscript

library <- function (...) { # It's a library! Shhh
   packages <- as.character(match.call(expand.dots = FALSE)[[2]])
   suppressWarnings(suppressMessages(lapply(packages, base::library, character.only = TRUE)))
   return(invisible())
}

library(optparse, GenomicFeatures, GenomicRanges, rtracklayer, magrittr)


option_list <- list(
  make_option(
    c('-p', '--peaks'), type='character', default=NULL, help='List of peak files to be analyzed; must be in CSV format', dest='peaks'
  ),
  make_option(
    c('-a', '--annotationFile'), type='character', default=NULL, help='GFF3 file with genome-specific annotations', dest='annotations'
  ),
  make_option(
    c('-o', '--outfile'), type='character', default=NULL, help='Prefix of output file for plots', dest='fh'
  ),
  make_option(
    c('-t', '--format'), type='character', default='svg', help='File format to save to (DEFAULT: svg)', dest='fout'
  )
)

opt_parser <- optparse::OptionParser(option_list=option_list)
opts <- optparse::parse_args(opt_parser)

peakFiles.list <- opts$peaks %>% 
  strsplit(., ",") %>%
  unlist %>%
  split(., ceiling(seq_along(.) / 5))

labels <- lapply(peakFiles.list, function(x) sub('\\.narrowPeak$', '', x))

tmpDim <- labels[[1]][1] %>%
  strsplit(., '_') %>%
  unlist %>%
  length()

labels %<>% lapply(
  function(x) strsplit(x, "_") %>%
    lapply(.,
           function(y) paste0(y[1:tmpDim-1], collapse=".")
           ) %>% unlist
)

testDim <- labels[[1]][1] %>%
  strsplit(., '/') %>%
  unlist %>%
  length

labels %<>% lapply(
  function(x) strsplit(x, '/') %>%
    lapply(., '[[', testDim) %>%
    unlist
)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
myGR.lol <- lapply(
  peakFiles.list, function(x) GRangesList(
    lapply(
    x,
    import,
    format = "BED",
    extraCols = extraCols_narrowPeak
    )
  )
)
for(i in 1:length(labels)){
  names(myGR.lol[[i]]) <- labels[[i]]
}

myTxDb <- suppressMessages(GenomicFeatures::makeTxDbFromGFF(opts$annotations))

genomicColors <- c(promoter = "#C40D00", geneDownstream = "#D55E00", geneBody = "#E66F00",
                   distalIntergenic = "#ffb38a", exon = "#000000", intron = "#666666", 
                   intergenic = "#DDDDDD", utr5 = "#0072B2", utr3 = "#56B4E9", 
                   CDS = "#0033BF", otherExon = "#009E73", undefined = "#FFFFFF")


lapply(myGR.lol, function(x) ChIPpeakAnno::genomicElementDistribution(
  peaks = x,  
  TxDb = myTxDb, 
  promoterRegion = c(upstream = 500, downstream = 200),
  labelColors = genomicColors, plot = F
  )$plot) -> myPlot.list

myPlot.list %<>% lapply(
  function(x) x + ggplot2::theme_classic(base_size = 16)+
    ggplot2::labs(x = 'Dataset', y = 'Percentage')
)

pdf(
  paste0(opts$fh, '.pdf'), width = 17.58, height = 7.94
)

for(p in myPlot.list) { print(p) }
dev.off()
