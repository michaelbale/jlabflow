#!/usr/bin/env Rscript

# Function to estimate RNA degredation from SE sequencing data
# Estimation is based on PEPPRO pipeline with all credit to Smith et al, and Sheffield (2021)
# 
# Function utilizes cutadapt report to calculate the ratio of small fragments to large fragments
# Print to STDOUT to redirect to a file. Can be parallelized at command-line

# TODO: make graphing function per sample-input

library <- function (...) { # It's a library! Shhh
   packages <- as.character(match.call(expand.dots = FALSE)[[2]])
   suppressWarnings(suppressMessages(lapply(packages, base::library, character.only = TRUE)))
   return(invisible())
}

calcRNADegrade <- function(report) {
	report %<>%
	 filter(expect / count < 0.01) %>%
	 mutate(invLength = max(length) - length) %>%
	 filter(invLength > 0)
	 
	if (20 %in% report$invLength) {
		degraded_upper <- 20
		degraded_lower <- 10
	} else {
		degraded_upper <- min(report$invLength) + 10
		degraded_lower <- max(1, degraded_upper - 10)
	}

	if (40 %in% report$invLength) {
		intact_upper <- 40
		intact_lower <- 30
	} else {
		intact_upper <- max(report$invLength)
		intact_lower <- max(1, intact_upper - 10)
	}

	degraded_count <- report %>%
		filter(invLength >= degraded_lower, invLength <= degraded_upper) %>%
		summarise(total = sum(count)) %>%
		pull(total)

	intact_count <- report %>%
		filter(invLength >= intact_lower, invLength <= intact_upper) %>%
		summarise(total = sum(count)) %>%
		pull(total)

	result <- degraded_count %>%
		divide_by(max(1, intact_count)) %>% 
		round(3) %>%
		format(nsmall = 3) %>%
		trimws
	return(result)
}

library(optparse, data.table, dplyr, magrittr)


option_list <- list(
  make_option(
    c('-r', '--report'), type='character', default=NULL, help='List of peak files to be analyzed; must be in CSV format', dest='report'
  )
)

opt_parser <- optparse::OptionParser(option_list=option_list)
opts <- optparse::parse_args(opt_parser)
myReports <- opts$report %>%
  strsplit(',') %>%
  unlist %>%
  set_names(
    strsplit(., '_') %>%
	lapply(
	  '[[',
	  1
	) %>%
	unlist
  ) %>%
  lapply(
    function(x) fread(x) %>% suppressWarnings
  )
  
myResults <- myReports %>%
  lapply(calcRNADegrade) %>%
  unlist %>%
  data.frame %>%
  set_colnames('mRNA Degradation')
  
write.table(myResults, '', quote = F, row.names = T, sep = '\t')




