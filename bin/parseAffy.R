#!/usr/bin/env Rscript
#' ---
#' Rscript definition of the parse method for the python AffymetrixParser class
#' ---

# Load libraries
library(optparse)
library(affy)
library(makecdfenv)

# Parse args
parser <- OptionParser()
parser <- add_option(parser, opt_str = c("--cel", "-i"),
                     type = "character", default = NULL,
                     help = "Path to affymetrix CEL file")
parser <- add_option(parser, opt_str = c("--cdf", "-c"),
                     type = "character", default = NULL,
                     help = "Path to affymetrix CDF file")
args <- parse_args(parser)

# Load cdf env
if (args$cdf) {
  gpl <- makecdfenv(args$cdf)
  data <- ReadAffy(args$cel, cdf_name = gpl)
} else {
  data <- ReadAffy(args$cel)
}

write.table(exprs(data), sep = "\t")