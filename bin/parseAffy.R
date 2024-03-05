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
parser <- add_option(parser, opt_str = c("--out", "-o"),
                     type = "character",
                     help("File path and name to print output matrix"))
args <- parse_args(parser)

# Load cdf env. When cdf is NULL, ReadAffy searches for the annotation package
gpl <- if (args$cdf) make.cdf.env(args$cdf) else NULL

# Parse data
data <- ReadAffy(args$cel, cdfname = "gpl")

# Annotate data using cdf environment 
data <- rma(ReadAffy)

write.table(exprs(data), sep = "\t")